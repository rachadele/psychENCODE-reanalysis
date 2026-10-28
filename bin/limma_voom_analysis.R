# limma-voom DE on a pseudobulk matrix, following Gemma's LinearModelAnalyzer for RNA-seq
# (DifferentialExpressionAnalysisConfig useWeights = true, moderateStatistics = true):
#   * log2 CPM = log2((count + 0.5) / (lib.size + 1) * 1e6), where lib.size is the sum over all
#     genes of the kept (unmasked) cells, computed before any gene filter;
#   * Gemma's DifferentialExpressionAnalysisFilter, in its order:
#       - samples with fewer than --min_cells_per_sample cells are dropped (Gemma default 100;
#         here 0 by default so the samples match the DESeq2 path),
#       - genes expressed in fewer than --min_cells_per_gene cells summed over the retained
#         samples are dropped (Gemma default 3; needs --gene_ncells),
#       - repetitive values: with >= 4 samples, genes with zero variance of the counts,
#       - genes with variance of log2 CPM below --min_variance (Gemma default 0.01);
#   * voom weights from the same design, lmFit with weights, eBayes (no trend, not robust);
#   * BH-adjusted p-values per coefficient.
# Same metadata handling, design and sample filters as bin/DESeq2_analysis.R (gemma mode).
# Output per coefficient: <coef>/<coef>_<cell_type>_results.tsv with limma columns (logFC,
# AveExpr, t, P.Value, adj.P.Val, B) plus DESeq2-named aliases (log2FoldChange = logFC,
# stat = t, pvalue = P.Value, padj = adj.P.Val, lfcSE = logFC / t, so lfc / lfcSE = moderated t)
# so the comparison steps read it unchanged.
suppressPackageStartupMessages({
  library(limma)
  library(edgeR)
  library(argparse)
  library(dplyr)
})
set.seed(42)

parser <- argparse::ArgumentParser(description = "limma-voom DE on a pseudobulk matrix (Gemma method)")
parser$add_argument("--pseudobulk_matrix", required = TRUE, help = "<ct>_pseudobulk_matrix.tsv.gz")
parser$add_argument("--metadata", default = "/space/grp/rschwartz/rschwartz/psychENCODE-reanalysis/gemma/metadata",
                    help = "directory of <cohort>_sample_meta.tsv files")
parser$add_argument("--cell_type", required = TRUE)
parser$add_argument("--gene_ncells", default = NULL,
                    help = "<ct>_gene_ncells.tsv.gz: expressing cells per gene x sample (same layout)")
parser$add_argument("--cell_counts", default = NULL, nargs = "+",
                    help = "per-study <study>_cell_counts.tsv files (sample, ct, n_kept, lib_size)")
parser$add_argument("--min_cells_per_sample", type = "integer", default = 0)
parser$add_argument("--min_cells_per_gene", type = "integer", default = 3)
parser$add_argument("--min_variance", type = "double", default = 0.01)
args <- parser$parse_args()
cell_type <- args$cell_type

# metadata (as in DESeq2_analysis.R, gemma mode) -------------------------------------------
metadata_files <- list.files(args$metadata, full.names = TRUE, pattern = ".tsv")
metadata <- bind_rows(lapply(metadata_files, function(x) {
  df <- read.table(x, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  df$Cohort <- gsub("_sample_meta.tsv", "", basename(x))
  df[] <- lapply(df, as.character)
  df
}))
metadata$Disorder[metadata$Disorder == "control"] <- "Control"
metadata$Age_death[metadata$Age_death == "NaN"] <- NA
metadata$Age_death <- as.numeric(gsub("\\+", "", metadata$Age_death))
rownames(metadata) <- metadata$Individual_ID

# counts --------------------------------------------------------------------------------------
read_mat <- function(path) {
  m <- read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  rownames(m) <- m$feature_name
  m <- m[, colnames(m) != "feature_name", drop = FALSE]
  colnames(m) <- gsub("^X", "", colnames(m))
  as.matrix(m)
}
counts <- read_mat(args$pseudobulk_matrix)

# library sizes: from the adapter's cell counts when given (exact), else column sums
lib_all <- colSums(counts)
n_cells <- NULL
if (!is.null(args$cell_counts)) {
  cc <- bind_rows(lapply(args$cell_counts, function(f) read.table(f, header = TRUE, sep = "\t",
                                                                   stringsAsFactors = FALSE,
                                                                   colClasses = c(sample = "character"))))
  cc <- cc[!is.na(cc$ct) & cc$ct == cell_type, ]
  if (anyDuplicated(cc$sample)) stop("duplicate sample names across studies for ", cell_type)
  rownames(cc) <- cc$sample
  # columns with no cells of this type (all-zero) have no cell_counts row: library size 0
  missing <- setdiff(colnames(counts), cc$sample)
  if (any(lib_all[missing] > 0)) stop("no cell counts for nonzero samples: ",
                                      paste(missing[lib_all[missing] > 0], collapse = ", "))
  lib_all <- setNames(ifelse(colnames(counts) %in% cc$sample, cc[colnames(counts), "lib_size"], 0), colnames(counts))
  n_cells <- setNames(ifelse(colnames(counts) %in% cc$sample, cc[colnames(counts), "n_kept"], 0), colnames(counts))
  rel <- max(abs(lib_all - colSums(counts)) / pmax(lib_all, 1))
  message("max relative difference, adapter library size vs matrix column sum: ", signif(rel, 3))
}

# sample filters (as in DESeq2_analysis.R) ----------------------------------------------------
bad <- names(lib_all)[lib_all == 0]
if (length(bad) > 0) message("Removed samples with 0 library sizes: ", paste(bad, collapse = ", "))
samples <- setdiff(colnames(counts), bad)
samples <- intersect(rownames(metadata), samples)
if (!is.null(n_cells) && args$min_cells_per_sample > 0) {
  few <- samples[n_cells[samples] < args$min_cells_per_sample]
  if (length(few) > 0) message("Removed ", length(few), " samples with < ", args$min_cells_per_sample,
                               " cells: ", paste(few, collapse = ", "))
  samples <- setdiff(samples, few)
}
message("Matching samples (", length(samples), "): ", paste(samples, collapse = ", "))
meta <- metadata[samples, ]
counts <- counts[, samples, drop = FALSE]
lib <- lib_all[samples]
meta$RIN <- as.numeric(meta$RIN)
meta$RIN[is.na(meta$RIN)] <- median(meta$RIN, na.rm = TRUE)
meta$PMI <- as.numeric(meta$PMI)
meta$PMI[is.na(meta$PMI)] <- median(meta$PMI, na.rm = TRUE)

# design (DESeq2 gemma-mode formula; factor levels as DESeq2 sets them) ------------------------
design_vars <- c("Disorder", "Age_death", "PMI", "Biological_Sex", "X1000G_ancestry")
if (any(is.na(meta[, design_vars]))) stop("missing values in design variables")
meta$Disorder <- relevel(factor(meta$Disorder), ref = "Control")
meta$Biological_Sex <- factor(meta$Biological_Sex)
meta$X1000G_ancestry <- factor(meta$X1000G_ancestry)
design <- model.matrix(~ Disorder + Age_death + PMI + Biological_Sex + X1000G_ancestry, data = meta)
# coefficient names as DESeq2 resultsNames() writes them
coef_names <- colnames(design)
deseq_names <- coef_names
deseq_names[coef_names == "(Intercept)"] <- "Intercept"
for (f in c("Disorder", "Biological_Sex", "X1000G_ancestry")) {
  ref <- levels(meta[[f]])[1]
  hit <- startsWith(coef_names, f) & !(coef_names %in% design_vars)
  deseq_names[hit] <- paste0(f, "_", make.names(substring(coef_names[hit], nchar(f) + 1)), "_vs_", make.names(ref))
}
keep_cols <- qr(design)$pivot[seq_len(qr(design)$rank)]
if (length(keep_cols) < ncol(design)) {
  message("design not full rank, dropping: ", paste(coef_names[-keep_cols], collapse = ", "))
  design <- design[, sort(keep_cols), drop = FALSE]
  deseq_names <- deseq_names[sort(keep_cols)]
}

# gene filters (Gemma's DifferentialExpressionAnalysisFilter) --------------------------------
n0 <- nrow(counts)
if (!is.null(args$gene_ncells) && args$min_cells_per_gene > 0) {
  nc <- read_mat(args$gene_ncells)
  nc <- nc[rownames(counts), samples, drop = FALSE]
  counts <- counts[rowSums(nc) >= args$min_cells_per_gene, , drop = FALSE]
}
n1 <- nrow(counts)
if (ncol(counts) >= 4) counts <- counts[apply(counts, 1, var) > 0, , drop = FALSE]
n2 <- nrow(counts)
logcpm <- t(log2(t(counts + 0.5) / (lib + 1) * 1e6))
counts <- counts[apply(logcpm, 1, var) >= args$min_variance, , drop = FALSE]
message(sprintf("genes: %d start, %d after min cells, %d after zero variance, %d after low variance",
                n0, n1, n2, nrow(counts)))

# voom + lmFit + eBayes --------------------------------------------------------------------------
png("voom_mean_variance.png", width = 800, height = 600)
v <- voom(counts, design, lib.size = lib, plot = TRUE)
dev.off()
fit <- eBayes(lmFit(v, design))
message("eBayes df.prior = ", signif(fit$df.prior, 4), ", s2.prior = ", signif(fit$s2.prior, 4))

for (i in seq_len(ncol(design))) {
  res_name <- deseq_names[i]
  tt <- topTable(fit, coef = i, number = Inf, sort.by = "P", adjust.method = "BH")
  df <- data.frame(gene = rownames(tt), tt, check.names = FALSE) %>%
    mutate(log2FoldChange = logFC, stat = t, pvalue = P.Value, padj = adj.P.Val,
           lfcSE = logFC / t, cell_type = cell_type)
  dir.create(res_name, showWarnings = FALSE, recursive = TRUE)
  write.table(df, file = file.path(res_name, paste0(res_name, "_", cell_type, "_results.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
}
