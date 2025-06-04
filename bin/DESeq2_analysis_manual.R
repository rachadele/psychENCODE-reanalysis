
library(DESeq2)
library(tibble)
library(argparse)
library(tidyr)
library(tibble)
library(dplyr)
library(ggplot2)
library(edgeR)
options(future.globals.maxSize = 20 * 1024^3)  # 5 GB
set.seed(42)  # for reproducibility
# make base theme for ggplot with white background and larger text
base_theme <- theme(
  panel.background = element_rect(fill = "white", color = NA),
  plot.background = element_rect(fill = "white", color = NA),
  panel.border = element_rect(color = "black", fill = NA),
  axis.line = element_line(color = "black"),
  axis.ticks = element_line(color = "black"),
  legend.background = element_rect(fill = "white", color = NA),
  legend.key = element_rect(fill = "white", color = NA)
)

parser = argparse::ArgumentParser(description = "Run DESeq2 analysis on pseudobulk matrix")
parser$add_argument("--pseudobulk_matrix", type = "character", default="/space/grp/rschwartz/rschwartz/psychENCODE-reanalysis/results/pseudobulks/manual/astrocyte/astrocyte_pseudobulk_matrix.tsv.gz",
					help = "Path to the pseudobulk matrix tsv gzipped file.")

parser$add_argument("--gemma_metadata", type = "character",
					default="/space/grp/rschwartz/rschwartz/psychENCODE-reanalysis/gemma/metadata",
					help = "Path to the gemma metadata directory")

args = parser$parse_args()
pseudobulk_matrix_path <- args$pseudobulk_matrix
gemma_metadata <- args$gemma_metadata


# metadata processing ---------------------------------------------------

# combine the gemma metadata files into one metadata file
metadata_files <- list.files(gemma_metadata, full.names = TRUE, pattern = ".tsv")
metadata_list <- lapply(metadata_files, function(x) {
  df <- read.table(x, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  df$Cohort <- gsub("_sample_meta.tsv", "", basename(x))  # extract cohort name from file path
  df[] <- lapply(df, as.character)  # convert all columns to character

  return(df)
})
# Combine while filling missing columns with NA
metadata <- bind_rows(metadata_list)
# for age death, remove +/- and convert to numeric

# make all controls the same
metadata$Disorder <- as.character(metadata$Disorder)
metadata$Disorder[metadata$Disorder == "control"] <- "Control"  # make sure control is capitalized

metadata$Age_death <- as.character(metadata$Age_death)
metadata$Age_death[metadata$Age_death == "NaN"] <- NA  # replace NaN with NA
metadata$Age_death <- as.numeric(gsub("\\+", "", metadata$Age_death))

# make sample rownames
rownames(metadata) <- metadata$Individual_ID


# pseudobulk processing ----------------------------------------------------

# read the pseudobulk matrix
pseudobulk_matrix <- read.table(pseudobulk_matrix_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

rownames(pseudobulk_matrix) <- pseudobulk_matrix$X
pseudobulk_matrix <- pseudobulk_matrix[, -which(colnames(pseudobulk_matrix)=="feature_name")]  # remove the feature_name column

# remove leading "X" from column names
colnames(pseudobulk_matrix) <- gsub("^X", "", colnames(pseudobulk_matrix))

# remove empty samples
libsizes <- colSums(pseudobulk_matrix)
# filter out samples with lib size = 0
bad_samples <- names(libsizes[libsizes == 0])                 
# print message
if (length(bad_samples) > 0) {
  pseudobulk_matrix <- pseudobulk_matrix[, !colnames(pseudobulk_matrix) %in% bad_samples]
  message("Removed samples with 0 library sizes: ", paste(bad_samples, collapse = ", "))
}


#get matching samples
matching_samples <- intersect(rownames(metadata), colnames(pseudobulk_matrix))

# put them in the same order
filtered_metadata <- metadata[matching_samples, ]
pseudobulk_matrix <- pseudobulk_matrix[, matching_samples]

# fill NA in metadata with "Unknown" for categorical variables
# fill with median for numeric variables
filtered_metadata$RIN <- as.numeric(as.character(filtered_metadata$RIN))
filtered_metadata$RIN[is.na(filtered_metadata$RIN)] <- median(filtered_metadata$RIN, na.rm = TRUE)

# do this for PMI
filtered_metadata$PMI <- as.numeric(as.character(filtered_metadata$PMI))
filtered_metadata$PMI[is.na(filtered_metadata$PMI)] <- median(filtered_metadata$PMI, na.rm = TRUE)


# ----------pre filtering of  pseudobulk matrix----------------------


# Filtering: CPM normalization to filter out lowly expressed genes (<0.5 in <30% of samples) and donors with <50 cells detected. 
# After filtering, cell types with <16 samples also removed

# Gene filter
cpm_mat <- cpm(pseudobulk_matrix, log = FALSE)
keep_genes <- rowSums(cpm_mat > 0.5) >= 0.3 * ncol(pseudobulk_matrix)
pseudobulk_matrix <- pseudobulk_matrix[keep_genes, ]


# create DESeq2 object ----------------
dds <- DESeqDataSetFromMatrix(
  countData = as.matrix(pseudobulk_matrix),
  colData = filtered_metadata,
  # need to add average UMI at earlier step
  design = ~Disorder  + Age_death + PMI + Biological_Sex + X1000G_ancestry # genotype missing
  # add other covariates from manuscript
)

# rachel changes ---------------------------------------------------
dds$Disorder <- relevel(dds$Disorder, ref = "Control")

dds <- estimateSizeFactors(dds, type = "poscounts")
dds <- DESeq(dds)

se <- SummarizedExperiment(log2(counts(dds, normalized=TRUE) + 1),
                           colData=colData(dds))
# the call to DESeqTransform() is needed to
# trigger our plotPCA method.

png("cohort_pca.png", width = 800, height = 800)
plotPCA( DESeqTransform(se), intgroup = c("Cohort") )
dev.off()

# -------model diagnostics and quality control plots ----------------

png("dispersion_plot.png", width = 800, height = 600)
plotDispEsts(dds)
dev.off()

cooks <- assays(dds)[["cooks"]]
max_cooks <- apply(cooks, 1, max, na.rm = TRUE)
png("cooks_distance_histogram.png", width = 800, height = 600)
hist(max_cooks, breaks = 100, main = "Max Cook’s Distance per Gene")
dev.off()


# ----------- differential expression plots ----------------

# Get all relevant result names (e.g., Disorder_X_vs_Control)
res_names <- resultsNames(dds)

# Create empty list to store results
res_list <- list()
df_list <- list()

# Loop through result names
for (res_name in res_names) {
  outdir <- res_name %>% gsub(".tsv", "", .) 
  dir.create(outdir, recursive = TRUE)
  # Get results
  res <- results(dds, name = res_name)
  res_list[[res_name]] <- res

  # Save MA plot
  png(file.path(outdir, "ma_plot.png"), width = 800, height = 600)
  DESeq2::plotMA(res, main = paste(res_name))
  dev.off()

  # Convert to data frame
  df <- as.data.frame(res)
  df <- df[order(df$pvalue), ]
  df_list[[res_name]] <- df

  # P-value histogram
  p <- ggplot(df, aes(x = pvalue)) +
    geom_histogram(bins = 50, fill = "gray", color = "black") +
    labs(title = paste("P-value Distribution:", res_name),
         x = "P-value", y = "Frequency") +
    base_theme

  ggsave(file.path(outdir, "pvalue_dist.png"), plot = p, width = 8, height = 6)

  # sort by p-value and save DF
  df <- df %>%
    arrange(pvalue) %>% 
    mutate(gene = rownames(df)) %>%
    select(gene, everything())  # move gene names to the first column
  write.table(df, file = file.path(outdir, "results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}

