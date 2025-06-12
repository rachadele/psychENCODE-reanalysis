library(GenomicRanges)
library(argparse)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(edgeR)


base_theme <- theme_bw() +
  theme(
	# increase base font size
	text = element_text(size = 20),
	# remove grid lines
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank()
	
  )

# Parse command line arguments
parser <- ArgumentParser(description = "Wrangle author pseudobulks")
parser$add_argument("--author_pseudobulk", type = "character", 
		default="/space/grp/rschwartz/rschwartz/psychENCODE-reanalysis/source_data/pseudobulks/pseudobulk_expr/Astro.expr.bed.gz", 
		help = "Author pseudobulk file in bed format")
parser$add_argument("--pavlab_pseudobulk", type = "character", 
		default="/space/grp/rschwartz/rschwartz/psychENCODE-reanalysis/results/ct_pseudobulks/manual/astrocyte/astrocyte_pseudobulk_matrix.tsv.gz", 
		help = "PabLab pseudobulk file in tsv format (either gemma or manually created)")
parser$add_argument("--gemma_metadata", type = "character", 
		default="/space/grp/rschwartz/rschwartz/psychENCODE-reanalysis/gemma/metadata", 
		help = "metadata for all experiments")
parser$add_argument("--author_metadata", type = "character", 
		default="/space/grp/rschwartz/rschwartz/psychENCODE-reanalysis/source_data/PEC2_sample_metadata.txt", 
		help = "Author metadata file in text format")

parser$add_argument("--author_cell_type", type = "character", 
		default="Astro", 
		help = "Cell type name in author pseudobulk file")
parser$add_argument("--pavlab_cell_type", type = "character",
		default="astrocyte", 
		help = "Cell type name in PabLab pseudobulk file")

args <- parser$parse_args()

#------------- get cell type names -------------------

author_cell_type <- args$author_cell_type
pavlab_cell_type <- args$pavlab_cell_type

#------------- read input files -------------------


# Read the input file
author_pseudobulk <- readr::read_tsv(args$author_pseudobulk, col_names=TRUE)
author_pseudobulk <- author_pseudobulk %>% select(-c("#chr","start","end","length","strand"))
rownames(author_pseudobulk) <- author_pseudobulk$gene  # set gene_id as rownames
author_pseudobulk <- author_pseudobulk %>% select(-gene)  # remove gene_id column
# Read the PabLab pseudobulk file
pavlab_pseudobulk <- readr::read_tsv(args$pavlab_pseudobulk, col_names=TRUE)
rownames(pavlab_pseudobulk) <- pavlab_pseudobulk$feature_name  # set gene_id as rownames
pavlab_pseudobulk <- pavlab_pseudobulk %>% select(-feature_name)  # remove feature_name column

#------------- process gemma metadata -------------------
metadata_files <- list.files(args$gemma_metadata, full.names = TRUE, pattern = ".tsv")
metadata_list <- lapply(metadata_files, function(x) {
  df <- read.table(x, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  df$Cohort <- gsub("_sample_meta.tsv", "", basename(x))  # extract cohort name from file path
  df[] <- lapply(df, as.character)  # convert all columns to character

  return(df)
})
# Combine while filling missing columns with NA
gemma_metadata <- bind_rows(metadata_list)



# --------------- map pavlab colnames to gemma metadata -------------------
pavlab_samples <- colnames(pavlab_pseudobulk) 
# find matching samples in gemma metadata
# create list mapping
matching_list <- list()
for (sample_id in pavlab_samples) {
	matching_list[[sample_id]] = gemma_metadata[gemma_metadata$sample_id == sample_id, ]$Individual_ID
}

colnames(pavlab_pseudobulk) <- sapply(pavlab_samples, function(x) {
	if (x %in% names(matching_list)) {
		return(matching_list[[x]])
	} else {
		return(NA)  # if no match found, return NA
	}
})


#------------- process author metadata -------------------
author_metadata <- read.table(args$author_metadata, header = TRUE, stringsAsFactors = FALSE, fill = NA)


#match sample names --------------------
matching_samples <- intersect(colnames(author_pseudobulk), colnames(pavlab_pseudobulk))
only_in_author <- setdiff(colnames(author_pseudobulk), colnames(pavlab_pseudobulk))
only_in_pavlab <- setdiff(colnames(pavlab_pseudobulk), colnames(author_pseudobulk))

# write report about sample matching, fill missing values with NA
# union of both sets
all_samples <- union(colnames(author_pseudobulk), colnames(pavlab_pseudobulk))
matching_samples_summary <- data.frame(Individual_ID = all_samples,
							   in_author = ifelse(all_samples %in% colnames(author_pseudobulk), "True", "False"),
							   in_pavlab = ifelse(all_samples %in% colnames(pavlab_pseudobulk), "True", "False"))

		

# join metadata with matching_samples_summary
matching_samples_table <- matching_samples_summary %>%
  left_join(author_metadata, by = "Individual_ID") 

filename = paste0("matching_samples_", author_cell_type, "_", pavlab_cell_type, ".tsv")

# write the matching samples table to a file
write.table(matching_samples_table, 
			file = filename,,
			sep = "\t",
			row.names = FALSE,
			quote = FALSE)
  
  
  count_table <- matching_samples_table %>% 
  group_by(Cohort) %>%
  summarise(
	n_samples = n(),
	n_matching = sum(in_author == "True" & in_pavlab == "True"),
	n_only_in_author = sum(in_author == "True" & in_pavlab == "False"),
	n_only_in_pavlab = sum(in_author == "False" & in_pavlab == "True")
  ) %>% ungroup()


  # make pivot table for plotting
  pivot_table <- count_table %>%
  pivot_longer(cols = c(n_only_in_author, n_only_in_pavlab, n_matching), 
			   names_to = "source", 
			   values_to = "count")


p <- ggplot(pivot_table, aes(x = Cohort, y = count, fill = source)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_text(aes(label = count), 
            position = position_dodge(width = 0.8), 
            hjust = -0.1, 
            size = 4) +
  coord_flip() +
  labs(
    title = paste("Sample Matching Summary for", author_cell_type, "and", pavlab_cell_type),
    x = "Cohort",
    y = "Sample Count",
    fill = "Source"
  ) +
  theme_bw(base_size = 20)

ggsave(
  filename = paste0("sample_matching_summary_", author_cell_type, "_", pavlab_cell_type, ".png"),
  plot = p,
  width = 10, height = 8
)




# ------------- compare pseudobulks -------------------
# take the log CPM values for pavlab pseudobulk
pavlab_pseudobulk_cpm <- cpm(pavlab_pseudobulk, log = TRUE, prior.count = 1)
rownames(pavlab_pseudobulk_cpm) <- rownames(pavlab_pseudobulk)

# compare total counts of author and pavlab pseudobulks
author_total_counts <- colSums(author_pseudobulk)
pavlab_total_counts <- colSums(pavlab_pseudobulk_cpm)

# createa data frame for comparison by sample
all_samples = union(names(author_total_counts), names(pavlab_total_counts))
comparison_df <- data.frame(
	  Individual_ID = all_samples,
	  Author_Counts = ifelse(all_samples %in% names(author_total_counts), author_total_counts[all_samples], NA),
	  Pavlab_Counts = ifelse(all_samples %in% names(pavlab_total_counts), pavlab_total_counts[all_samples], NA)
)	%>%
  left_join(author_metadata %>% select(Individual_ID, Cohort), 
            by = "Individual_ID") %>%
  replace_na(list(Author_Counts = 0, Pavlab_Counts = 0))

# write the comparison data frame to a file
write.table(comparison_df, 
			file = paste0("pseudobulk_comparison_", author_cell_type, "_", pavlab_cell_type, ".tsv"),
			sep = "\t",
			row.names = FALSE,
			quote = FALSE)


comparison_df <- comparison_df %>%
  mutate(
    avg = (Author_Counts + Pavlab_Counts) / 2,
    diff = Author_Counts - Pavlab_Counts
  )

comparison_plot <- ggplot(comparison_df, aes(x = avg, y = diff, color = Cohort)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
  geom_point(alpha = 0.7, size = 3) +
  labs(
    title = paste("Bland-Altman Plot for", author_cell_type, "vs", pavlab_cell_type),
    x = "Average log CPM",
    y = "Author - Pavlab log CPM",
    color = "Cohort"
  ) +
  theme_bw(base_size = 20) +
  coord_cartesian(ylim = quantile(comparison_df$diff, c(0.01, 0.99)))  # trim outliers

# save the comparison plot
ggsave(
  filename = paste0("pseudobulk_comparison_plot_", author_cell_type, "_", pavlab_cell_type, ".png"),
  plot = comparison_plot,
  width = 10, height = 15
)

#-------------- gene correlations-------------------
# calculate gene-wise correlations between author and pavlab pseudobulks
# put them in the same order



common_genes <- intersect(rownames(author_pseudobulk), rownames(pavlab_pseudobulk_cpm))
common_samples <- intersect(colnames(author_pseudobulk), colnames(pavlab_pseudobulk_cpm))

author_pseudobulk_common <- author_pseudobulk[common_genes, common_samples]
pavlab_pseudobulk_common <- pavlab_pseudobulk[common_genes, common_samples]


cor_results <- sapply(colnames(author_pseudobulk_common), function(sample) {
  cor(author_pseudobulk_common[, sample], pavlab_pseudobulk_common[, sample], use = "pairwise.complete.obs", method = "pearson")
})

# 2. Convert to matrix for heatmap plotting
cor_mat <- matrix(cor_results, ncol = 1)
rownames(cor_mat) <- names(cor_results)
colnames(cor_mat) <- "Correlation"