library(GenomicRanges)
library(argparse)
library(dplyr)
library(tidyr)
library(readr)

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

args <- parser$parse_args()
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

#------------- process author metadata -------------------
author_metadata <- read.table(args$author_metadata, header = TRUE, stringsAsFactors = FALSE, fill = NA)


#match sample names --------------------
matching_samples <- intersect(colnames(author_pseudobulk), metadata$Individual_ID)
only_in_author <- setdiff(colnames(author_pseudobulk), metadata$Individual_ID)
only_in_pavlab <- setdiff(metadata$Individual_ID, colnames(author_pseudobulk))

# write report about sample matching, fill missing values with NA
# union of both sets
all_samples <- union(colnames(author_pseudobulk), metadata$Individual_ID)
matching_samples_summary <- data.frame(Individual_ID = all_samples,
							   in_author = ifelse(all_samples %in% colnames(author_pseudobulk), "True", "False"),
							   in_pavlab = ifelse(all_samples %in% metadata$Individual_ID, "True", "False"))

# join metadata with matching_samples_summary
matching_samples_table <- matching_samples_summary %>%
  left_join(author_metadata, by = "Individual_ID") %>% 
  group_by(Cohort) %>%
  summarise(
	n_samples = n(),
	n_matching = sum(in_author == "True" & in_pavlab == "True"),
	n_only_in_author = sum(in_author == "True" & in_pavlab == "False"),
	n_only_in_pavlab = sum(in_author == "False" & in_pavlab == "True")
  ) %>% ungroup()



