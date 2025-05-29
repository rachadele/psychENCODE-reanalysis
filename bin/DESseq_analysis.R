
library(DESeq2)
library(tibble)
library(argparse)
library(tidyr)
library(tibble)
library(dplyr)
options(future.globals.maxSize = 20 * 1024^3)  # 5 GB


parser = argparse::ArgumentParser(description = "Run DESeq2 analysis on pseudobulk matrix")
parser$add_argument("--pseudobulk_matrix", type = "character", default="/space/grp/rschwartz/rschwartz/psychENCODE-reanalysis/results/pseudobulks/astrocyte/astrocyte_pseudobulk_matrix.tsv.gz",
					help = "Path to the pseudobulk matrix tsv gzipped file.")

parser$add_argument("--gemma_metadata", type = "character",
					default="/space/grp/rschwartz/rschwartz/psychENCODE-reanalysis/gemma/metadata",
					help = "Path to the gemma metadata directory")

args = parser$parse_args()
pseudobulk_matrix <- args$pseudobulk_matrix
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

metadata$Age_death <- as.character(metadata$Age_death)
metadata$Age_death[metadata$Age_death == "NaN"] <- NA  # replace NaN with NA
metadata$Age_death <- as.numeric(gsub("\\+", "", metadata$Age_death))

# make sample rownames
rownames(metadata) <- metadata$Individual_ID


# pseudobulk processing ----------------------------------------------------

# read the pseudobulk matrix
pseudobulk_matrix <- read.table(pseudobulk_matrix, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

rownames(pseudobulk_matrix) <- pseudobulk_matrix$feature_name
pseudobulk_matrix <- pseudobulk_matrix[, -which(colnames(pseudobulk_matrix)=="feature_name")]  # remove the feature_name column

# remove leading "X" from column names
colnames(pseudobulk_matrix) <- gsub("^X", "", colnames(pseudobulk_matrix))

#get matching samples
matching_samples <- intersect(rownames(metadata), colnames(pseudobulk_matrix))

# put them in the same order
filtered_metadata <- metadata[matching_samples, ]
pseudobulk_matrix <- pseudobulk_matrix[, matching_samples]

# find NAs in pseudobulk matrix


# create DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = as.matrix(pseudobulk_matrix),
  colData = filtered_metadata,
  design = ~Disorder # add other covariates from manuscript
)

# Laureen's work goes here ------------------------------
# make sure to deal with missing information