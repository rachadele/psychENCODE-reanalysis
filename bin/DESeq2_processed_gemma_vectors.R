library(Seurat)
library(dplyr)
library(ggplot2)
library(data.table)
library(reticulate)
use_condaenv("~/anaconda3/envs/r4.3/")
library(BiocParallel)
library(stringr)
library(argparse)
set.seed(42)
library(ggExtra)
library(patchwork)
library(DESeq2)

parser <- ArgumentParser(description = "DESeq2 analysis for PsychENCODE reanalysis manuscript")
parser$add_argument("--pseudobulk_matrix_dir", type="character",  help="Path to the pseudobulk matrix files from Gemma", default="/space/grp/rschwartz/rschwartz/psychENCODE-reanalysis/gemma/unfiltered_pseudobulks")
parser$add_argument("--metadata_files", type="character", help="Path to the metadata files from Gemma", default = "/space/grp/rschwartz/rschwartz/psychENCODE-reanalysis/gemma/metadata")

args <- parser$parse_args()
pseudobulk_matrix_dir <- args$pseudobulk_matrix_dir
metadata_files <- args$metadata_files

matrix_files <- list.files(pseudobulk_matrix_dir, full.names = TRUE)
#pseudobulks <- lapply(matrix_files, fread, header = TRUE, sep = "\t", data.table = TRUE)
pseudobulks <- lapply(matrix_files, function(x) {
  fread(x, header = TRUE, sep = "\t", data.table = TRUE)
})

names(pseudobulks) <- gsub(".txt", "", basename(matrix_files)) %>%
  gsub("_expmat.unfilt.data.gz", "", .)

