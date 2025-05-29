
library(Seurat)
library(reticulate)
use_condaenv("/home/rschwartz/anaconda3/envs/r4.3/")
library(sceasy)
library(argparse)
library(tidyr)
options(future.globals.maxSize = 10 * 1024^3)  # 5 GB


parser = argparse::ArgumentParser(description = "Convert H5AD to H5Seurat.")
parser$add_argument("--rds_files", type="character", help="Paths to rds files", nargs='+')

args = parser$parse_args()
rds_files = args$rds_files

all_seurat_objects <- list()


for (file in rds_files) {
	name <- sub("\\.rds$", "", basename(file))
	seurat_obj <- readRDS(file)
	seurat_obj@meta.data$cohort <- strsplit(name, "_")[[1]][1]

	# need to keep original library size somewhere in here


	all_seurat_objects[[name]] <- seurat_obj
}

# Combine all Seurat objects into one
combined_seurat <- Reduce(function(x, y) merge(x, y, merge.data = TRUE), all_seurat_objects)
# Save the combined Seurat object
output_file <- "combined_seurat.rds"
saveRDS(combined_seurat, file = output_file)
cat("Combined Seurat object saved to", output_file, "\n")

