
library(Seurat)
library(reticulate)
use_condaenv("/home/rschwartz/anaconda3/envs/r4.3/")
library(sceasy)
library(argparse)
library(tidyr)
options(future.globals.maxSize = 5 * 1024^3)  # 5 GB

rename_features <- function(seurat_obj) {
    counts <- seurat_obj[["RNA"]]@counts
    data <- seurat_obj[["RNA"]]@data 

    feature_meta <- seurat_obj[["RNA"]][[]]
    feature_meta$feature_id <- rownames(feature_meta)
    rownames(feature_meta) <- feature_meta$feature_name

    rownames(counts) <- rownames(feature_meta)
    rownames(data) <- rownames(feature_meta)

    newRNA <- CreateAssayObject(counts=counts)
    newRNA$data <- data
    newRNA[[]] <- feature_meta

   seurat_obj[["RNA"]] <- newRNA
   DefaultAssay(seurat_obj) <- "RNA"
   return(seurat_obj)
}


parser = argparse::ArgumentParser(description = "Convert H5AD to H5Seurat.")
parser$add_argument("--h5ad_file", type="character", help="Path to H5AD file.", 
        default = "/space/grp/rschwartz/rschwartz/get_gemma_data.nf/null_author_false_sample_split_false/homo_sapiens/DevBrain.h5ad")

args = parser$parse_args()
h5ad_file = args$h5ad_file


ad <- import("anndata", convert = FALSE)
adata <- ad$read_h5ad(h5ad_file)
counts <- py_to_r(adata$X)
counts <- counts %>% as.matrix() %>% t()
rownames(counts) <- py_to_r(adata$var$feature_name$to_list())
colnames(counts) <- py_to_r(adata$obs_names$to_list())


# remove NaN rownames
#counts <- counts[-grep("-", rownames(counts)), ]
valid_rows <- !is.na(rownames(counts)) & rownames(counts) != "NaN"
counts <- counts[valid_rows, ]
seurat_obj <- CreateSeuratObject(counts = counts)

# Get var and obs as data frames
var <- py_to_r(adata$var) %>% as.data.frame()
obs <- py_to_r(adata$obs) %>% as.data.frame()
# Set rownames of var and obs to match Seurat
var$feature_id <- rownames(var)
# Subset var to match valid gene rows
var <- var[valid_rows, ]
rownames(var) <- var$feature_name

seurat_obj@assays$RNA[[]] <- var
seurat_obj@meta.data <- obs

#if ("feature_id" %in% colnames(seurat_obj@assays$RNA[[]])) {
  #seurat_obj <- rename_features(seurat_obj, column_name="feature_id")
#}

saveRDS(seurat_obj, file = gsub(".h5ad",".rds", h5ad_file))
