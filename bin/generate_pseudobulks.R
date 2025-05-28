
library(Seurat)
library(reticulate)
use_condaenv("/home/rschwartz/anaconda3/envs/r4.3/")
library(sceasy)
library(argparse)
library(tidyr)
options(future.globals.maxSize = 5 * 1024^3)  # 5 GB


parser = argparse::ArgumentParser(description = "Convert H5AD to H5Seurat.")
parser$add_argument("--rds_file", type="character", help="Paths to rds file", nargs='+')

args = parser$parse_args()
rds_files = args$rds_file
