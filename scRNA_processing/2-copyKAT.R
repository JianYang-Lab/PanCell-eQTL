#!/usr/bin/env Rscript

#' CopyKAT Analysis Pipeline for scRNA-seq Data
#'
#' This script processes single-cell RNA sequencing data by:

#' 1. Loading h5ad files processed with scanpy

#' 2. Converting to Seurat objects

#' 3. Running CopyKAT to detect copy number alterations
#'
#' Usage: Rscript copyKAT.R <sample_name>

# Load required libraries -----------------------------------------

suppressPackageStartupMessages({
  library(reticulate)
  library(Seurat)
  library(anndata)
  library(copykat)
})

# Setup Python environment ----------------------------------------
use_condaenv('scanpy_env')
sc <- import("scanpy")

# Parse command line arguments ------------------------------------
args <- commandArgs(TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript copyKAT.R <tissue> <sample_name>")

}
tissue <- args[1]
sample_name <- args[2]  # example: PRJNA782639_TD1; PRJNA634159_SRR11821880

# Define paths and settings ---------------------------------------
h5ad_dir <- paste0("/path/to/your/project/",tissue,"/h5ad/")

output_dir <- paste0("/path/to/your/project/",tissue,"/h5ad_merge/tumor_copykat_prediction")

h5ad_file <- file.path(h5ad_dir, paste0(sample_name, "_rawQC.h5ad"))

# Define reference cell types (immune cells) for CopyKAT
# Change this according to immune cells available in your data
reference_cell_types <- c(
  'Alveolar Mph CCL3+', 'Alveolar Mph MT-positive', 'Alveolar Mph proliferating', 
  'Alveolar macrophages', 'B cells', 'CD4 T cells', 'CD8 T cells', 
  'Classical monocytes', 'DC1', 'DC2', 'Migratory DCs', 'Mast cells',
  'Monocyte-derived Mph', 'NK cells', 'Non-classical monocytes', 
  'Plasma cells', 'Plasmacytoid DCs', 'T cells proliferating'
)

# Load and process data ------------------------------------------
message("Reading h5ad file:", h5ad_file)
h5ad <- sc$read_h5ad(h5ad_file)

# Create Seurat object

message("Converting to Seurat object...")
counts <- Matrix::Matrix(t(as.matrix(h5ad$X)), sparse = TRUE)
seurat_obj <- CreateSeuratObject(counts, meta.data = h5ad$obs)
rownames(seurat_obj) <- rownames(counts)
colnames(seurat_obj) <- colnames(counts)

# Clean up metadata
seurat_obj@meta.data <- seurat_obj@meta.data[, !colnames(seurat_obj@meta.data) %in% "barcode"]

# Get normal cells for reference
message("Identifying reference cells...")
normal_cell_names <- colnames(seurat_obj)[seurat_obj@meta.data[, "celltypist_pred"] %in% reference_cell_types]

if (length(normal_cell_names) == 0) {
  warning("No reference cells found matching the specified cell types!")

}

# Run CopyKAT ---------------------------------------------------
message("Running CopyKAT analysis...")
setwd(output_dir)

# Prepare count matrix
exp_matrix <- as.matrix(seurat_obj@assays$RNA$counts)

# Execute CopyKAT
copykat_results <- copykat(
  rawmat = exp_matrix, 
  id.type = "S", 
  ngene.chr = 5, 
  win.size = 25, 
  KS.cut = 0.1,
  sam.name = sample_name, 
  distance = "euclidean", 
  norm.cell.names = normal_cell_names,
  output.seg = FALSE, 
  plot.genes = TRUE, 
  genome = "hg20", 
  n.cores = 4
)

message("CopyKAT analysis completed for sample:", sample_name)

# pred.test <- data.frame(copykat_results$prediction)
# pred.test <- pred.test[which(pred.test$copykat.pred %in% c("aneuploid","diploid")), ]
# CNA.test <- data.frame(copykat_results$CNAmat)

print("Session completed")