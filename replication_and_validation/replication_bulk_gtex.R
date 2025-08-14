################################################################################
#
# --- Match and Standardize sc-eQTL and Bulk GTEx eQTLs ---
#
# Description:
# This script identifies significant single-cell eQTLs (sc-eQTLs) and matches
# them to corresponding bulk tissue eQTLs from the GTEx project. It handles
# different allele encodings (ref/alt flipping) and calculates standardized
# effect sizes for both datasets to enable direct comparison.
#
# The script takes a specific tissue, its corresponding GTEx dataset name, and
# an ancestry as inputs.
#
################################################################################


## 1. SETUP
# ------------------------------------------------------------------------------
print("--- Section 1: Setup ---")

# Load required libraries
library(tidyverse)
library(data.table)

# --- User-Defined Parameters ---

# Define base paths for input and output
# Users should modify these paths to match their environment.
base_dir <- "/path/to/"
bulk_path <- "/path/to/gtex_resources_besd/eQTL_raw/"
output_dir <- file.path(base_dir, "replication/")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Get arguments from the command line
args <- commandArgs(TRUE)
tissue <- args[1]
bulk_name <- args[2]
ancestry <- args[3]
# --- Example for interactive use ---
# tissue <- "Lung"
# bulk_name <- "Lung"
# ancestry <- "EUR"

if (is.na(tissue) || is.na(bulk_name) || is.na(ancestry)) {
    stop("Usage: Rscript script.R <tissue_name> <gtex_bulk_name> <ancestry>")
}


## 2. HELPER FUNCTIONS
# ------------------------------------------------------------------------------
print("--- Section 2: Defining Helper Functions ---")

#' Calculate standardized beta and standard error from Z-score.
#'
#' @param z Z-score (slope / slope_se).
#' @param p Minor allele frequency.
#' @param n Sample size.
#' @return A data frame with standardized beta (`std_b_hat`) and SE (`std_se`).
calcu_std_b_se <- function(z, p, n) {
    std_b_hat <- z / sqrt(2 * p * (1 - p) * (n + z^2))
    std_se <- 1 / sqrt(2 * p * (1 - p) * (n + z^2))
    res <- data.frame(std_b_hat, std_se)
    return(res)
}

#' Load and preprocess GTEx bulk eQTL data.
#'
#' @param bulk_path Path to the directory containing GTEx files.
#' @param bulk_name The name of the GTEx tissue dataset.
#' @return A data.table with the loaded and formatted bulk eQTL data.
load_bulk <- function(bulk_path, bulk_name) {
    print(paste0("Loading ", bulk_name, " GTEx bulk data..."))
    file_to_load <- file.path(bulk_path, paste0("GTEx_Analysis_v8_QTLs_GTEx_Analysis_v8_eQTL_all_associations_", bulk_name, ".allpairs.txt.gz"))
    bulk_eqtl <- fread(file_to_load)
    
    # Pre-process variant and gene IDs
    bulk_eqtl[, c("Chr", "BP", "A2", "A1") := tstrsplit(variant_id, "_", keep = 1:4)]
    bulk_eqtl <- bulk_eqtl[Chr != "chrX"]
    bulk_eqtl[, phenotype_id := tstrsplit(gene_id, "[.]", keep = 1)]
    return(bulk_eqtl)
}

#' Match sc-eQTLs to bulk eQTLs and standardize effect sizes.
#'
#' @param sceqtl_path Path to the significant sc-eQTL results.
#' @param tissue The name of the tissue being analyzed.
#' @param ancestry The ancestry group being analyzed.
#' @param celltype The specific cell type to match.
#' @param bulk_eqtl The pre-loaded GTEx data.table.
#' @return A data.table of matched pairs with standardized betas.
match_df_bulk <- function(sceqtl_path, tissue, ancestry, celltype, bulk_eqtl) {
    print(paste0("Matching SNP-gene pairs for GTEx bulk and ", celltype, " sc-eQTL..."))
    
    # Load significant sc-eQTLs for the specific cell type
    df <- fread(paste0(sceqtl_path, celltype, "_", ancestry, "_pc5_sig.tsv.gz"), select = c(1, 2, 4, 6:9))
    
    # Pre-process sc-eQTL data
    df[, c("Chr", "BP", "A2", "A1") := tstrsplit(variant_id, "_", keep = 1:4)]
    df[, Chr := paste0("chr", Chr)]
    df[, maf := ifelse(af > 0.5, 1 - af, af)]

    # Match pairs based on gene, chromosome, position, and alleles
    # Case 1: Exact allele match (A1=A1, A2=A2)
    exact_match <- merge(df, bulk_eqtl,
        by = c("phenotype_id", "Chr", "BP", "A1", "A2"),
        all.x = FALSE, suffixes = c("_df", "_bulk")
    )[, .(phenotype_id, Chr, BP, A1, A2, maf_df, pval_nominal_df, slope_df, slope_se_df, ma_count_df, maf_bulk, pval_nominal_bulk, slope_bulk, slope_se_bulk, ma_count_bulk)]

    # Case 2: Flipped allele match (A1=A2, A2=A1)
    flip_match <- merge(df, bulk_eqtl,
        by.x = c("phenotype_id", "Chr", "BP", "A1", "A2"),
        by.y = c("phenotype_id", "Chr", "BP", "A2", "A1"),
        all.x = FALSE, suffixes = c("_df", "_bulk")
    )[, .(phenotype_id, Chr, BP, A1, A2, maf_df, pval_nominal_df, slope_df, slope_se_df, ma_count_df, maf_bulk, pval_nominal_bulk, slope_bulk, slope_se_bulk, ma_count_bulk)]
    
    # Correct the effect size sign for flipped matches
    flip_match[, slope_bulk := -slope_bulk]

    # Combine exact and flipped matches
    df_comb <- rbind(exact_match, flip_match)

    # Standardize effect sizes for both datasets
    # Standardize sc-eQTL beta
    df_comb[, c("slope_std", "slope_std_se") := {
        z <- slope_df / slope_se_df
        p <- maf_df
        n <- round(ma_count_df / maf_df / 2)
        res <- calcu_std_b_se(z, p, n)
        list(res$std_b_hat, res$std_se)
    }]
    # Standardize bulk eQTL beta
    df_comb[, c("beta_std", "beta_std_se") := {
        z <- slope_bulk / slope_se_bulk
        p <- maf_bulk
        n <- round(ma_count_bulk / maf_bulk / 2)
        res <- calcu_std_b_se(z, p, n)
        list(res$std_b_hat, res$std_se)
    }]

    # Add metadata columns
    df_comb[, tissue := tissue]
    df_comb[, celltype := celltype]

    return(df_comb)
}


## 3. MAIN ANALYSIS
# ------------------------------------------------------------------------------
print("--- Section 3: Main Analysis ---")

# --- Step 3.1: Load GTEx Bulk Data ---
bulk_eqtl <- load_bulk(bulk_path, bulk_name)

# --- Step 3.2: Identify sc-eQTL cell types to process ---
sceqtl_path <- file.path(base_dir, tissue, "sc-eQTL/results/sig/")
summary_file <- file.path(sceqtl_path, "eQTL_summary_5en8.tsv")
celltypes_to_process <- fread(summary_file)[ancestry == !!ancestry]$cell_type

print(paste0("Found ", length(celltypes_to_process), " cell types to process for ", tissue, " in ", ancestry, " ancestry."))

# --- Step 3.3: Match sc-eQTLs with Bulk Data for each cell type ---
df_comb_list <- lapply(celltypes_to_process, function(ct) {
    match_df_bulk(sceqtl_path, tissue, ancestry, ct, bulk_eqtl)
})
matched_sceqtl_bulk_dt <- rbindlist(df_comb_list)

# The original script saved the result in a list, so we replicate that structure
matched_sceqtl_bulk <- list()
matched_sceqtl_bulk[[tissue]] <- matched_sceqtl_bulk_dt

# --- Step 3.4: Save the final matched data ---
output_file <- file.path(output_dir, paste0("matched_sceqtl_bulk_", tissue, "_", ancestry, ".rds"))
saveRDS(matched_sceqtl_bulk, output_file)

print(paste("✅ Analysis complete. Matched data saved to:", output_file))