################################################################################
#
# --- Cross-Tissue eQTL Heterogeneity Analysis ---
#
# Description:
# This script performs a meta-analysis to assess the heterogeneity of eQTL 
# effects for a specific cell type across multiple tissues. For eQTLs found to 
# be significant in at least one tissue, it calculates Cochran's Q and the I^2 
# statistic to quantify the level of heterogeneity in effect sizes.
#
# Input:
# 1. Pre-computed eQTL summary statistics (.parquet files) for each tissue.
# 2. A list of significant eQTL pairs from each tissue.
# 3. A mapping file defining cell type names across tissues.
# 4. A list of overlapping variants across all datasets.
#
# Output:
# An RDS file containing the heterogeneity statistics (I^2, Q, QEp) for each 
# significant gene-variant pair.
#
################################################################################


## 1. SETUP
# ------------------------------------------------------------------------------
# Load required libraries
library(metafor)      # For meta-analysis (rma function)
library(tidyverse)    # For data manipulation (dplyr, etc.)
library(data.table)   # For efficient data handling (fread, dcast)
library(arrow)        # For reading parquet files
library(parallel)     # For parallel processing (mclapply)

# --- User-Defined Parameters ---

# Set the main project directory
# All other paths will be relative to this
base_dir <- "/path/to/your/scRNA_project_folder/" 
setwd(base_dir)

# Define the cell type for analysis. 
# This can be set manually or read from command line arguments.
# ct <- "T_cells" 
ct <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(ct)) {
  stop("Cell type not provided. Please specify a cell type to analyze.")
}

# Define the list of tissues to include in the analysis
tissue_list <- c("Blood", "Lung", "Skin", "Colon", "Liver")

# Set the number of cores for parallel processing
num_cores <- detectCores() - 2 # Using 2 fewer cores than available

# Define paths for input and output files
# Path to the directory for mashR-related intermediate files
mashr_dir <- "mashR/cross_tissue/"
# Path for final output
output_dir <- "cross_tissue_het/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)


## 2. HELPER FUNCTIONS
# ------------------------------------------------------------------------------

#' Convert p-value to Z-score
#' Assigns direction based on the sign of the effect size (beta).
ConvertP2Z <- function(pval, beta) {
    z <- abs(qnorm(pval / 2))
    z[beta < 0] <- -1 * z[beta < 0]
    return(z)
}

#' Handle missing or NaN effect sizes (beta/slope)
handle_nan_b <- function(x) {
    x[is.nan(x) | is.na(x)] <- 0
    return(x)
}

#' Handle missing, zero, or infinite standard errors
#' Replaces problematic values with a large number (1E3) to down-weight their influence.
handle_nan_s <- function(x) {
    x[is.nan(x) | is.infinite(x) | is.na(x) | x == 0] <- 1E3
    return(x)
}

#' Handle missing Z-scores based on beta and SE values
handle_nan_z <- function(z, b, s) {
    z[is.nan(s) | is.infinite(s) | is.na(s) | s == 0] <- 0
    z[is.nan(b) | is.na(b) | b == 0] <- 0
    z[is.nan(z) | is.na(z)] <- 0
    return(z)
}

#' Calculate Z-scores and clean summary statistics
GetSS <- function(dat) {
    dat$z <- ConvertP2Z(dat$pval_nominal, dat$slope)
    dat[['z']] <- handle_nan_z(dat[['z']], dat[['slope']], dat[['slope_se']])
    dat[['slope_se']] <- handle_nan_s(dat[['slope_se']])
    dat[['slope']] <- handle_nan_b(dat[['slope']])
    return(dat)
}

#' Process a single tissue-ancestry-chromosome combination
#' Reads and filters eQTL summary statistics from a parquet file.
process_combinations <- function(combo, geno_overlap, gene_keep) {
    tissue <- combo$tissue
    ct_raw <- combo$ct_raw
    anc <- combo$anc
    chrom <- combo$chrom
    
    file_path <- file.path(tissue, "sc-eQTL/results/raw", 
                           paste0(ct_raw, "_", anc, "_pc5.cis_qtl_pairs.", chrom, ".parquet"))
    
    if (file.exists(file_path)) {
        df <- read_parquet(file_path)
        df <- df[df$phenotype_id %in% gene_keep & df$variant_id %in% geno_overlap, ]
        df$ancestry <- anc
        df$tissue <- tissue
        return(df)
    }
    return(NULL) # Return NULL if file doesn't exist
}

#' Load and format summary statistics for a given chromosome
load_summary_statistics <- function(ct, chrom, geno_overlap, gene_keep, ct_cross_tissue) {
    print(paste0("Reading raw data for ", ct, " chr", chrom, "."))

    # Create a list of combinations to process
    combinations <- list()
    anc <- "EUR" # Analysis is focused on EUR ancestry
    for (ts in tissue_list) {
        ct_raw_name <- ct_cross_tissue %>% 
            dplyr::filter(tissue == ts & celltype == ct) %>% 
            pull(celltype_raw)
        
        if (length(ct_raw_name) > 0) {
            combinations[[length(combinations) + 1]] <- list(
                anc = anc,
                ct_raw = ct_raw_name,
                tissue = ts,
                chrom = chrom
            )
        }
    }

    # Process combinations in parallel
    df_list <- mclapply(combinations, process_combinations, 
                        geno_overlap = geno_overlap, gene_keep = gene_keep, mc.cores = num_cores)

    # Filter out NULL results and get cleaned summary stats
    df_list <- Filter(Negate(is.null), df_list)
    if (length(df_list) == 0) return(NULL)
    
    df_ss_list <- mclapply(df_list, GetSS, mc.cores = num_cores)
    df_all_ss <- rbindlist(df_ss_list)

    # Reshape data from long to wide format
    df_slope_wide <- dcast(df_all_ss, phenotype_id + variant_id ~ tissue, value.var = "slope")
    df_se_wide <- dcast(df_all_ss, phenotype_id + variant_id ~ tissue, value.var = "slope_se")

    return(list(df_all_slope_wide = df_slope_wide, df_all_se_wide = df_se_wide))
}


## 3. MAIN ANALYSIS
# ------------------------------------------------------------------------------

# --- A. Load cross-tissue cell type mapping and significant eQTLs ---
print("Step 1: Loading prerequisite data.")
ct_cross_tissue <- read_delim("cross_tissue_celltype.txt", delim = "\t")

sig_list <- list()
anc <- "EUR"
for (ts in tissue_list) {
    ct_raw_name <- ct_cross_tissue %>% 
        dplyr::filter(tissue == ts & celltype == ct) %>% 
        pull(celltype_raw)
    
    if (length(ct_raw_name) > 0) {
        sig_file <- paste0(ts, "/sc-eQTL/results/sig/", ct_raw_name, "_", anc, "_pc5_sig_5en8.tsv.gz")
        if (file.exists(sig_file)) {
            df <- fread(sig_file)
            df <- df[pval_nominal < 5e-8, .(phenotype_id, variant_id)]
            sig_list[[length(sig_list) + 1]] <- df
        }
    }
}
sig_pairs <- distinct(rbindlist(sig_list))
print(paste("Found", nrow(sig_pairs), "unique significant eQTL pairs across all tissues."))

# --- B. Define common set of variants and genes for analysis ---
geno_overlap <- fread(paste0(mashr_dir, "step0_overlap_EUR.bim"))$V2
gene_list <- readRDS(paste0(mashr_dir, "step1_gene_list_", ct, ".rds"))

# Keep only genes testable in all tissues to ensure a fair comparison
count <- table(unlist(gene_list))
n_max <- max(count)
gene_keep <- names(count)[count == n_max]
print(paste("Identified", length(gene_keep), "genes testable in all", n_max, "tissues."))

# --- C. Iterate through chromosomes and perform heterogeneity test ---
i2_results_list <- list()
for (chrom in 1:22) {
    print(paste0("--- Processing Chromosome ", chrom, " ---"))

    summary_stats <- load_summary_statistics(ct, chrom, geno_overlap, gene_keep, ct_cross_tissue)
    
    if (is.null(summary_stats) || nrow(summary_stats$df_all_slope_wide) == 0) {
        print(paste("No overlapping significant pairs found on chr", chrom, ". Skipping."))
        next
    }
    
    # Filter for significant pairs on the current chromosome
    df_slope_wide <- merge(summary_stats$df_all_slope_wide, sig_pairs, by = c("phenotype_id", "variant_id"))
    df_se_wide <- merge(summary_stats$df_all_se_wide, sig_pairs, by = c("phenotype_id", "variant_id"))
    
    if (nrow(df_slope_wide) == 0) {
        print(paste("No significant pairs to process on chr", chrom, ". Skipping."))
        next
    }

    df_pairs <- df_slope_wide[, .(phenotype_id, variant_id)]
    effect_sizes <- as.matrix(df_slope_wide[, -c(1:2)])
    std_errors <- as.matrix(df_se_wide[, -c(1:2)])
    
    # Clean up memory
    rm(summary_stats, df_slope_wide, df_se_wide)
    gc()
    
    print(paste("Calculating I^2 for", nrow(df_pairs), "pairs..."))
    start_time <- proc.time()

    # Fit a random-effects model for each eQTL pair to estimate heterogeneity
    # We use mclapply for parallel execution of this step.
    het_stats <- mclapply(1:nrow(effect_sizes), function(i) {
        # NA values mean the eQTL was not tested/available in that tissue; rma handles this.
        res <- tryCatch({
            rma.uni(yi = effect_sizes[i, ],
                    sei = std_errors[i, ],
                    method = "REML", # REML is a standard method for random-effects models
                    control = list(maxiter = 1000))
        }, error = function(e) NULL)
        
        if (is.null(res)) {
            return(list(I2 = NA, Q = NA, QEp = NA))
        } else {
            return(list(I2 = res$I2, Q = res$QE, QEp = res$QEp))
        }
    }, mc.cores = num_cores)
    
    end_time <- proc.time()
    print(paste("Time elapsed:", round((end_time - start_time)[3], 2), "seconds"))
    
    # Append results
    df_pairs$I2 <- sapply(het_stats, function(x) x$I2)
    df_pairs$Q <- sapply(het_stats, function(x) x$Q)
    df_pairs$QEp <- sapply(het_stats, function(x) x$QEp)
    
    i2_results_list[[chrom]] <- df_pairs
}


## 4. SAVE RESULTS
# ------------------------------------------------------------------------------
print("Finalizing and saving results.")
all_i2_results <- rbindlist(i2_results_list)
output_file <- file.path(output_dir, paste0("I2_REML_", ct, ".rds"))

saveRDS(all_i2_results, output_file)

print(paste("Analysis complete. Heterogeneity results saved to:", output_file))