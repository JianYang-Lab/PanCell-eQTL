################################################################################
#
# --- Apply Pre-trained mashR Model to All eQTLs by Chromosome ---
#
# Description:
# This script applies a pre-trained mashR model to all SNP-gene pairs for a 
# given cell type, processing the data one chromosome at a time. It first
# loads and formats the complete summary statistics for the target chromosome,
# then uses the fitted model (from a previous step) to compute posterior 
# summaries (e.g., effect sizes, lfsr) for every eQTL test.
#
# This script is designed to be run as a batch job for each chromosome.
#
# Input:
# 1. A pre-trained mashR model object (`m`).
# 2. A file containing overlapping SNPs to use for filtering.
# 3. A file containing lists of tested genes for each condition.
# 4. Full eQTL summary statistics (.parquet files) for each tissue.
#
# Output:
# A mashR results object (`m2`) containing posterior summaries for all 
# eQTLs on the specified chromosome, saved as an .rds file.
#
################################################################################


## 1. SETUP
# ------------------------------------------------------------------------------
# Load required libraries
library(tidyverse)
library(data.table)
library(arrow)
library(parallel)
library(mashr)

# --- User-Defined Parameters ---

# Set the main project directory. All other paths will be relative to this.
base_dir <- "/path/to/your/scRNA_project_folder/"
setwd(base_dir)

# Define the cell type and chromosome for analysis.
# These are typically passed as command-line arguments for batch jobs.
# ct <- "T_cells" 
# chrom <- 22
args <- commandArgs(trailingOnly = TRUE)
ct <- args[1]
chrom <- as.numeric(args[2])

if (is.na(ct) || is.na(chrom)) {
  stop("Cell type or chromosome not provided. Usage: Rscript script.R <cell_type> <chromosome>")
}
print(paste0("🔬 Applying mashR model to ", ct, " on chromosome ", chrom))

# Define the list of tissues and the ancestry to analyze
tissue_list <- c("Blood", "Lung", "Skin", "Colon", "Liver")
ancestry <- "EUR"

# Set the number of cores for parallel processing
num_cores <- detectCores() - 2

# Define paths for mashR-related intermediate files and outputs
mashr_dir <- "mashr/cross_tissue/"
dir.create(mashr_dir, showWarnings = FALSE, recursive = TRUE)


## 2. HELPER FUNCTIONS
# ------------------------------------------------------------------------------

#' Convert p-value to Z-score with direction from beta
ConvertP2Z <- function(pval, beta) {
    z <- abs(qnorm(pval / 2))
    z[beta < 0] <- -1 * z[beta < 0]
    return(z)
}

#' Handle missing effect sizes (beta)
handle_nan_b <- function(x) {
    x[is.nan(x) | is.na(x)] <- 0
    return(x)
}

#' Handle problematic standard errors (SE)
handle_nan_s <- function(x) {
    x[is.nan(x) | is.infinite(x) | is.na(x) | x == 0] <- 1E3
    return(x)
}

#' Handle missing Z-scores based on beta and SE
handle_nan_z <- function(z, b, s) {
    z[is.nan(s) | is.infinite(s) | is.na(s) | s == 0] <- 0
    z[is.nan(b) | is.na(b) | b == 0] <- 0
    z[is.nan(z) | is.na(z)] <- 0
    return(z)
}

#' Process a data.table of summary statistics
GetSS <- function(dat) {
    dat$"z" <- ConvertP2Z(dat$pval_nominal, dat$slope)
    dat[['z']] <- handle_nan_z(dat[['z']], dat[['slope']], dat[['slope_se']])
    dat[['slope_se']] <- handle_nan_s(dat[['slope_se']])
    dat[['slope']] <- handle_nan_b(dat[['slope']])
    return(dat)
}

#' Load and filter summary statistics for one tissue-chromosome combination
process_combinations <- function(combo, gene_keep, geno_overlap) {
    tissue <- combo$tissue
    ct_raw <- combo$ct_raw
    anc <- combo$anc
    chrom <- combo$chrom
    
    file_path <- file.path(tissue, "sc-eQTL/results/raw", 
                           paste0(ct_raw, "_", anc, "_pc5.cis_qtl_pairs.", chrom, ".parquet"))
    
    if (file.exists(file_path)) {
        print(paste("... loading data for", tissue, ct_raw))
        df <- read_parquet(file_path)
        # Filter for overlapping genes/SNPs and reasonable allele frequency
        df <- df[df$phenotype_id %in% gene_keep & df$variant_id %in% geno_overlap, ]
        df <- df[df$af >= 0.01 & df$af <= 0.99, ]
        df$ancestry <- anc
        df$tissue <- tissue
        return(df)
    }
    return(NULL)
}


## 3. LOAD AND PREPARE SUMMARY STATISTICS
# ------------------------------------------------------------------------------
summary_stats_file <- file.path(mashr_dir, paste0("summary_stats_", ct, "_chr", chrom, ".rds"))

if (!file.exists(summary_stats_file)) {
    print("Step 3.1: Loading and processing summary statistics from scratch.")
    
    # Load prerequisite files for filtering
    geno_overlap <- fread(file.path(mashr_dir, paste0("step_overlap_", ancestry, ".bim")), header = F)$V1
    gene_list <- readRDS(file.path(mashr_dir, paste0("step_gene_list_", ct, ".rds")))
    ct_cross_tissue <- read_delim("celltype_mappings.txt", delim = "\t")
    
    # Identify genes that were tested in all conditions for comparability
    count <- table(unlist(gene_list))
    gene_keep <- names(count)[count == max(count)]

    # Create a list of all tissue/cell type combinations to process
    combinations <- list()
    for (ts in tissue_list) {
        ct_raw_name <- ct_cross_tissue %>% 
            dplyr::filter(tissue == ts & celltype == ct) %>% 
            pull(celltype_raw)
        
        if (length(ct_raw_name) > 0) {
            combinations[[length(combinations) + 1]] <- list(
                anc = ancestry,
                ct_raw = ct_raw_name,
                tissue = ts,
                chrom = chrom
            )
        }
    }
    
    # Load and filter data in parallel
    df_all_list <- mclapply(combinations, process_combinations, 
                            gene_keep = gene_keep, geno_overlap = geno_overlap, mc.cores = num_cores)

    df_all_list <- Filter(Negate(is.null), df_all_list)
    if (length(df_all_list) == 0) {
        stop("No data loaded. Check paths and file existence.")
    }

    # Clean and combine summary statistics
    df_all_ss <- mclapply(df_all_list, GetSS, mc.cores = num_cores)
    df_all_ss <- rbindlist(df_all_ss)

    # Reshape data from long to wide format (tests x conditions)
    bhat_wide <- dcast(df_all_ss, phenotype_id + variant_id ~ tissue, value.var = "slope")
    shat_wide <- dcast(df_all_ss, phenotype_id + variant_id ~ tissue, value.var = "slope_se")

    # Save the processed data for future runs
    saveRDS(list(bhat = bhat_wide, shat = shat_wide), summary_stats_file)
    summary_stats <- list(bhat = bhat_wide, shat = shat_wide)
    
} else {
    print("Step 3.1: Loading pre-processed summary statistics.")
    summary_stats <- readRDS(summary_stats_file)
}


## 4. APPLY FITTED MASH MODEL
# ------------------------------------------------------------------------------
print("Step 4.1: Setting up data for mashR.")

# Load the previously fitted mash model
fitted_model_file <- file.path(mashr_dir, paste0("step_mash_fit_", ct, ".rds"))
if (!file.exists(fitted_model_file)) {
    stop("Fitted mash model not found. Please run the model training script first.")
}
m <- readRDS(fitted_model_file)

# Load the random data used for training to estimate the null correlation (V)
# This ensures the null correlation structure is identical to the one used for model fitting.
random_data_file <- file.path(mashr_dir, paste0("step_strong_random_for_mash_", ct, ".rds"))
out <- readRDS(random_data_file)
data_temp <- mash_set_data(Bhat = as.matrix(out$random_b), Shat = as.matrix(out$random_s))
Vhat <- estimate_null_correlation_simple(data_temp)
rm(data_temp, out)

# Set up the new data object for the current chromosome
# We provide the estimated null correlation matrix Vhat.
bhat_matrix <- as.matrix(summary_stats$bhat[, -c(1:2)])
shat_matrix <- as.matrix(summary_stats$shat[, -c(1:2)])
data_all_pairs <- mash_set_data(Bhat = bhat_matrix, Shat = shat_matrix, V = Vhat)

print("Step 4.2: Computing posterior summaries.")
# Apply the fitted model from 'm' to the new data.
# `get_fitted_g(m)` extracts the learned mixture of covariance matrices.
# `fixg = TRUE` tells mash to use this existing model instead of re-fitting it.
m2_posteriors <- mash(data_all_pairs, g = get_fitted_g(m), fixg = TRUE)

# Save the final results object, which contains all posterior information
output_file <- file.path(mashr_dir, paste0("mash_posteriors_", ct, "_chr", chrom, ".rds"))
saveRDS(m2_posteriors, output_file)

print(paste("✅ Analysis complete! Posterior results for chr", chrom, "saved to:", output_file))