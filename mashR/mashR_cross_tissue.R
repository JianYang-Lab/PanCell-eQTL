################################################################################
#
# --- Cross-Tissue eQTL Sharing Analysis with mashR ---
#
# Description:
# This script uses 'mashR' to analyze the sharing patterns of single-cell eQTLs
# for a specific cell type across multiple tissues, focusing on EUR ancestry.
#
################################################################################


## 1. SETUP
# ------------------------------------------------------------------------------
print("--- Section 1: Setup ---")

# Load required libraries
library(tidyverse)
library(data.table)
library(mashr)
library(arrow)
library(parallel)
library(ashr)

# --- User-Defined Parameters ---

# Define base paths
base_path <- "/path/to/your/project/"
mashr_path <- file.path(base_path, "mashR/cross_tissue/")

# Define the list of tissues to analyze
tissue_list <- c("Blood", "Lung", "Skin", "Colon", "Liver")
ancestry <- "EUR"

# Define the cell type for analysis from a command line argument
args <- commandArgs(TRUE)
ct <- args[1]
# ct <- "T_cells" # Manual override for interactive use

if (is.na(ct)) {
    stop("Cell type not provided. Please specify a cell type to analyze.")
}
print(paste0("Conducting mashR for ", ct))

# Load cross-tissue cell type mapping file
ct_cross_tissue <- read_delim(file.path(base_path, "cross_tissue_celltype.txt"), delim = "\t")


## 2. HELPER FUNCTIONS
# ------------------------------------------------------------------------------
print("--- Section 2: Defining Helper Functions ---")

ConvertP2Z <- function(pval, beta) {
    z <- abs(qnorm(pval / 2))
    z[which(beta < 0)] <- -1 * z[which(beta < 0)]
    return(z)
}

handle_nan_b <- function(x) {
    x[which(is.nan(x) | is.na(x))] = 0
    return(x)
}

handle_nan_s <- function(x) {
    x[which(is.nan(x) | is.infinite(x) | is.na(x) | x == 0)] = 1E3
    return(x)
}

handle_nan_z <- function(z, b, s) {
    z[which(is.nan(s) | is.infinite(s) | is.na(s) | s == 0)] = 0
    z[which(is.nan(b) | is.na(b) | b == 0)] = 0
    z[which(is.nan(z) | is.na(z))] = 0
    return(z)
}

GetSS <- function(dat) {
    dat$"z" <- ConvertP2Z(dat$"pval_nominal", dat$"slope")
    dat[['z']] <- handle_nan_z(dat[['z']], dat[['slope']], dat[['slope_se']])
    dat[['slope_se']] <- handle_nan_s(dat[['slope_se']])
    dat[['slope']] <- handle_nan_b(dat[['slope']])
    return(dat)
}


## 3. DATA PREPARATION FOR MASH
# ------------------------------------------------------------------------------
print("--- Section 3: Preparing Data for mashR ---")

# --- Step 3.1: Find overlapping SNPs across all tissues ---
print("Step 3.1: Find overlapping SNPs across tissues.")
if (!file.exists(paste0(mashr_path, "step0_overlap_EUR.bim"))) {
    geno_list <- list()
    for (ts in tissue_list) {
        print(ts)
        tissue_path <- paste0(base_path, ts, "/sc-eQTL/results/")
        geno_list[[length(geno_list) + 1]] <- read.table(paste0(tissue_path, "../genotype/", ts, "_", ancestry, "_indivpruned.bim"))$V2
    }
    geno_all <- unlist(geno_list)
    geno_all_tab <- table(geno_all)
    geno_overlap <- names(geno_all_tab)[geno_all_tab == length(geno_list)]
    write.table(geno_overlap, paste0(mashr_path, "step0_overlap_EUR.bim"), sep = '\n', col = F, row = F, quo = F)
} else {
    geno_overlap <- read.table(paste0(mashr_path, "step0_overlap_EUR.bim"))$V1
}

# --- Step 3.2: Find top eQTLs per tissue and all tested genes ---
print("Step 3.2: Find top eQTLs within each tissue.")
if (!file.exists(paste0(mashr_path, "step1_top_all_", ct, ".rds"))) {
    top_list <- list()
    gene_list <- list()
    for (ts in tissue_list) {
        tissue_path <- paste0(base_path, ts, "/sc-eQTL/results/")
        ct_raw <- ct_cross_tissue %>% dplyr::filter(tissue == ts & celltype == ct) %>% pull(celltype_raw)
        
        if (file.exists(file.path(tissue_path, "raw", paste0(ct_raw, "_", ancestry, "_pc5.cis_qtl_pairs.22.parquet")))) {
            print(paste0(ts, "_", ancestry, "_", ct))
            df_top <- fread(file.path(tissue_path, "sig", paste0(ct_raw, "_", ancestry, "_pc5_sig_5en8.tsv.gz")))
            gene_list[[length(gene_list) + 1]] <- fread(file.path(tissue_path, "raw", "old", paste0(ct_raw, "_", ancestry, "_pc5_perm.tsv")), select = "phenotype_id")$phenotype_id
            df_top <- df_top[variant_id %in% geno_overlap, .(phenotype_id, variant_id, slope, slope_se, pval_nominal)]
            df_top <- df_top[, .SD[which.min(pval_nominal)], by = phenotype_id]
            df_top$tissue <- ts
            df_top$ancestry <- ancestry
            top_list[[(length(top_list) + 1)]] <- df_top
        }
    }
    top_all <- rbindlist(top_list)
    saveRDS(top_all, paste0(mashr_path, "step1_top_all_", ct, ".rds"))
    saveRDS(gene_list, paste0(mashr_path, "step1_gene_list_", ct, ".rds"))
} else {
    top_all <- readRDS(paste0(mashr_path, "step1_top_all_", ct, ".rds"))
    gene_list <- readRDS(paste0(mashr_path, "step1_gene_list_", ct, ".rds"))
}

# --- Step 3.3: Define final "strong" and "random" sets ---
print("Step 3.3: Find the strong set (top eQTLs across all tissues).")
count <- table(unlist(gene_list))
n_max <- max(count)
gene_keep <- names(count)[count == n_max]
if (!file.exists(paste0(mashr_path, "step2_top_all_cond_", ct, ".rds"))) {
    top_all_cond <- top_all[phenotype_id %in% gene_keep, .SD[which.min(pval_nominal)], by = phenotype_id]
    top_all_cond <- top_all_cond[pval_nominal < 5e-8]
    saveRDS(top_all_cond, paste0(mashr_path, "step2_top_all_cond_", ct, ".rds"))
} else {
    top_all_cond <- readRDS(paste0(mashr_path, "step2_top_all_cond_", ct, ".rds"))
}

# --- Step 3.4: Prepare mashR input data (strong and random sets) ---
if (!file.exists(paste0(mashr_path, "step4_strong_random_for_mash_", ct, ".rds"))) {
    print("Step 3.4: Find the random set.")
    process_chromosome <- function(chromosome) {
        print(chromosome)
        tissue_path <- paste0(base_path, "Lung/sc-eQTL/results/")
        ct_raw <- ct_cross_tissue %>% dplyr::filter(tissue == "Lung" & celltype == ct) %>% pull(celltype_raw)
        df <- read_parquet(file.path(tissue_path, "raw", paste0(ct_raw, "_EUR_pc5.cis_qtl_pairs.", chromosome, ".parquet"))) %>% as.data.table()
        df <- df[phenotype_id %in% gene_keep & variant_id %in% geno_overlap,]
        process_gene <- function(gene) {
            df_sub <- df[phenotype_id == gene,]
            if (nrow(df_sub) > 0) {
                df_sub <- df_sub[sample(1:nrow(df_sub), min(50, nrow(df_sub))), .(phenotype_id, variant_id)]
                return(df_sub)
            }
        }
        df_random_list <- lapply(unique(df$phenotype_id), process_gene)
        return(df_random_list)
    }
    df_random_list_all <- mclapply(1:22, process_chromosome, mc.cores = 2)
    df_random_list <- rbindlist(unlist(df_random_list_all, recursive = FALSE))
    # The full list of lists is saved here in the original script
    saveRDS(df_random_list_all, paste0(mashr_path, "step4_random_list_", ct, ".rds"))
    rm(df_random_list_all)

    print("... fetching summary statistics for random and strong sets (this is slow).")
    df_random_all_list <- list()
    df_strong_all_list <- list()
    for (ts in tissue_list) {
        tissue_path <- paste0(base_path, ts, "/sc-eQTL/results/")
        ct_raw <- ct_cross_tissue %>% dplyr::filter(tissue == ts & celltype == ct) %>% pull(celltype_raw)
        if (file.exists(file.path(tissue_path, "raw", paste0(ct_raw, "_", ancestry, "_pc5.cis_qtl_pairs.22.parquet")))) {
            print(paste0(ancestry, "_", ct, " in ", ts))
            df_list <- mclapply(1:22, function(i) {
                print(i)
                df <- read_parquet(file.path(tissue_path, "raw", paste0(ct_raw, "_", ancestry, "_pc5.cis_qtl_pairs.", i, ".parquet"))) %>% as.data.table()
                df_random_select <- merge(df, df_random_list, by = c("phenotype_id", "variant_id"))[, .(phenotype_id, variant_id, slope, slope_se, pval_nominal)]
                df_strong_select <- merge(df, distinct(top_all_cond[, .(phenotype_id, variant_id)]), by = c("phenotype_id", "variant_id"))[, .(phenotype_id, variant_id, slope, slope_se, pval_nominal)]
                return(list(random = df_random_select, strong = df_strong_select))
            }, mc.cores = 2)

            df_list_random <- rbindlist(lapply(df_list, function(x) x$random))
            df_list_strong <- rbindlist(lapply(df_list, function(x) x$strong))
            df_list_random$ancestry <- ancestry
            df_list_random$tissue <- ts
            df_list_strong$ancestry <- ancestry
            df_list_strong$tissue <- ts
            df_random_all_list[[length(df_random_all_list) + 1]] <- df_list_random
            df_strong_all_list[[length(df_strong_all_list) + 1]] <- df_list_strong
        }
    }

    print("... cleaning and reshaping data.")
    df_random_all_ss <- rbindlist(mclapply(df_random_all_list, GetSS, mc.cores = 2))
    df_strong_all_ss <- rbindlist(mclapply(df_strong_all_list, GetSS, mc.cores = 2))

    df_random_all_slope_wide <- dcast(df_random_all_ss, phenotype_id + variant_id ~ tissue + ancestry, value.var = "slope")
    df_random_all_se_wide <- dcast(df_random_all_ss, phenotype_id + variant_id ~ tissue + ancestry, value.var = "slope_se")
    df_strong_all_slope_wide <- dcast(df_strong_all_ss, phenotype_id + variant_id ~ tissue + ancestry, value.var = "slope")
    df_strong_all_se_wide <- dcast(df_strong_all_ss, phenotype_id + variant_id ~ tissue + ancestry, value.var = "slope_se")
    
    # The original script also creates Z-score matrices which are not used later, but preserved here
    df_random_all_z_wide <- dcast(df_random_all_ss, phenotype_id+variant_id~tissue+ancestry, value.var=c("z"))
    df_strong_all_z_wide <- dcast(df_strong_all_ss, phenotype_id+variant_id~tissue+ancestry, value.var=c("z"))

    # Finalize mashR input object
    ncond <- ncol(df_random_all_slope_wide)
    out <- list(
        random_pair = as.data.frame(df_random_all_z_wide)[, 1:2],
        random_z = as.data.frame(df_random_all_z_wide)[, 3:ncond],
        random_b = as.data.frame(df_random_all_slope_wide[, 3:ncond]),
        random_s = as.data.frame(df_random_all_se_wide[, 3:ncond]),
        strong_pair = as.data.frame(df_strong_all_z_wide)[, 1:2],
        strong_z = as.data.frame(df_strong_all_z_wide)[, 3:ncond],
        strong_b = as.data.frame(df_strong_all_slope_wide[, 3:ncond]),
        strong_s = as.data.frame(df_strong_all_se_wide[, 3:ncond])
    )
    saveRDS(out, paste0(mashr_path, "step4_strong_random_for_mash_", ct, ".rds"))
} else {
    out <- readRDS(paste0(mashr_path, "step4_strong_random_for_mash_", ct, ".rds"))
}


## 4. MASH ANALYSIS
# ------------------------------------------------------------------------------
print("--- Section 4: Running mashR Analysis ---")

# --- Step 4.1: Estimate null correlation and set up data ---
print("Step 4.1: Estimate null correlation structure.")
data.temp <- mash_set_data(Bhat = as.matrix(out$random_b), Shat = as.matrix(out$random_s))
Vhat <- estimate_null_correlation_simple(data.temp)
rm(data.temp)

data_RANDOM <- mash_set_data(as.matrix(out$random_b), as.matrix(out$random_s), V = Vhat)
data_STRONG <- mash_set_data(as.matrix(out$strong_b), as.matrix(out$strong_s), V = Vhat)

# --- Step 4.2: Set up covariance matrices ---
print("Step 4.2: Set up data-driven and canonical covariance matrices.")
U.pca <- cov_pca(data_STRONG, 3)
U.ed <- cov_ed(data_STRONG, U.pca)
U.c <- cov_canonical(data_RANDOM)
Ulist <- c(U.ed, U.c)

# --- Step 4.3: Fit the mash model ---
print("Step 4.3: Fit mash model on the random set.")
m <- mash(data_RANDOM, Ulist = Ulist, outputlevel = 1)
saveRDS(m, paste0(mashr_path, "step5_m_", ct, ".rds"))

# --- Step 4.4: Compute posterior summaries ---
print("Step 4.4: Compute posterior summaries on the strong set.")
m2 <- mash(data_STRONG, g = get_fitted_g(m), fixg = TRUE)
saveRDS(m2, paste0(mashr_path, "step5_m2_", ct, ".rds"))

# --- Step 4.5: Calculate pairwise sharing ---
print("Step 4.5: Assess pairwise sharing.")
m.pairwise_PM <- get_pairwise_sharing(m2, lfsr_thresh = 0.05, factor = 0.5)
saveRDS(m.pairwise_PM, paste0(mashr_path, "step5_pairwise_PM_", ct, ".rds"))

print(paste("✅ Analysis complete! Results for", ct, "saved in:", mashr_path))