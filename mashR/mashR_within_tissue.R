################################################################################
#
# --- Intra-Tissue eQTL Sharing Analysis with mashR ---
#
# Description:
# This script uses the 'mashR' package to analyze sharing patterns of sc-eQTLs
# across multiple cell types and ancestries within a single specified tissue.
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
library(ComplexHeatmap)
library(UpSetR)
library(circlize)
library(corrplot)
library(rmeta)
library(ashr)

# --- User-Defined Parameters ---

# Define base paths
mashr_path <- "/path/to/mashR/intra_tissue/"
base_path <- "/path/to/"

# Define the tissue for analysis (from command line argument)
args <- commandArgs(TRUE)
tissue <- args[1]
# tissue <- "Blood" # Manual override for interactive use

if (is.na(tissue)) {
  stop("Tissue not provided. Please specify a tissue to analyze.")
}

tissue_path <- paste0(base_path, tissue, "/sc-eQTL/results/")
print(paste0("Conducting mashR for ", tissue))


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


## 3. LOAD METADATA & DEFINE CONDITIONS
# ------------------------------------------------------------------------------
print("--- Section 3: Loading Metadata & Defining Conditions ---")

# Define lists of tissues, ancestries, and cell types
tissue_list <- c("Blood", "Lung", "Skin", "Colon", "Liver")
ancestry_list <- c("EUR", "EAS", "AFR", "AMR")
ancestries_all_list <- list(
  Blood = c("EUR", "EAS", "AFR", "AMR"), Lung = c("EUR", "EAS"),
  Skin = c("EUR", "EAS"), Colon = c("EUR"), Liver = c("EUR", "EAS")
)
celltypes_all_list <- list(
  Blood = c('CD4TNC','CD4TEM','Treg','CD8TNC','CD8TEM','CD8TEMRA','MAIT','NKp','NKn','BIN','BMem','Plasma','Plasmablasts','MonoC','MonoNC','Nph','DC','pDC'),
  Liver = c('CD4T',"CD8T","Treg","MAIT","Th","gdT",'Circulating_NK','Resident_NK','B','Plasma','Monocytes','Macrophages','DC','pDC','Neutrophils','Basophils','Hepatocytes','Cholangiocytes','Endothelium','Fibroblasts'),
  Skin = c("Th","Tc", "Treg","NK","DC1","DC2","MigDC","Macro1","Macro2","MonoMac","Mast","KCdiff","KCundiff","Melanocyte","VE1","VE2","VE3","LE1","Pericyte1","Pericyte2","F2"),
  Lung = c("CD4T","CD8T","NK","B","Monocytes","Mast","Macrophages","DC","AT1","AT2","Airway_Epi_Multiciliated","Airway_Epi_Basal","Airway_Epi_Secretory","EC_arterial","EC_capillary","EC_venous","LEC_mature","LEC_diff","Fibroblasts","Smooth_muscle","Mesothelium","Rare"),
  Colon = c("CD4T","CD8T","Treg","Th","gdT","NKT","NK","ILC3","BIN","BMem","Bcyc","Plasma","Mono","Macro","cDC2","Mast","Colonocyte","GLoblet","Tuft","TA", "EEC" ,"ECcap","ECven","Stromal1","Stromal2","Myofibroblast","Glia")
)

# Load metadata for cell type names and sample sizes
cell_type_cat_colors <- read_delim(file.path(base_path, "cell_type_colors_new.txt"), delim = '\t')
ts_ct_anc <- read_delim(file.path(base_path, "tissue_celltype_ancestry_ss_new.txt"), delim = '\t')

# Create a dictionary to map long cell type names to short names
cell_type_short_dict <- lapply(tissue_list, function(x) {
  tmp <- cell_type_cat_colors$cell_type_short[cell_type_cat_colors$tissue == x]
  names(tmp) <- cell_type_cat_colors$cell_type[cell_type_cat_colors$tissue == x]
  tmp
})
names(cell_type_short_dict) <- tissue_list

# Load eQTL summary counts from all tissues
smr_list <- list()
for (t in tissue_list) {
  smr <- read_table(file.path(base_path, t, "sc-eQTL/results/sig/eQTL_summary_5en8_maf01.tsv"))
  smr$tissue <- t
  smr_list[[t]] <- smr
}

# Combine all summaries and filter for high-quality conditions
smr_list_all <- list_rbind(smr_list) %>%
  left_join(ts_ct_anc, by = c("tissue", "cell_type", "ancestry")) %>%
  filter(sample_size >= 50) %>%
  inner_join(cell_type_cat_colors, by = c("tissue", "cell_type")) %>%
  mutate(
    ancestry = factor(ancestry, levels = ancestry_list),
    ct_anc = paste0(cell_type_short, "_", ancestry),
    total_count = round(sample_size * mean_count)
  ) %>%
  filter(total_count >= 10000)


## 4. DATA PREPARATION FOR MASH
# ------------------------------------------------------------------------------
print("--- Section 4: Preparing Data for mashR ---")

# --- Step 4.1: Find overlapping SNPs ---
print("Step 4.1: Find overlapping SNPs across ancestries.")
if (!file.exists(paste0(mashr_path, "step3_", tissue, "_overlap.bim"))) {
  geno_list <- list()
  for (anc in ancestries_all_list[[tissue]]) {
    print(anc)
    geno_list[[length(geno_list) + 1]] <- read.table(paste0(tissue_path, "../genotype/", tissue, "_", anc, "_indivpruned.bim"))$V2
  }
  geno_all <- unlist(geno_list)
  if (length(ancestries_all_list[[tissue]]) == 1) {
    geno_overlap <- geno_all
    rm(geno_all)
  } else {
    geno_all_tab <- table(geno_all)
    geno_overlap <- names(geno_all_tab)[geno_all_tab == length(ancestries_all_list[[tissue]])]
  }
  write.table(geno_overlap, paste0(mashr_path, "step3_", tissue, "_overlap.bim"), sep = '\n', col = F, row = F, quo = F)
} else {
  geno_overlap <- read.table(paste0(mashr_path, "step3_", tissue, "_overlap.bim"))$V1
}

# --- Step 4.2: Find top eQTLs per condition and all tested genes ---
print("Step 4.2: Find top eQTLs within each condition.")
if (!file.exists(paste0(mashr_path, "step1_top_all_", tissue, ".rds"))) {
  top_list <- list()
  gene_list <- list()
  for (anc in ancestries_all_list[[tissue]]) {
    for (ct in celltypes_all_list[[tissue]]) {
      ct_anc_tmp <- paste0(cell_type_short_dict[[tissue]][ct], "_", anc)
      if (ct_anc_tmp %in% smr_list_all$ct_anc[smr_list_all$tissue == tissue]) {
        if (file.exists(file.path(tissue_path, "sig", paste0(ct, "_", anc, "_pc5_sig_5en8.tsv.gz")))) {
          print(paste0(anc, "_", ct))
          df_top <- fread(file.path(tissue_path, "sig", paste0(ct, "_", anc, "_pc5_sig_5en8.tsv.gz")))
          gene_list[[length(gene_list) + 1]] <- fread(file.path(tissue_path, "raw", paste0(ct, "_", anc, "_pc5_perm.tsv")), select = "phenotype_id")$phenotype_id
          df_top <- df_top[variant_id %in% geno_overlap, .(phenotype_id, variant_id, slope, slope_se, pval_nominal)]
          df_top <- df_top[, .SD[which.min(pval_nominal)], by = phenotype_id]
          df_top$ancestry <- anc
          df_top$celltype <- ct
          top_list[[(length(top_list) + 1)]] <- df_top
        }
      }
    }
  }
  top_all <- rbindlist(top_list)
  saveRDS(top_all, paste0(mashr_path, "step1_top_all_", tissue, ".rds"))
  saveRDS(gene_list, paste0(mashr_path, "step1_gene_list_", tissue, ".rds"))
} else {
  top_all <- readRDS(paste0(mashr_path, "step1_top_all_", tissue, ".rds"))
  gene_list <- readRDS(paste0(mashr_path, "step1_gene_list_", tissue, ".rds"))
}

# --- Step 4.3: Define final "strong" and "random" sets ---
print("Step 4.3: Find the strong set (top eQTLs across all conditions).")
count <- table(unlist(gene_list))
n_max <- max(count)
gene_keep <- names(count)[count >= n_max - length(unique(top_all$ancestry))]
if (!file.exists(paste0(mashr_path, "step2_top_all_cond_", tissue, ".rds"))) {
  top_all_cond <- top_all[phenotype_id %in% gene_keep, .SD[which.min(pval_nominal)], by = phenotype_id]
  saveRDS(top_all_cond, paste0(mashr_path, "step2_top_all_cond_", tissue, ".rds"))
} else {
  top_all_cond <- readRDS(paste0(mashr_path, "step2_top_all_cond_", tissue, ".rds"))
}

# --- Step 4.4: Prepare mashR input data (strong and random sets) ---
if (!file.exists(paste0(mashr_path, "step4_strong_random_for_mash_", tissue, ".rds"))) {
    print("Step 4.4: Find the random set.")
    process_chromosome <- function(chromosome, ct) {
        print(chromosome)
        df <- read_parquet(file.path(tissue_path, "raw", paste0(ct, "_EUR_pc5.cis_qtl_pairs.", chromosome, ".parquet"))) %>% as.data.table()
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
    if (!file.exists(paste0(mashr_path, "step4_random_list_", tissue, ".rds"))) {
        tmp <- table(top_all_cond$celltype)
        ct <- names(tmp)[which.max(tmp)]
        df_random_list_all <- mclapply(1:22, function(x) {process_chromosome(x, ct)}, mc.cores = 2)
        df_random_list <- rbindlist(unlist(df_random_list_all, recursive = FALSE))
        saveRDS(df_random_list, paste0(mashr_path, "step4_random_list_", tissue, ".rds"))
        rm(df_random_list_all)
    } else {
        df_random_list <- readRDS(paste0(mashr_path, "step4_random_list_", tissue, ".rds"))
    }

    print("... fetching summary statistics for random and strong sets (this is slow).")
    df_random_all_list <- list()
    df_strong_all_list <- list()
    for (anc in ancestries_all_list[[tissue]]) {
        for (ct in celltypes_all_list[[tissue]]) {
            ct_anc_tmp <- paste0(cell_type_short_dict[[tissue]][ct], "_", anc)
            if (ct_anc_tmp %in% smr_list_all$ct_anc[smr_list_all$tissue == tissue]) {
                if (file.exists(file.path(tissue_path, "raw", paste0(ct, "_", anc, "_pc5.cis_qtl_pairs.22.parquet")))) {
                    print(paste0(anc, "_", ct))
                    df_list <- mclapply(1:22, function(i) {
                        print(i)
                        df <- read_parquet(file.path(tissue_path, "raw", paste0(ct, "_", anc, "_pc5.cis_qtl_pairs.", i, ".parquet"))) %>% as.data.table()
                        df_random_select <- merge(df, df_random_list, by = c("phenotype_id", "variant_id"))[, .(phenotype_id, variant_id, slope, slope_se, pval_nominal)]
                        df_strong_select <- merge(df, top_all_cond[, .(phenotype_id, variant_id)], by = c("phenotype_id", "variant_id"))[, .(phenotype_id, variant_id, slope, slope_se, pval_nominal)]
                        return(list(random = df_random_select, strong = df_strong_select))
                    }, mc.cores = 2)

                    df_list_random <- rbindlist(lapply(df_list, function(x) x$random))
                    df_list_strong <- rbindlist(lapply(df_list, function(x) x$strong))
                    df_list_random$ancestry <- anc
                    df_list_random$celltype <- ct
                    df_list_strong$ancestry <- anc
                    df_list_strong$celltype <- ct
                    df_random_all_list[[length(df_random_all_list) + 1]] <- df_list_random
                    df_strong_all_list[[length(df_strong_all_list) + 1]] <- df_list_strong
                }
            }
        }
    }

    print("... cleaning and reshaping data.")
    df_random_all_ss <- rbindlist(mclapply(df_random_all_list, GetSS, mc.cores = 6))
    saveRDS(df_random_all_ss, paste0(mashr_path, "step4_random_all_ss_", tissue, ".rds"))
    df_strong_all_ss <- rbindlist(mclapply(df_strong_all_list, GetSS, mc.cores = 6))
    saveRDS(df_strong_all_ss, paste0(mashr_path, "step4_strong_all_ss_", tissue, ".rds"))

    df_random_all_z_wide <- dcast(df_random_all_ss, phenotype_id + variant_id ~ celltype + ancestry, value.var = "z")
    df_random_all_slope_wide <- dcast(df_random_all_ss, phenotype_id + variant_id ~ celltype + ancestry, value.var = "slope")
    df_random_all_se_wide <- dcast(df_random_all_ss, phenotype_id + variant_id ~ celltype + ancestry, value.var = "slope_se")
    df_strong_all_z_wide <- dcast(df_strong_all_ss, phenotype_id + variant_id ~ celltype + ancestry, value.var = "z")
    df_strong_all_slope_wide <- dcast(df_strong_all_ss, phenotype_id + variant_id ~ celltype + ancestry, value.var = "slope")
    df_strong_all_se_wide <- dcast(df_strong_all_ss, phenotype_id + variant_id ~ celltype + ancestry, value.var = "slope_se")

    ncond <- ncol(df_random_all_z_wide)
    for (i in 3:ncond) {
        df_random_all_slope_wide[[i]] <- handle_nan_b(df_random_all_slope_wide[[i]])
        df_random_all_se_wide[[i]] <- handle_nan_s(df_random_all_se_wide[[i]])
        df_random_all_z_wide[[i]] <- handle_nan_z(df_random_all_z_wide[[i]], df_random_all_slope_wide[[i]], df_random_all_se_wide[[i]])
        df_strong_all_slope_wide[[i]] <- handle_nan_b(df_strong_all_slope_wide[[i]])
        df_strong_all_se_wide[[i]] <- handle_nan_s(df_strong_all_se_wide[[i]])
        df_strong_all_z_wide[[i]] <- handle_nan_z(df_strong_all_z_wide[[i]], df_strong_all_slope_wide[[i]], df_strong_all_se_wide[[i]])
    }

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
    saveRDS(out, paste0(mashr_path, "step4_strong_random_for_mash_", tissue, ".rds"))
} else {
    out <- readRDS(paste0(mashr_path, "step4_strong_random_for_mash_", tissue, ".rds"))
}


## 5. MASH ANALYSIS
# ------------------------------------------------------------------------------
print("--- Section 5: Running mashR Analysis ---")

# --- Step 5.1: Estimate null correlation and set up data ---
print("Step 5.1: Estimate null correlation structure.")
data.temp <- mash_set_data(Bhat = as.matrix(out$random_b), Shat = as.matrix(out$random_s))
Vhat <- estimate_null_correlation_simple(data.temp)

data_RANDOM <- mash_update_data(data.temp, V = Vhat)
data_STRONG <- mash_set_data(as.matrix(out$strong_b), as.matrix(out$strong_s), V = Vhat)

# --- Step 5.2: Set up covariance matrices ---
print("Step 5.2: Set up data-driven and canonical covariance matrices.")
U.pca <- cov_pca(data_STRONG, 5)
U.ed <- cov_ed(data_STRONG, U.pca)
U.c <- cov_canonical(data_RANDOM)
Ulist <- c(U.ed, U.c)

# --- Step 5.3: Fit the mash model ---
print("Step 5.3: Fit mash model on the random set.")
m <- mash(data_RANDOM, Ulist = Ulist, outputlevel = 1)
saveRDS(m, paste0(mashr_path, "step5_m_", tissue, ".rds"))

# --- Step 5.4: Compute posterior summaries ---
print("Step 5.4: Compute posterior summaries on the strong set.")
m2 <- mash(data_STRONG, g = get_fitted_g(m), fixg = TRUE)
saveRDS(m2, paste0(mashr_path, "step5_m2_", tissue, ".rds"))

# --- Step 5.5: Calculate pairwise sharing ---
print("Step 5.5: Assess pairwise sharing.")
m.pairwise_PM <- get_pairwise_sharing(m2, lfsr_thresh = 0.1, factor = 0.5)
saveRDS(m.pairwise_PM, paste0(mashr_path, "step5_pairwise_PM_", tissue, ".rds"))

print(paste("✅ Analysis complete! Results for", tissue, "saved in:", mashr_path))