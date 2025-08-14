################################################################################
#
# --- Apply Pre-trained Intra-Tissue mashR Model by Chromosome ---
#
# Description:
# This script applies a pre-trained mashR model (specific to a single tissue)
# to all SNP-gene pairs on a given chromosome. It loads the complete summary
# statistics for that chromosome, formats the data, and then uses the fitted
# model to compute posterior summaries for every eQTL test.
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
base_path <- "/path/to/"
mashr_path <- file.path(base_path, "crossTissue_analysis/mashR/")

# Define the tissue and chromosome for analysis from command line arguments
args <- commandArgs(TRUE)
tissue <- args[1]
chrom <- args[2]
# tissue <- "Blood"   # Manual override for interactive use
# chrom <- "22"       # Manual override for interactive use

if (is.na(tissue) || is.na(chrom)) {
  stop("Tissue or chromosome not provided. Usage: Rscript script.R <tissue> <chromosome_number>")
}

tissue_path <- file.path(base_path, tissue, "sc-eQTL/results/")
print(paste0("🔬 Applying mashR model to ", tissue, " on chromosome ", chrom))


## 2. HELPER FUNCTIONS
# ------------------------------------------------------------------------------
print("--- Section 2: Defining Helper Functions ---")

ConvertP2Z <- function(pval, beta) {
    z <- abs(qnorm(pval / 2))
    z[beta < 0] <- -1 * z[beta < 0]
    return(z)
}

handle_nan_b <- function(x) {
    x[which(is.nan(x) | is.na(x))] <- 0
    return(x)
}

handle_nan_s <- function(x) {
    x[which(is.nan(x) | is.infinite(x) | is.na(x) | x == 0)] <- 1E3
    return(x)
}

handle_nan_z <- function(z, b, s) {
    z[which(is.nan(s) | is.infinite(s) | is.na(s) | s == 0)] <- 0
    z[which(is.nan(b) | is.na(b) | b == 0)] <- 0
    z[which(is.nan(z) | is.na(z))] <- 0
    return(z)
}

GetSS <- function(dat) {
    dat$"z" <- ConvertP2Z(dat$"pval_nominal", dat$"slope")
    dat[['z']] <- handle_nan_z(dat[['z']], dat[['slope']], dat[['slope_se']])
    dat[['slope_se']] <- handle_nan_s(dat[['slope_se']])
    dat[['slope']] <- handle_nan_b(dat[['slope']])
    return(dat)
}

# This function loads and filters data for one combination.
# NOTE: It relies on several global variables defined later in the script:
# `smr_list_all`, `cell_type_short_dict`, `tissue_path`, `gene_keep`, `geno_overlap`.
process_combinations <- function(ct_anc_combo) {
  ct <- ct_anc_combo$ct
  anc <- ct_anc_combo$anc
  tissue <- ct_anc_combo$tissue
  chrom <- ct_anc_combo$chrom

  ct_anc_tmp <- paste0(cell_type_short_dict[[tissue]][ct], "_", anc)
  if (ct_anc_tmp %in% smr_list_all$ct_anc[smr_list_all$tissue == tissue]) {
    file_path <- file.path(tissue_path, "raw", paste0(ct, "_", anc, "_pc5.cis_qtl_pairs.", chrom, ".parquet"))

    if (file.exists(file_path)) {
      print(paste0(anc, "_", ct))
      df <- read_parquet(file_path)
      df <- df[df$phenotype_id %in% gene_keep & df$variant_id %in% geno_overlap, ]
      df <- df[df$af >= 0.01 & df$af <= 0.99, ]
      df$ancestry <- anc
      df$celltype <- ct
      return(df)
    }
  }
  return(NULL) # Return NULL if the combination does not meet criteria
}


## 3. LOAD PREREQUISITE METADATA
# ------------------------------------------------------------------------------
print("--- Section 3: Loading Prerequisite Metadata ---")

# The `process_combinations` function requires this metadata to filter for valid conditions.
tissue_list <- c("Blood", "Lung", "Skin", "Colon", "Liver")
ancestry_list <- c("EUR", "EAS", "AFR", "AMR")
cell_type_cat_colors <- read_delim(file.path(base_path, "crossTissue_analysis/cell_type_colors_new.txt"), delim = '\t')
ts_ct_anc <- read_delim(file.path(base_path, "crossTissue_analysis/tissue_celltype_ancestry_ss_new.txt"), delim = '\t')

cell_type_short_dict <- lapply(tissue_list, function(x) {
  tmp <- cell_type_cat_colors$cell_type_short[cell_type_cat_colors$tissue == x]
  names(tmp) <- cell_type_cat_colors$cell_type[cell_type_cat_colors$tissue == x]
  tmp
})
names(cell_type_short_dict) <- tissue_list

smr_list <- list()
for (t in tissue_list) {
  smr <- read_table(file.path(base_path, t, "sc-eQTL/results/sig/eQTL_summary_5en8_maf01.tsv"))
  smr$tissue <- t
  smr_list[[t]] <- smr
}

smr_list_all <- list_rbind(smr_list) %>%
  left_join(ts_ct_anc, by = c("tissue", "cell_type", "ancestry")) %>%
  filter(!(tissue == "Blood" & cell_type %in% c("T", "B", "Eryth"))) %>%
  filter(sample_size >= 50) %>%
  inner_join(cell_type_cat_colors, by = c("tissue", "cell_type")) %>%
  mutate(
    ancestry = factor(ancestry, levels = ancestry_list),
    ct_anc = paste0(cell_type_short, "_", ancestry),
    total_count = round(sample_size * mean_count)
  ) %>%
  filter(total_count >= 10000)

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


## 4. LOAD AND PREPARE SUMMARY STATISTICS FOR THE TARGET CHROMOSOME
# ------------------------------------------------------------------------------
print(paste0("--- Section 4: Loading Summary Statistics for Chr ", chrom, " ---"))

summary_stats_file <- file.path(mashr_path, paste0("step5_summary_statistics_", tissue, "_chr", chrom, ".rds"))

if (!file.exists(summary_stats_file)) {
  # Load files needed for filtering within the helper function
  geno_overlap <- read.table(paste0(mashr_path, "step3_", tissue, "_overlap.bim"))$V1
  top_all <- readRDS(paste0(mashr_path, "step1_top_all_", tissue, ".rds"))
  gene_list <- readRDS(paste0(mashr_path, "step1_gene_list_", tissue, ".rds"))
  
  # Define genes to keep based on presence across conditions
  count <- table(unlist(gene_list))
  n_max <- max(count)
  gene_keep <- names(count)[count >= n_max - length(unique(top_all$ancestry))]
  rm(top_all)

  print(paste0("Reading in raw data for ", tissue, " chr", chrom, "."))

  # Create a list of all possible ancestry and cell type combinations to check
  ct_anc_combinations <- list()
  for (anc in ancestries_all_list[[tissue]]) {
    for (ct in celltypes_all_list[[tissue]]) {
      ct_anc_combinations[[length(ct_anc_combinations) + 1]] <- list(
        anc = anc, ct = ct, tissue = tissue, chrom = chrom
      )
    }
  }

  # Process all combinations in parallel
  df_all_list <- mclapply(ct_anc_combinations, process_combinations, mc.cores = 2)
  df_all_list <- Filter(Negate(is.null), df_all_list)

  # Clean and combine summary statistics
  df_all_ss <- rbindlist(mclapply(df_all_list, GetSS, mc.cores = 2))

  # Reshape data from long to wide format
  df_all_slope_wide <- dcast(df_all_ss, phenotype_id + variant_id ~ celltype + ancestry, value.var = "slope")
  df_all_se_wide <- dcast(df_all_ss, phenotype_id + variant_id ~ celltype + ancestry, value.var = "slope_se")
  df_all_z_wide <- dcast(df_all_ss, phenotype_id + variant_id ~ celltype + ancestry, value.var = "z")

  # Save intermediate files
  df_all_pair <- df_all_z_wide[, 1:2]
  saveRDS(df_all_pair, paste0(mashr_path, "step5_all_pairs_", tissue, "_chr", chrom, ".rds"))
  saveRDS(list(df_all_slope_wide = df_all_slope_wide, df_all_se_wide = df_all_se_wide), summary_stats_file)
  
} else {
  print("Loading pre-processed summary statistics from file.")
  summary_stats <- readRDS(summary_stats_file)
  df_all_slope_wide <- summary_stats$df_all_slope_wide
  df_all_se_wide <- summary_stats$df_all_se_wide
}

ncond <- ncol(df_all_slope_wide)


## 5. APPLY FITTED MASH MODEL
# ------------------------------------------------------------------------------
print("--- Section 5: Applying Fitted mashR Model ---")

# Load the random data set to re-calculate the null correlation matrix
out <- readRDS(paste0(mashr_path, "step4_strong_random_for_mash_", tissue, ".rds"))
data.temp <- mash_set_data(Bhat = as.matrix(out$random_b), Shat = as.matrix(out$random_s))
Vhat <- estimate_null_correlation_simple(data.temp)
rm(data.temp)

# Prepare the data for the current chromosome with the estimated null correlation
data_STRONG <- mash_set_data(
  as.matrix(df_all_slope_wide[, 3:ncond]),
  as.matrix(df_all_se_wide[, 3:ncond]),
  V = Vhat
)

# Load the pre-trained mash model
print("Read in fitted model.")
m <- readRDS(paste0(mashr_path, "step5_m_", tissue, ".rds"))

# Compute posterior summaries for all pairs on the chromosome
print("Computing posterior summaries.")
m2 <- mash(data_STRONG, g = get_fitted_g(m), fixg = TRUE)

# Save the final posterior results object
saveRDS(m2, paste0(mashr_path, "step5_m2_all_posterior_", tissue, "_chr", chrom, ".rds"))

print(paste("✅ Analysis complete for", tissue, "chr", chrom))