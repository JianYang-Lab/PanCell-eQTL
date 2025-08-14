################################################################################
#
# Script: Replication of Indep1K sc-eQTLs in OneK1K Dataset
#
# Description:
# This script compares single-cell eQTL results from Indep1K with data from the 
# OneK1K project. The analysis is focused on corresponding blood cell types.
#
# The script performs the following steps:
#   1. Iterates through mapped cell types from both studies.
#   2. Loads data, performs liftover on OneK1K coordinates (hg19 to hg38).
#   3. Matches eQTLs between the two datasets, handling allele flipping.
#   4. Calculates standardized effect sizes for both datasets.
#   5. Saves the combined data and generates plots to visualize the results,
#      including an effect size scatter plot and a sign concordance bar plot.
#
################################################################################


## 1. SETUP & GLOBAL PARAMETERS
# ------------------------------------------------------------------------------
print("--- Section 1: Setup & Global Parameters ---")

# Load required libraries
library(tidyverse)
library(data.table)
library(rtracklayer)
library(GenomicRanges)
library(cowplot)

# --- User-Defined Paths and Parameters ---

# Set the base directory for the project. All paths will be relative to this.
base_dir <- "/path/to/your/scRNA_project_folder/"
setwd(base_dir)

# Define input/output directory for plots and data files
output_dir <- "replication_plots/onek1k/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Define paths to external datasets and resources
# Users should modify these paths to match their environment.
path_to_onek1k_data <- "/path/to/public_data/OneK1K/raw/"
path_to_liftover_chain <- "/path/to/liftover/b37ToHg38.over.chain"
path_to_sceqtl_data <- file.path(base_dir, "Blood/sc-eQTL/results/")

# --- Load Prerequisite Metadata ---
# This metadata is required for ordering and coloring plots.
cell_type_cat_colors <- read_delim("metadata/cell_type_colors_new.txt", delim = '\t')

# Create dictionaries for cell type names and colors
cell_type_short_dict <- setNames(
    cell_type_cat_colors$cell_type_short[cell_type_cat_colors$tissue == "Blood"],
    cell_type_cat_colors$cell_type[cell_type_cat_colors$tissue == "Blood"]
)
cell_type_color_dict <- cell_type_cat_colors$cell_type_colors
names(cell_type_color_dict) <- cell_type_cat_colors$cell_type_short

# --- Define Cell Type Mapping ---
# Map cell type names between our study (Indep1K) and the OneK1K project
celltype_onek1k_list <- c("cd4nc", "cd4et", "cd8nc", "cd8et", "nk", "bin", "bmem", "plasma", "monoc", "mononc", "dc")
celltype_indep1k_list <- c("CD4TNC", "CD4TEM", "CD8TNC", "CD8TEM", "NKp", "BIN", "BMem", "Plasma", "MonoC", "MonoNC", "DC")


## 2. HELPER FUNCTIONS
# ------------------------------------------------------------------------------
print("--- Section 2: Defining Helper Functions ---")

#' Calculate standardized beta and standard error from Z-score.
calcu_std_b_se <- function(z, p, n) {
  std_b_hat <- z / sqrt(2 * p * (1 - p) * (n + z^2))
  std_se <- 1 / sqrt(2 * p * (1 - p) * (n + z^2))
  res <- data.frame(std_b_hat, std_se)
  return(res)
}

#' Lift over genomic coordinates from hg19/b37 to hg38.
liftover <- function(df) {
  chain_file <- import.chain(path_to_liftover_chain)
  setnames(df, c("Chr", "BP"), c("chr_b37", "pos_b37"))
  
  # Add 'chr' prefix if missing, which is common for b37 coordinates
  snp_set <- distinct(df[, .(chr_b37, pos_b37)])
  snp_set$chr_b37 <- paste0("chr", snp_set$chr_b37)
  
  snp_gr <- GRanges(
    seqnames = snp_set$chr_b37,
    ranges = IRanges(start = snp_set$pos_b37, end = snp_set$pos_b37)
  )
  lifted_snps <- liftOver(snp_gr, chain_file)
  
  lifted_positions <- as.data.table(lifted_snps)
  old_positions <- as.data.table(snp_gr)
  old_positions[, group := .I]
  
  lifted_positions <- merge(
      old_positions[, .(chr_b37 = sub("chr", "", seqnames), pos_b37 = start, group)],
      lifted_positions[, .(chr_b38 = seqnames, pos_b38 = start, group)],
      by = "group"
  )
  
  lifted_positions[, chr_b37 := as.character(chr_b37)]
  df[, chr_b37 := as.character(chr_b37)] # Ensure join key is character
  
  df <- lifted_positions[, .(chr_b37, pos_b37, chr_b38, pos_b38)][df, on = c("chr_b37", "pos_b37")]
  return(df)
}


## 3. MAIN ANALYSIS: DATA MATCHING AND STANDARDIZATION
# ------------------------------------------------------------------------------
print("--- Section 3: Matching and Standardizing eQTL Data ---")

matched_rds_file <- file.path(output_dir, "matched_indep1k_onek1k.rds")

if (!file.exists(matched_rds_file)) {
  df_comb_all_list <- list()
  
  for (i in 1:length(celltype_indep1k_list)) {
    ct_onek1k <- celltype_onek1k_list[i]
    ct_indep1k <- celltype_indep1k_list[i]
    print(paste("Processing:", ct_indep1k, "vs.", ct_onek1k))

    # --- Load and process OneK1K data ---
    df_onek1k_ct <- fread(file.path(path_to_onek1k_data, paste0(ct_onek1k, "_eqtl_table.tsv.gz")))
    df_onek1k_ct <- df_onek1k_ct[ROUND == 1]
    setnames(df_onek1k_ct, c("GENE_ID", "CHR", "POS", "A1", "A2", "A2_FREQ_ONEK1K", "SPEARMANS_RHO", "FDR", "P_VALUE"), c("phenotype_id", "Chr", "BP", "A2", "A1", "af", "beta", "fdr", "pval_nominal"))
    
    # Perform liftover from hg19 to hg38
    df_onek1k_ct <- liftover(df_onek1k_ct)
    df_onek1k_ct[, Chr := sub("chr", "", chr_b38)]
    df_onek1k_ct[, BP := as.character(pos_b38)]
    df_onek1k_ct <- df_onek1k_ct[, .(phenotype_id, Chr, BP, A1, A2, af, beta, pval_nominal, fdr)]
    df_onek1k_ct$n <- 982

    # --- Load Indep1K (our) data ---
    df_indep1k <- fread(file.path(path_to_sceqtl_data, "sig", paste0(ct_indep1k, "_EUR_pc5_sub_sig.tsv.gz")))
    df_indep1k[, n := round(ma_count / ifelse(af > 0.5, 1 - af, af) / 2)]
    df_indep1k[, c("Chr", "BP", "A2", "A1") := tstrsplit(variant_id, "_", keep = 1:4)]

    # --- Match datasets ---
    exact_match <- merge(df_indep1k, df_onek1k_ct, by = c("phenotype_id", "Chr", "BP", "A1", "A2"), all.x = FALSE, suffixes = c("_df", "_1k1k"))
    flip_match <- merge(df_indep1k, df_onek1k_ct, by.x = c("phenotype_id", "Chr", "BP", "A1", "A2"), by.y = c("phenotype_id", "Chr", "BP", "A2", "A1"), all.x = FALSE, suffixes = c("_df", "_1k1k"))
    flip_match[, beta := -beta]
    
    df_comb <- rbind(exact_match, flip_match)
    
    # --- Standardize effect sizes ---
    df_comb[, c("slope_std", "slope_std_se") := {
      res <- calcu_std_b_se(slope / slope_se, af_df, n_df)
      list(res$std_b_hat, res$std_se)
    }]
    df_comb[, c("betae_std", "beta_std_se") := {
      # Recalculate Z-score from p-value for OneK1K as SE is not provided
      z <- sqrt(qchisq(pval_nominal_1k1k, 1, lower.tail = FALSE)) * sign(beta)
      res <- calcu_std_b_se(z, af_1k1k, n_1k1k)
      list(res$std_b_hat, res$std_se)
    }]
    
    df_comb[, celltype := ct_indep1k]
    df_comb_all_list[[length(df_comb_all_list) + 1]] <- df_comb
    gc() # Clean up memory
  }
  saveRDS(df_comb_all_list, matched_rds_file)
} else {
  print("Loading pre-computed matched data from file.")
  df_comb_all_list <- readRDS(matched_rds_file)
}


## 4. PLOTTING AND SUMMARIZATION
# ------------------------------------------------------------------------------
print("--- Section 4: Generating Plots and Summaries ---")

# Combine list into a single data.table for plotting
df_celltype_all <- rbindlist(df_comb_all_list) %>%
  mutate(sig = ifelse(fdr < 0.05, "Significant", "Non-significant"))

# Set factor levels for consistent plotting order
ct_order <- cell_type_short_dict[celltype_indep1k_list]
df_celltype_all$cell_type_short <- factor(cell_type_short_dict[df_celltype_all$celltype], levels = ct_order)
color_pal <- c(cell_type_color_dict[levels(df_celltype_all$cell_type_short)], grey = "#bebebe83")

# --- 4.1: Effect Size Scatter Plot ---
p_scatter <- df_celltype_all %>%
  mutate(cell_type_color = ifelse(sig == "Significant", as.character(cell_type_short), "grey")) %>%
  ggplot(aes(x = slope_std, y = betae_std, color = cell_type_color)) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(values = color_pal) +
  facet_wrap(~cell_type_short, ncol = 4) +
  labs(
    x = "Standardized effect size of sc-eQTLs in Indep1K",
    y = "Standardized effect size in OneK1K",
    color = ""
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 20, family = "Helvetica"),
    axis.text = element_text(size = 20, color = 'black'),
    axis.title = element_text(size = 20, color = 'black'),
    panel.grid = element_blank(),
    strip.background = element_rect(colour = "white"),
    legend.position = 'none',
    strip.text = element_text(size = 20)
  )
ggsave(file.path(output_dir, "scatter_indep1k_vs_onek1k.png"), plot = p_scatter, width = 12, height = 10)

# --- 4.2: Sign Concordance Analysis ---
# Concordance for all matched significant sc-eQTLs
sign_df <- df_celltype_all %>%
  rowwise() %>%
  mutate(sign_concordance = (sign(slope_std) == sign(beta))) %>%
  ungroup() %>%
  group_by(cell_type_short) %>%
  summarise(Concordance = mean(sign_concordance, na.rm = TRUE), Disconcordance = 1 - mean(sign_concordance, na.rm = TRUE), n = n())

p_concordance <- sign_df %>%
  pivot_longer(cols = c(Concordance, Disconcordance), names_to = "Concordance_type", values_to = "value") %>%
  mutate(Concordance_type = factor(Concordance_type, levels = c("Disconcordance", "Concordance"))) %>%
  ggplot(aes(x = cell_type_short, y = value * 100, fill = interaction(cell_type_short, Concordance_type))) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = c(
    setNames(rep("#d1d1d1", length(unique(sign_df$cell_type_short))), paste0(unique(sign_df$cell_type_short), ".Disconcordance")),
    setNames(cell_type_color_dict[unique(sign_df$cell_type_short)], paste0(unique(sign_df$cell_type_short), ".Concordance"))
  )) +
  labs(fill = "", y = "Percentage of sign concordance (%)", x = "") +
  theme_classic() +
  theme(
    text = element_text(size = 20, family = "Helvetica"), 
    legend.position = "none",
    axis.text = element_text(size = 18, color = "black"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title = element_text(size = 18, color = "black")
  )
ggsave(file.path(output_dir, "sign_concordance_indep1k_vs_onek1k.png"), plot = p_concordance, width = 8, height = 6)

# Calculate and print concordance for eQTLs significant in both studies
sign_df2 <- df_celltype_all %>%
  filter(sig == "Significant") %>%
  rowwise() %>%
  mutate(sign_concordance = (sign(slope_std) == sign(beta))) %>%
  ungroup() %>%
  group_by(cell_type_short) %>%
  summarise(Concordance = mean(sign_concordance, na.rm = TRUE), n = n())

print("Sign concordance for eQTLs significant in both studies:")
print(as.data.frame(sign_df2))
