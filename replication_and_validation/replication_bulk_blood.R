################################################################################
#
# Script: Replication and Comparison of Blood sc-eQTLs
#
# Description:
# This script performs a series of replication and comparison analyses for
# single-cell eQTLs discovered in Blood. It compares the sc-eQTL results
# against several large-scale bulk eQTL datasets from different ancestries:
#   1. eQTLGen (primarily European ancestry)
#   2. JCTF (Japanese, East Asian ancestry)
#   3. AFGR (African ancestry)
#   4. MAGE (Multi-ancestry)
#
# The script loads data, matches SNP-gene pairs, handles allele flipping,
# standardizes effect sizes, and generates summary plots.
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
library(ggvenn)
library(cowplot)
library(arrow)

# Set global ggplot theme
theme_set(
  theme_minimal() +
    theme(text = element_text(family = "Helvetica"))
)

# --- User-Defined Paths and Parameters ---

# Set the base directory for the project
base_dir <- "/path/to/scRNA/"
setwd(base_dir)

# Define output directory for plots and RDS files
output_dir <- "replication_plots/blood/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Define paths to external datasets and resources
# Users should modify these paths to match their environment.
path_to_gtex <- "/path/to/eQTL/gtex_resources_besd/eQTL_raw"
path_to_eqtlgen <- "/path/to/eQTL/eQTLGen/cis-eQTL/FULL/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz"
path_to_jctf <- "/path/to/eQTL/JCTF_Japanese_blood_eQTL/raw/"
path_to_afgr <- "/path/to/eQTL/AFGR/"
path_to_mage <- "/path/to/eQTL/MAGE/QTL_results/eQTL_results/eQTL_nominal_results/"
path_to_liftover_chain <- "/path/to/liftover/b37ToHg38.over.chain"
path_to_gene_gtf <- file.path(base_dir, "Blood/resource/genes.gtf")

# --- Global Variables & Metadata Loading ---

tissue_list <- c("Blood", "Lung", "Skin", "Colon", "Liver")
ancestry_list <- c("EUR", "EAS", "AFR", "AMR")
ancestry_color_list <- c(EUR = "#66C2A5", EAS = "#FC8D62", AFR = "#8DA0CB", AMR = "#E78AC3")

celltypes_all_list <- list(
  Blood = c('CD4TNC', 'CD4TEM', 'Treg', 'CD8TNC', 'CD8TEM', 'CD8TEMRA', 'MAIT', 'NKp', 'NKn', 'BIN', 'BMem', 'Plasma', 'Plasmablasts', 'MonoC', 'MonoNC', 'Nph', 'DC', 'pDC')
  # Other tissues omitted for this script's focus
)

# Load cell type color and category information
cell_type_cat_colors <- read_delim("metadata/cell_type_colors_new.txt", delim = '\t') %>%
  mutate(lineage = case_when(
    cell_type_category == "T cells" ~ "Immune Lymphoid",
    cell_type_category == "B cells" ~ "Immune Lymphoid",
    cell_type_category == "ILC" ~ "Immune Lymphoid",
    cell_type_category == "Myeloid" ~ "Immune Myeloid",
    cell_type_category == "Endothelial" ~ "Endothelial",
    cell_type_category == "Epithelial" ~ "Epithelial",
    cell_type_category == "Fibroblasts/muscle cells" ~ "Mesenchymal",
    cell_type == "Mesothelium" ~ "Mesenchymal",
    .default = "Other"
  ))

# Create dictionaries for cell type names and colors
cell_type_short_dict <- lapply(tissue_list, function(x) {
  tmp <- cell_type_cat_colors$cell_type_short[cell_type_cat_colors$tissue == x]
  names(tmp) <- cell_type_cat_colors$cell_type[cell_type_cat_colors$tissue == x]
  tmp
})
names(cell_type_short_dict) <- tissue_list
cell_type_color_dict <- cell_type_cat_colors$cell_type_colors
names(cell_type_color_dict) <- cell_type_cat_colors$cell_type_short
lineage_colors <- c(`Immune Lymphoid` = "#8AB9D8", `Immune Myeloid` = "#e5a03d", Epithelial = "#8d4c28", Endothelial = "#b84e3c", Mesenchymal = "#c76b85", Other = "#c3e0e5")

# Load and process summary statistics to define high-quality conditions
smr_list <- list()
for (t in tissue_list) {
  smr_file <- file.path(t, "sc-eQTL", "results", "sig", "eQTL_summary_5en8_maf01.tsv")
  if (file.exists(smr_file)) {
    smr <- read_table(smr_file)
    smr$tissue <- t
    smr_list[[t]] <- smr
  }
}
ts_ct_anc <- read_delim("metadata/tissue_celltype_ancestry_ss_new.txt", delim = '\t')

smr_list_all <- list_rbind(smr_list) %>%
  dplyr::filter(!(tissue == "Blood" & cell_type %in% c("T", "B"))) %>%
  left_join(ts_ct_anc, by = c("tissue", "cell_type", "ancestry")) %>%
  dplyr::filter(sample_size >= 50) %>%
  inner_join(cell_type_cat_colors, by = c("tissue", "cell_type")) %>%
  mutate(
    ancestry = factor(ancestry, levels = ancestry_list),
    ct_anc = paste0(cell_type_short, "-", ancestry),
    total_count = round(sample_size * mean_count)
  ) %>%
  dplyr::filter(total_count > 10000) %>%
  dplyr::rename(neGenes = neGenes_5en8, neQTLs = neQTLs_5en8) %>%
  distinct()


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
  snp_set <- distinct(df[, .(chr_b37, pos_b37)])
  snp_gr <- GRanges(
    seqnames = snp_set$chr_b37,
    ranges = IRanges(start = snp_set$pos_b37, end = snp_set$pos_b37)
  )
  lifted_snps <- liftOver(snp_gr, chain_file)
  lifted_positions <- as.data.table(lifted_snps)
  old_positions <- as.data.table(snp_gr)
  old_positions[, group := .I]
  lifted_positions <- merge(old_positions[, .(chr_b37 = seqnames, pos_b37 = start, group)], lifted_positions[, .(chr_b38 = seqnames, pos_b38 = start, group)], by = "group")
  lifted_positions[, chr_b38 := as.character(chr_b38)]
  lifted_positions[, chr_b37 := as.integer(as.character(sub("chr", "", chr_b37)))]
  df <- lifted_positions[, .(chr_b37, pos_b37, chr_b38, pos_b38)][df, on = c("chr_b37", "pos_b37")]
  return(df)
}


## 3. ANALYSIS: BLOOD (EUR) vs. eQTLGEN
# ------------------------------------------------------------------------------
print("--- Section 3: Analysis - Blood (EUR) vs. eQTLGen ---")

# --- 3.1: eGene Overlap and Expression Analysis ---
print("Step 3.1: Comparing eGene overlap with eQTLGen.")
sceqtl_path_blood <- file.path(base_dir, "Blood/sc-eQTL/results/")
eGene_rds_file <- file.path(output_dir, "Blood_eQTLGen_eGene_processed.rds")

if (!file.exists(eGene_rds_file)) {
  EUR_bulk_eqtl_eQTLGen <- fread(path_to_eqtlgen)
  egene_bulk <- unique(EUR_bulk_eqtl_eQTLGen[FDR < 0.05]$Gene)

  egene_df_list <- list()
  for (ct in celltypes_all_list$Blood) {
    perm_file <- file.path(sceqtl_path_blood, "raw", paste0(ct, "_EUR_pc5_perm.tsv"))
    if (file.exists(perm_file)) {
      df <- fread(perm_file, select = c("phenotype_id", "pval_nominal", "pval_nominal_threshold", "af", "qval"))
      sig_df <- unique(df[pval_nominal < pval_nominal_threshold & af >= 0.01 & af <= 0.99 & qval <= 0.05, ]$phenotype_id)
      egene_df_list[[ct]] <- sig_df
    }
  }
  egene_df_all <- unique(unlist(egene_df_list))

  gene_expr_df <- fread(file.path(base_dir, "expression_decile/Blood_EUR_cross_celltype.txt"))
  gene_gtf <- fread(path_to_gene_gtf)
  setnames(gene_gtf, "gene_id", "phenotype_id")
  gene_expr_df <- gene_gtf[, .(phenotype_id, gene_name)][gene_expr_df, on = 'gene_name']
  gene_expr_df <- gene_expr_df[, -13]

  gene_expr_df_count <- data.frame(phenotype_id = gene_expr_df$phenotype_id, ct_count = apply(gene_expr_df[, -c(1:2)], 1, function(x) sum(!is.na(x))))
  gene_expr_df_max <- data.table(phenotype_id = gene_expr_df$phenotype_id, expr_max = apply(gene_expr_df[, -c(1:2)], 1, function(x) max(x, na.rm = TRUE)))
  gene_expr_df_max$ct_count <- gene_expr_df_count$ct_count

  egene_df_all_sig_in_bulk <- data.table(phenotype_id = egene_df_all, egene_df_all_sig_in_bulk = (egene_df_all %in% egene_bulk))

  out <- egene_df_all_sig_in_bulk[gene_expr_df_max, on = "phenotype_id"]
  out <- out[!is.na(egene_df_all_sig_in_bulk)]
  out[, expr_decile := cut(expr_max, breaks = quantile(out$expr_max, probs = seq(0, 1, 0.1)), labels = 1:10, include.lowest = TRUE)]
  out[is.na(expr_decile), "expr_decile"] <- 1
  saveRDS(out, eGene_rds_file)
} else {
  out <- readRDS(eGene_rds_file)
}

# Plot 1: eGene count by expression decile
p1 <- out %>%
  group_by(egene_df_all_sig_in_bulk, expr_decile) %>%
  summarise(count = n(), .groups = 'drop') %>%
  mutate(
    group = ifelse(egene_df_all_sig_in_bulk, "eGene in eQTLGen", "Not an eGene in eQTLGen"),
    expr_decile = factor(expr_decile, labels = paste0(seq(10, 100, 10), "%"))
  ) %>%
  ggplot(aes(x = expr_decile, y = count, fill = group)) +
  geom_bar(stat = 'identity', position = 'dodge', width = 0.8) +
  theme_classic() +
  scale_y_continuous(labels = scales::label_comma()) +
  scale_fill_manual(values = c("#f1bd69", "#2A629A")) +
  theme(
    text = element_text(size = 18, family = "Helvetica"),
    axis.text = element_text(size = 18, color = "black"),
    axis.title = element_text(size = 18, color = "black"),
    legend.position = "none"
  ) +
  labs(fill = "", x = "Max expression decile", y = "eGene (n)")
ggsave(file.path(output_dir, "Blood_eQTLGen_eGene_by_expression.png"), plot = p1, width = 8, height = 6)

# Plot 2: eGene count by number of cell types expressed
p2 <- out %>%
  group_by(egene_df_all_sig_in_bulk, ct_count) %>%
  summarise(count = n(), .groups = 'drop') %>%
  mutate(group = ifelse(egene_df_all_sig_in_bulk, "eGene in eQTLGen", "Not an eGene in eQTLGen")) %>%
  ggplot(aes(x = ct_count, y = count, fill = group)) +
  geom_bar(stat = 'identity', width = 0.8, position = 'dodge') +
  scale_x_continuous(breaks = 1:17) +
  scale_fill_manual(values = c("#f1bd69", "#2A629A")) +
  theme_classic() +
  scale_y_continuous(labels = scales::label_comma()) +
  theme(
    text = element_text(size = 18, family = "Helvetica"),
    axis.text = element_text(size = 18, color = "black"),
    axis.title = element_text(size = 18, color = "black"),
    legend.position = c(0.25, 0.8)
  ) +
  labs(fill = "", x = "Number of cell types gene is tested in", y = "eGene (n)")
ggsave(file.path(output_dir, "Blood_eQTLGen_eGene_by_celltype_count.png"), plot = p2, width = 8, height = 6)

# Combined plot
plot_grid(p1, p2, nrow = 2)
ggsave(file.path(output_dir, "Blood_eQTLGen_eGene_stats_combined.png"), width = 8, height = 8)


# --- 3.2: Effect Size Comparison ---
print("Step 3.2: Comparing effect sizes with eQTLGen.")
matched_rds_file <- file.path(output_dir, "matched_sceqtl_bulk_Blood_EUR.rds")

if (!file.exists(matched_rds_file)) {
  # Load and prepare eQTLGen data (requires liftover)
  EUR_bulk_eqtl_eQTLGen <- fread(path_to_eqtlgen)
  setnames(EUR_bulk_eqtl_eQTLGen, c("Gene", "SNPChr", "SNPPos", "AssessedAllele", "OtherAllele", "NrSamples", "Pvalue"), c("phenotype_id", "Chr", "BP", "A1", "A2", "N", "pval_nominal"))
  EUR_bulk_eqtl_eQTLGen <- liftover(EUR_bulk_eqtl_eQTLGen)
  EUR_bulk_eqtl_eQTLGen[, Chr := sub("chr", "", chr_b38)]
  EUR_bulk_eqtl_eQTLGen[, BP := as.character(pos_b38)]

  df_comb_all_list <- list()
  for (i in 1:length(celltypes_all_list$Blood)) {
    ct <- celltypes_all_list$Blood[i]
    sig_file <- file.path(sceqtl_path_blood, "sig", paste0(ct, "_EUR_pc5_sig.tsv.gz"))
    if (file.exists(sig_file)) {
      print(paste("... matching", ct))
      df <- fread(sig_file)
      df[, n := round(ma_count / ifelse(af > 0.5, 1 - af, af) / 2)]
      df[, c("Chr", "BP", "A2", "A1") := tstrsplit(variant_id, "_", keep = 1:4)]
      
      exact_match <- merge(df, EUR_bulk_eqtl_eQTLGen, by = c("phenotype_id", "Chr", "BP", "A1", "A2"), all.x = FALSE, suffixes = c("_df", "_bulk"))
      flip_match <- merge(df, EUR_bulk_eqtl_eQTLGen, by.x = c("phenotype_id", "Chr", "BP", "A1", "A2"), by.y = c("phenotype_id", "Chr", "BP", "A2", "A1"), all.x = FALSE, suffixes = c("_df", "_bulk"))
      flip_match[, Zscore := -Zscore]
      
      df_comb <- rbind(exact_match, flip_match)
      df_comb[, c("slope_std", "slope_std_se") := { res <- calcu_std_b_se(slope / slope_se, af, n); list(res$std_b_hat, res$std_se) }]
      df_comb[, c("beta_std", "beta_std_se") := { res <- calcu_std_b_se(Zscore, af, N); list(res$std_b_hat, res$std_se) }]
      df_comb[, celltype := ct]
      df_comb_all_list[[length(df_comb_all_list) + 1]] <- df_comb
    }
  }
  saveRDS(df_comb_all_list, matched_rds_file)
} else {
  df_comb_all_list <- readRDS(matched_rds_file)
}

# Plotting effect size correlations
df_celltype_all <- rbindlist(df_comb_all_list) %>%
  mutate(sig = ifelse(FDR < 0.05, "Significant", "Non-significant"))
ct_order <- cell_type_short_dict$Blood[celltypes_all_list$Blood]
df_celltype_all$cell_type_short <- factor(df_celltype_all$celltype, levels = ct_order)
color_pal <- c(cell_type_color_dict[levels(df_celltype_all$cell_type_short)], grey = "#bebebe83")

df_celltype_all %>%
  filter(cell_type_short != "Plasmablasts") %>%
  mutate(cell_type_color = ifelse(sig == "Significant", as.character(cell_type_short), "grey")) %>%
  ggplot(aes(x = slope_std, y = beta_std, color = cell_type_color)) +
  geom_point(size = 3) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  scale_color_manual(values = color_pal) +
  facet_wrap(~cell_type_short, ncol = 4) +
  labs(x = "Standardized effect size in our study", y = "Standardized effect size in eQTLGen", color = "") +
  theme_classic() +
  theme(text = element_text(size = 24), panel.grid = element_blank(), strip.background = element_rect(colour = "white"), legend.position = 'none')
ggsave(file.path(output_dir, "Blood_EUR_eQTLGen_effect_size_scatter.png"), width = 20, height = 18)

# Plotting sign concordance
sign_df <- df_celltype_all %>%
  rowwise() %>%
  mutate(sign_concordance = (sign(slope_std) == sign(beta_std))) %>%
  ungroup() %>%
  group_by(cell_type_short) %>%
  summarise(Concordance = mean(sign_concordance, na.rm = TRUE), Disconcordance = 1 - mean(sign_concordance, na.rm = TRUE), n = n())

sign_df %>%
  filter(cell_type_short != "Plasmablasts") %>%
  pivot_longer(cols = c(Concordance, Disconcordance), names_to = "Concordance_type", values_to = "value") %>%
  mutate(Concordance_type = factor(Concordance_type, levels = c("Disconcordance", "Concordance"))) %>%
  ggplot(aes(x = value * 100, y = fct_rev(cell_type_short), fill = interaction(cell_type_short, Concordance_type))) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = c(
    setNames(rep("#d1d1d1", length(unique(sign_df$cell_type_short))), paste0(unique(sign_df$cell_type_short), ".Disconcordance")),
    setNames(cell_type_color_dict[unique(sign_df$cell_type_short)], paste0(unique(sign_df$cell_type_short), ".Concordance"))
  )) +
  labs(fill = "", y = "Cell type", x = "Percentage of sign concordance (%)") +
  theme_classic() +
  theme(text = element_text(size = 18), legend.position = "none")
ggsave(file.path(output_dir, "Blood_EUR_eQTLGen_sign_concordance.png"), width = 6, height = 10)


## 4. ANALYSIS: BLOOD (EAS) vs. JCTF
# ------------------------------------------------------------------------------
print("--- Section 4: Analysis - Blood (EAS) vs. JCTF ---")

matched_rds_file_eas <- file.path(output_dir, "matched_sceqtl_bulk_Blood_EAS.rds")
if(!file.exists(matched_rds_file_eas)) {
  EAS_bulk_eqtl <- fread(gzfile(file.path(path_to_jctf, "eqtl_sumstats_wang_qs_et_al_2024_japan_covid19_taskforce.tsv.gz")))
  EAS_bulk_eqtl <- EAS_bulk_eqtl[!grepl("chrX|chrY", variant_id_hg38)]
  EAS_bulk_eqtl[, phenotype_id := tstrsplit(gene_id, "[.]", keep = 1)]
  EAS_bulk_eqtl[, c("chr", "snp_pos", "ref", "alt") := tstrsplit(variant_id_hg38, "[:]")]
  EAS_bulk_eqtl[, chr := sub("chr", "", chr)]

  df_celltype_list <- list()
  for (celltype in celltypes_all_list$Blood) {
    sig_file <- file.path(sceqtl_path_blood, "sig", paste0(celltype, "_EAS_pc5_sig.tsv.gz"))
    if (file.exists(sig_file)) {
      print(paste("... matching", celltype))
      df <- fread(sig_file)
      df[, c("chr", "snp_pos", "ref", "alt") := tstrsplit(variant_id, "_", keep = 1:4)]
      df[, n := round(ma_count / ifelse(af > 0.5, 1 - af, af) / 2)]
      
      exact_match <- merge(df, EAS_bulk_eqtl, by = c("phenotype_id", "chr", "snp_pos", "ref", "alt"), all.x = FALSE, suffixes = c("_df", "_bulk"))
      flip_match <- merge(df, EAS_bulk_eqtl, by.x = c("phenotype_id", "chr", "snp_pos", "ref", "alt"), by.y = c("phenotype_id", "chr", "snp_pos", "alt", "ref"), all.x = FALSE, suffixes = c("_df", "_bulk"))
      flip_match[, slope_bulk := -slope_bulk]
      flip_match[, maf := 1 - maf]
      
      df_comb <- rbind(exact_match, flip_match)
      df_comb[, c("slope_std", "slope_std_se") := { res <- calcu_std_b_se(slope_df / slope_se_df, af, n); list(res$std_b_hat, res$std_se) }]
      df_comb[, c("beta_std", "beta_std_se") := { res <- calcu_std_b_se(slope_bulk / slope_se_bulk, maf, 1405); list(res$std_b_hat, res$std_se) }]
      df_comb[, ancestry := "EAS"]
      df_comb[, celltype := celltype]
      df_celltype_list[[paste0("EAS_", celltype)]] <- df_comb[!is.na(df_comb$slope_df), ]
    }
  }
  saveRDS(df_celltype_list, matched_rds_file_eas)
} else {
  df_celltype_list <- readRDS(matched_rds_file_eas)
}

# Plotting for EAS
df_celltype_all <- rbindlist(df_celltype_list) %>%
  filter(pval_nominal_df < 5e-8) %>%
  mutate(sig = ifelse(pval_nominal_bulk < 5e-8, "Significant", "Non-significant"))
df_celltype_all$cell_type_short <- factor(df_celltype_all$celltype, levels = ct_order)

df_celltype_all %>%
  filter(cell_type_short != "Plasmablasts") %>%
  mutate(cell_type_color = ifelse(sig == "Significant", as.character(cell_type_short), "grey")) %>%
  ggplot(aes(x = slope_std, y = beta_std, color = cell_type_color)) +
  geom_point(size = 3) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  scale_color_manual(values = color_pal) +
  facet_wrap(~cell_type_short, ncol = 4) +
  labs(x = "Standardized effect size in our study (EAS)", y = "Standardized effect size in JCTF", color = "") +
  theme_classic() + theme(text = element_text(size = 24), legend.position = 'none')
ggsave(file.path(output_dir, "Blood_EAS_JCTF_effect_size_scatter.png"), width = 20, height = 18)


## 5. ANALYSIS: BLOOD (AFR) vs. AFGR
# ------------------------------------------------------------------------------
print("--- Section 5: Analysis - Blood (AFR) vs. AFGR ---")

matched_rds_file_afr <- file.path(output_dir, "matched_sceqtl_bulk_Blood_AFR.rds")
if(!file.exists(matched_rds_file_afr)) {
  AFR_bulk_eqtl <- fread(file.path(path_to_afgr, "AFGR.sorted.dist.hwe.af.AFR_META.eQTL.nominal.hg38a.txt.gz"))
  AFR_bulk_eqtl[, phenotype_id := tstrsplit(feature, "[.]", keep = 1)]
  AFR_bulk_eqtl <- AFR_bulk_eqtl[, .(chr, snp_pos, ref, alt, effect_af_eqtl, phenotype_id, pvalue, beta, se)]

  df_celltype_list <- list()
  for (celltype in celltypes_all_list$Blood) {
    sig_file <- file.path(sceqtl_path_blood, "sig", paste0(celltype, "_AFR_pc5_sig.tsv.gz"))
    if (file.exists(sig_file)) {
        print(paste("... matching", celltype))
        df <- fread(sig_file)
        df[, c("chr", "snp_pos", "ref", "alt") := tstrsplit(variant_id, "_", keep = 1:4)]
        df[, n := round(ma_count / ifelse(af > 0.5, 1 - af, af) / 2)]
        df$chr <- as.numeric(df$chr)
        df$snp_pos <- as.numeric(df$snp_pos)

        exact_match <- merge(df, AFR_bulk_eqtl, by = c("phenotype_id", "chr", "snp_pos", "ref", "alt"), all.x = FALSE)
        flip_match <- merge(df, AFR_bulk_eqtl, by.x = c("phenotype_id", "chr", "snp_pos", "ref", "alt"), by.y = c("phenotype_id", "chr", "snp_pos", "alt", "ref"), all.x = FALSE)
        flip_match[, beta := -beta]
        
        df_comb <- rbind(exact_match, flip_match)
        df_comb[, c("slope_std", "slope_std_se") := { res <- calcu_std_b_se(slope / slope_se, af, n); list(res$std_b_hat, res$std_se) }]
        df_comb[, c("beta_std", "beta_std_se") := { res <- calcu_std_b_se(beta / se, ifelse(is.na(effect_af_eqtl), af, effect_af_eqtl), 599); list(res$std_b_hat, res$std_se) }]
        df_comb[, ancestry := "AFR"]
        df_comb[, celltype := celltype]
        df_celltype_list[[paste0("AFR_", celltype)]] <- df_comb[!is.na(df_comb$slope), ]
    }
  }
  saveRDS(df_celltype_list, matched_rds_file_afr)
}


## 6. ANALYSIS: BLOOD (AMR) vs. MAGE
# ------------------------------------------------------------------------------
print("--- Section 6: Analysis - Blood (AMR) vs. MAGE ---")

matched_rds_file_amr <- file.path(output_dir, "matched_sceqtl_bulk_Blood_AMR.rds")
if(!file.exists(matched_rds_file_amr)) {
  AMR_bulk_eqtl <- fread(file.path(path_to_mage, "eQTL_FastQTL_results.nominal_pass.allAssociations.MAGE.v1.0.txt.gz"))
  setnames(AMR_bulk_eqtl, c("variantChrom", "variantPosition", "variantRef", "variantAlt"), c("chr", "snp_pos", "ref", "alt"))
  AMR_bulk_eqtl <- AMR_bulk_eqtl[!grepl("chrX|chrY", chr)]
  AMR_bulk_eqtl[, phenotype_id := tstrsplit(ensemblID, "[.]", keep = 1)]
  AMR_bulk_eqtl[, chr := sub("chr", "", chr)]
  AMR_bulk_eqtl <- AMR_bulk_eqtl[, .(chr, snp_pos, ref, alt, phenotype_id, ma_samples, ma_count, maf, pval_nominal, slope, slope_se)]
  AMR_bulk_eqtl[, n := round(ma_count / maf / 2)]

  df_celltype_list <- list()
  for (celltype in celltypes_all_list$Blood) {
    sig_file <- file.path(sceqtl_path_blood, "sig", paste0(celltype, "_AMR_pc5_sig.tsv.gz"))
    if (file.exists(sig_file)) {
        print(paste("... matching", celltype))
        df <- fread(sig_file)
        df[, c("chr", "snp_pos", "ref", "alt") := tstrsplit(variant_id, "_", keep = 1:4)]
        df[, n := round(ma_count / ifelse(af > 0.5, 1 - af, af) / 2)]
        df$snp_pos <- as.numeric(df$snp_pos)
        
        exact_match <- merge(df, AMR_bulk_eqtl, by = c("phenotype_id", "chr", "snp_pos", "ref", "alt"), suffixes = c("_df", "_bulk"))
        flip_match <- merge(df, AMR_bulk_eqtl, by.x = c("phenotype_id", "chr", "snp_pos", "ref", "alt"), by.y = c("phenotype_id", "chr", "snp_pos", "alt", "ref"), suffixes = c("_df", "_bulk"))
        flip_match[, slope_df := -slope_df] # In AMR section, the flip was on slope_df, not slope_bulk
        
        df_comb <- rbind(exact_match, flip_match)
        df_comb[, c("slope_std", "slope_std_se") := { res <- calcu_std_b_se(slope_df / slope_se_df, af, n_df); list(res$std_b_hat, res$std_se) }]
        df_comb[, c("beta_std", "beta_std_se") := { res <- calcu_std_b_se(slope_bulk / slope_se_bulk, ifelse(is.na(maf), af, maf), n_bulk); list(res$std_b_hat, res$std_se) }]
        df_comb[, ancestry := "AMR"]
        df_comb[, celltype := celltype]
        df_celltype_list[[paste0("AMR_", celltype)]] <- df_comb[!is.na(df_comb$slope_df), ]
    }
  }
  saveRDS(df_celltype_list, matched_rds_file_amr)
}
