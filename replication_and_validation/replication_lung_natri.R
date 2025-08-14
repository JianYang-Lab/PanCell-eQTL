################################################################################
#
# Script: Replication and Comparison of Lung sc-eQTLs with Natri et al. Data
#
# Description:
# This script performs a replication analysis for single-cell eQTLs (sc-eQTLs)
# discovered in East Asian (EAS) Lung samples. It compares the sc-eQTL
# results against the relevant EUR lung sc-eQTL dataset from Natri et al.
#
# The script performs three main analyses:
#   1. Compares the standardized effect sizes of matched eQTLs.
#   2. Calculates and summarizes the sign concordance of effects.
#   3. Compares the number of eGenes and eQTLs discovered in both studies.
#
################################################################################


## 1. SETUP & GLOBAL PARAMETERS
# ------------------------------------------------------------------------------
print("--- Section 1: Setup & Global Parameters ---")

# Load required libraries
library(tidyverse)
library(data.table)
library(cowplot)
library(parallel)

# --- User-Defined Paths and Parameters ---

# Set the base directory for the project. All paths will be relative to this.
base_dir <- "/path/to/your/scRNA_project_folder/"
setwd(base_dir)

# Define input/output directory for plots and data files
output_dir <- "replication_plots/lung_natri/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Define path to the external Natri et al. Lung dataset
path_to_natri_data <- file.path(base_dir, "public_data/Natri_Lung/")

# Define parameters for the analysis
TISSUE <- "Lung"
ANCESTRY <- "EAS"
sceqtl_path <- file.path(base_dir, TISSUE, "sc-eQTL/results/sig/")


# --- Load Prerequisite Metadata ---
# This metadata is required for filtering, ordering, and coloring plots.
cell_type_cat_colors <- read_delim("metadata/cell_type_colors_new.txt", delim = '\t')

# Create dictionaries for cell type names and colors
cell_type_short_dict <- setNames(
    cell_type_cat_colors$cell_type_short[cell_type_cat_colors$tissue == TISSUE],
    cell_type_cat_colors$cell_type[cell_type_cat_colors$tissue == TISSUE]
)
cell_type_color_dict <- cell_type_cat_colors$cell_type_colors
names(cell_type_color_dict) <- cell_type_cat_colors$cell_type_short

# Define the list of lung cell types and the mapping to Natri et al. cell types
lung_cell_types <- c("CD4T","CD8T","NK","B","Monocytes","Mast","Macrophages","DC","AT1","AT2","Airway_Epi_Multiciliated","Airway_Epi_Basal","Airway_Epi_Secretory","EC_arterial","EC_capillary","EC_venous","LEC_mature","LEC_diff","Fibroblasts","Smooth_muscle","Mesothelium")

natri_celltypes_select <- c(
    EC_arterial = "endothelial_arteriole", LEC_diff = "endothelial_Lymphatic", EC_venous = "endothelial_venule",
    LEC_mature = "endothelial_Lymphatic", EC_capillary = "endothelial_gCap", NK = "immune_NK",
    Monocytes = "immune_Monocyte", Mast = "immune_Mast", Macrophages = "immune_Monocyte-derivedmacrophage",
    DC = "immune_moDC", CD4T = "immune_CD4", B = "immune_Bcells", CD8T = "immune_CD8NKT",
    AT1 = "epithelial_AT1", AT2 = "epithelial_AT2", Airway_Epi_Basal = "epithelial_Basal",
    Airway_Epi_Multiciliated = "epithelial_Ciliated", Airway_Epi_Secretory = "epithelial_Secretory-SCGB1A1+SCGB3A2+",
    Fibroblasts = "mesenchymal_AdventitialFB", Smooth_muscle = "mesenchymal_SMC", Mesothelium = "mesenchymal_Mesothelial"
)


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

#' Load and preprocess data from Natri et al. for a specific cell type.
load_natri <- function(natri_path, natri_ct) {
  print(paste0("Loading Natri et al. data for: ", natri_ct))
  file_path <- file.path(natri_path, paste0(natri_ct, "_qtl_results_all.txt"))
  if (!file.exists(file_path)) {
      warning(paste("File not found:", file_path)); return(NULL)
  }
  df1 <- fread(file_path)
  df1[, Chr := paste0("chr", feature_chromosome)]
  setnames(df1, c("feature_id", "snp_position", "assessed_allele", "beta", "beta_se"), c("phenotype_id", "BP", "A1", "slope", "slope_se"))
  df1[, BP := as.character(BP)]
  df1 <- df1[, .(phenotype_id, Chr, BP, A1, maf, slope, slope_se, p_value, n_e_samples)]
  return(df1)
}

#' Match sc-eQTLs with Natri eQTLs and standardize effect sizes.
match_df_natri <- function(sceqtl_path, celltype, celltype1, bulk_eqtl) {
  print(paste0("Matching: ", celltype, " (our study) vs. ", celltype1, " (Natri)"))
  sc_file <- file.path(sceqtl_path, paste0(celltype, "_", ANCESTRY, "_pc5_sig_5en8.tsv.gz"))
  if (!file.exists(sc_file)) {
      warning(paste("sc-eQTL file not found:", sc_file)); return(NULL)
  }
  df <- fread(sc_file, select = c(1, 2, 4, 6:9))
  df[, c("Chr", "BP", "A2", "A1") := tstrsplit(variant_id, "_", keep = 1:4)]
  df[, Chr := paste0("chr", Chr)]
  df[, maf_df := ifelse(af > 0.5, 1 - af, af)]

  exact_match <- merge(df, bulk_eqtl, by = c("phenotype_id", "Chr", "BP", "A1"), suffixes = c("_df", "_bulk"))
  flip_match <- merge(df, bulk_eqtl, by.x = c("phenotype_id", "Chr", "BP", "A2"), by.y = c("phenotype_id", "Chr", "BP", "A1"), suffixes = c("_df", "_bulk"))
  flip_match[, slope_bulk := -slope_bulk]

  df_comb <- rbind(exact_match, flip_match)
  df_comb[, c("slope_std", "slope_std_se") := {
    res <- calcu_std_b_se(slope_df / slope_se_df, maf_df, round(ma_count / maf_df / 2))
    list(res$std_b_hat, res$std_se)
  }]
  df_comb[, c("beta_std", "beta_std_se") := {
    res <- calcu_std_b_se(slope_bulk / slope_se_bulk, maf_bulk, n_e_samples)
    list(res$std_b_hat, res$std_se)
  }]
  df_comb[, celltype := celltype]
  df_comb[, celltype1 := celltype1]
  return(df_comb)
}


## 3. ANALYSIS 1: EFFECT SIZE COMPARISON
# ------------------------------------------------------------------------------
print("--- Section 3: Analysis 1 - Effect Size Comparison ---")
matched_rds_file <- file.path(output_dir, "matched_sceqtl_natri_Lung_EAS.rds")

if (!file.exists(matched_rds_file)) {
    df_comb_list <- list()
    celltypes_to_process <- intersect(names(natri_celltypes_select), lung_cell_types)
    
    for (ct in celltypes_to_process) {
        natri_ct <- natri_celltypes_select[ct]
        bulk_eqtl <- load_natri(path_to_natri_data, natri_ct)
        if (!is.null(bulk_eqtl)) {
            matched_data <- match_df_natri(sceqtl_path, ct, natri_ct, bulk_eqtl)
            if (!is.null(matched_data) && nrow(matched_data) > 0) {
                 df_comb_list[[length(df_comb_list) + 1]] <- matched_data
            }
        }
    }
    saveRDS(df_comb_list, matched_rds_file)
} else {
    print("Loading pre-computed matched data from file.")
    df_comb_list <- readRDS(matched_rds_file)
}

# --- Plotting Effect Size Comparison ---
df_comb_all <- rbindlist(df_comb_list) %>%
  filter(pval_nominal < 5e-8) %>%
  mutate(sig = ifelse(p_value < 5e-8, "Significant", "Non-significant")) %>%
  left_join(cell_type_cat_colors %>% filter(tissue == TISSUE) %>% rename(celltype = cell_type) %>% select(celltype, cell_type_short), by = "celltype")

df_comb_all$cell_type_short <- factor(df_comb_all$cell_type_short, levels = cell_type_short_dict)
color_pal <- c(cell_type_color_dict[levels(df_comb_all$cell_type_short)], grey = "#bebebe83")

p_scatter <- df_comb_all %>%
  mutate(cell_type_color = ifelse(sig == "Significant", as.character(cell_type_short), "grey")) %>%
  ggplot(aes(x = slope_std, y = beta_std, color = cell_type_color)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(values = color_pal) +
  facet_wrap(~cell_type_short, ncol = 4) +
  labs(x = "Standardized effect size in our study (EAS)", y = "Standardized effect size in Natri et al. (EUR)", color = "") +
  theme_classic() +
  theme(text = element_text(size = 30), panel.grid = element_blank(), strip.background = element_rect(colour = "white"), legend.position = "none")

ggsave(file.path(output_dir, "scatter_Lung_EAS_vs_Natri.png"), plot = p_scatter, width = 20, height = 18)

# Calculate and print sign concordance
sign_df <- df_comb_all %>%
  filter(!is.na(beta_std) & !is.na(slope_std)) %>%
  rowwise() %>%
  mutate(sign_concordance = (sign(slope_std) == sign(beta_std))) %>%
  ungroup() %>%
  group_by(cell_type_short) %>%
  summarise(Concordance = mean(sign_concordance, na.rm = TRUE), n = n())

print("Sign concordance for all significant sc-eQTLs:")
print(as.data.frame(sign_df))


## 4. ANALYSIS 2: EGENE AND EQTL COUNT COMPARISON
# ------------------------------------------------------------------------------
print("--- Section 4: Analysis 2 - eGene and eQTL Count Comparison ---")
stats_rds_file <- file.path(output_dir, "stats_eGene_eQTL_counts_Lung_EAS.rds")

if (!file.exists(stats_rds_file)) {
    stats_list <- mclapply(names(natri_celltypes_select), function(ct) {
        natri_ct <- natri_celltypes_select[ct]
        bulk_eqtl <- load_natri(path_to_natri_data, natri_ct)
        if (is.null(bulk_eqtl)) return(NULL)

        sc_file <- file.path(sceqtl_path, paste0(ct, "_", ANCESTRY, "_pc5_sig_5en8.tsv.gz"))
         if (!file.exists(sc_file)) return(NULL)
        sc_eqtl <- fread(sc_file)

        negene_natri <- length(unique(bulk_eqtl[p_value < 5e-8, ]$phenotype_id))
        neqtl_natri <- nrow(bulk_eqtl[p_value < 5e-8, ])
        negene_df <- length(unique(sc_eqtl[pval_nominal < 5e-8, ]$phenotype_id))
        neqtl_df <- nrow(sc_eqtl[pval_nominal < 5e-8, ])

        data.frame(celltype_df = ct, celltype_natri = natri_ct, negene_natri, neqtl_natri, negene_df, neqtl_df)
    }, mc.cores = detectCores() - 1)

    stats <- rbindlist(Filter(Negate(is.null), stats_list))
    saveRDS(stats, stats_rds_file)
} else {
    print("Loading pre-computed eGene/eQTL counts from file.")
    stats <- readRDS(stats_rds_file)
}

# --- Plotting eGene and eQTL Counts ---
stats$cell_type_short <- factor(cell_type_short_dict[stats$celltype_df], levels = cell_type_short_dict)

p1 <- stats %>%
  ggplot(aes(x = negene_natri, y = negene_df, color = cell_type_short)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_color_manual(values = cell_type_color_dict[levels(stats$cell_type_short)]) +
  labs(x = "No. eGenes in Natri et al. (EUR)", y = "No. eGenes in our study (EAS)", color = "Cell type") +
  geom_point(size = 5) +
  theme_classic() +
  theme(text = element_text(size = 20, family = "Helvetica"), axis.text = element_text(size = 18, color = "black"), axis.title = element_text(size = 18, color = "black"))

p2 <- stats %>%
  ggplot(aes(x = neqtl_natri, y = neqtl_df, color = cell_type_short)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_color_manual(values = cell_type_color_dict[levels(stats$cell_type_short)]) +
  labs(x = "No. eQTLs in Natri et al. (EUR)", y = "No. eQTLs in our study (EAS)", color = "Cell type") +
  scale_y_continuous(labels = scales::label_comma()) +
  geom_point(size = 5) +
  theme_classic() +
  theme(text = element_text(size = 20, family = "Helvetica"), axis.text = element_text(size = 18, color = "black"), axis.title = element_text(size = 18, color = "black"))

# Save individual plots
ggsave(file.path(output_dir, "count_comparison_eGenes.png"), plot = p1, width = 7, height = 6)
ggsave(file.path(output_dir, "count_comparison_eQTLs.png"), plot = p2, width = 7, height = 6)

# Save combined plot
p_combined <- plot_grid(
    p1 + theme(legend.position = "none", plot.margin = margin(0.5, 0.5, 0.5, 0.5, unit = "cm")),
    p2 + theme(legend.position = "none", plot.margin = margin(0.5, 0.5, 0.5, 0.5, unit = "cm")),
    nrow = 2, axis = "lr", align = "v"
)
ggsave(file.path(output_dir, "count_comparison_combined.png"), plot = p_combined, width = 5, height = 8)
