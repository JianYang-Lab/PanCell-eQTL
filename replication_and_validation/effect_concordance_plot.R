################################################################################
#
# Script: Generate Effect Size Concordance Plots
#
# Description:
# This script loads pre-matched single-cell and bulk eQTL summary statistics
# to assess the concordance of effect sizes between them. It performs the
# following steps:
#   1. Iterates through specified tissues and ancestries.
#   2. Loads the matched eQTL data for each combination.
#   3. Generates and saves a scatter plot comparing standardized effect sizes.
#   4. Calculates and saves the sign concordance for two sets of eQTLs:
#      a) All significant sc-eQTLs.
#      b) eQTLs that are significant in both sc-eQTL and bulk studies.
#   5. Aggregates the concordance results and generates a final summary plot
#      comparing concordance across all analyzed tissue-ancestry pairs.
#
# Original functionality from the source script is fully preserved.
#
################################################################################


## 1. SETUP & GLOBAL PARAMETERS
# ------------------------------------------------------------------------------
print("--- Section 1: Setup & Global Parameters ---")

# Load required libraries
library(tidyverse)
library(data.table)
library(cowplot)

# --- User-Defined Paths and Parameters ---

# Set the base directory for the project. All paths will be relative to this.
base_dir <- "/path/to/your/scRNA_project_folder/"
setwd(base_dir)

# Define input/output directory for plots and intermediate data
# This script assumes the matched RDS files from previous steps are here.
data_dir <- "replication_data/"
output_dir <- "replication_plots/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- Load Prerequisite Metadata ---
# This metadata is required for filtering, ordering, and coloring plots.

tissue_list <- c("Blood", "Lung", "Skin", "Colon", "Liver")
ancestry_list <- c("EUR", "EAS", "AFR", "AMR")

# Load cell type color and category information
cell_type_cat_colors <- read_delim("metadata/cell_type_colors_new.txt", delim = '\t')

# Create dictionaries for cell type names and colors
cell_type_short_dict <- lapply(tissue_list, function(x) {
  tmp <- cell_type_cat_colors$cell_type_short[cell_type_cat_colors$tissue == x]
  names(tmp) <- cell_type_cat_colors$cell_type[cell_type_cat_colors$tissue == x]
  tmp
})
names(cell_type_short_dict) <- tissue_list
cell_type_color_dict <- cell_type_cat_colors$cell_type_colors
names(cell_type_color_dict) <- cell_type_cat_colors$cell_type_short

# Define cell type lists per tissue (used for ordering factors)
celltypes_all_list <- list(
  Blood = c('CD4TNC','CD4TEM','Treg','CD8TNC','CD8TEM','CD8TEMRA','MAIT','NKp','NKn','BIN','BMem','Plasma','Plasmablasts','MonoC','MonoNC','Nph','DC','pDC'),
  Liver = c('CD4T',"CD8T","Treg","MAIT","Th","gdT",'Circulating_NK','Resident_NK','B','Plasma','Monocytes','Macrophages','DC','pDC','Neutrophils','Basophils','Hepatocytes','Cholangiocytes','Endothelium','Fibroblasts'),
  Skin = c("Th","Tc", "Treg","NK","DC1","DC2","MigDC","Macro1","Macro2","MonoMac","Mast","KCdiff","KCundiff","Melanocyte","VE1","VE2","VE3","LE1","Pericyte1","Pericyte2","F2"),
  Lung = c("CD4T","CD8T","NK","B","Monocytes","Mast","Macrophages","DC","AT1","AT2","Airway_Epi_Multiciliated","Airway_Epi_Basal","Airway_Epi_Secretory","EC_arterial","EC_capillary","EC_venous","LEC_mature","LEC_diff","Fibroblasts","Smooth_muscle","Mesothelium","Rare"),
  Colon = c("CD4T","CD8T","Treg","Th","gdT","NKT","NK","ILC3","BIN","BMem","Bcyc","Plasma","Mono","Macro","cDC2","Mast","Colonocyte","GLoblet","Tuft","TA", "EEC" ,"ECcap","ECven","Stromal1","Stromal2","Myofibroblast","Glia")
)

# Load and process summary statistics to define high-quality conditions
# (This is used later to filter the concordance data)
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


## 2. MAIN ANALYSIS: CALCULATE CONCORDANCE AND GENERATE PLOTS
# ------------------------------------------------------------------------------
print("--- Section 2: Calculating Concordance and Generating Plots ---")

sign_df_list <- list()
sign_df1_list <- list()

# Loop over each tissue and ancestry to process matched data
for (tissue in tissue_list) {
  for (anc in ancestry_list) {
    
    matched_file <- file.path(data_dir, paste0("matched_sceqtl_bulk_", tissue, "_", anc, ".rds"))
    
    if (file.exists(matched_file)) {
      print(paste0("Processing: ", tissue, " - ", anc))
      
      # Load and format data (handling slight differences in structure for Blood)
      if (tissue != "Blood") {
        df_comb_all <- readRDS(matched_file)[[tissue]]
        df_comb_all[, ancestry := anc]
        df_comb_all <- df_comb_all[maf_df >= 0.01 & maf_df <= 0.99]
      } else {
        df_celltype_list <- readRDS(matched_file)
        df_comb_all <- rbindlist(df_celltype_list)
        df_comb_all <- df_comb_all[af >= 0.01 & af <= 0.99]
        df_comb_all$tissue <- tissue
        df_comb_all$ancestry <- anc
      }
      
      # Filter for significant sc-eQTLs and join with metadata
      df_comb_all <- df_comb_all %>%
        dplyr::filter(pval_nominal_df < 5e-8) %>%
        dplyr::mutate(sig = ifelse(pval_nominal_bulk < 5e-8, "Significant", "Non-significant")) %>%
        left_join(cell_type_cat_colors %>% rename(celltype = cell_type) %>% select(tissue, celltype, cell_type_short), by = c("tissue", "celltype")) %>%
        inner_join(smr_list_all %>% filter(ancestry == anc) %>% select(tissue, ancestry, cell_type, cell_type_short), by = c("tissue", "ancestry", "cell_type_short"))
      
      df_comb_all$cell_type_short <- factor(df_comb_all$cell_type_short, levels = cell_type_short_dict[[tissue]][celltypes_all_list[[tissue]]])
      
      # --- Generate Scatter Plot ---
      color_pal <- c(cell_type_color_dict[levels(df_comb_all$cell_type_short)], grey = "#bebebe83")
      
      p_scatter <- df_comb_all %>%
        mutate(cell_type_color = ifelse(sig == "Significant", as.character(cell_type_short), "grey")) %>%
        ggplot(aes(x = slope_std, y = beta_std, color = cell_type_color)) +
        geom_point(alpha = 0.7) +
        geom_hline(yintercept = 0) +
        geom_vline(xintercept = 0) +
        scale_color_manual(values = color_pal) +
        facet_wrap(~cell_type_short, ncol = 4) +
        labs(x = "Standardized effect size of sc-eQTLs in our study", y = "Standardized effect size in bulk study") +
        theme_classic() +
        theme(text = element_text(size = 24), panel.grid = element_blank(), strip.background = element_rect(colour = "white"), legend.position = "none")
      
      ggsave(
        file.path(output_dir, paste0("scatter_sceqtl_bulk_", tissue, "_", anc, ".png")),
        plot = p_scatter,
        width = 16,
        height = ceiling(length(unique(df_comb_all$cell_type_short)) / 4) * 4
      )
      
      # --- Calculate Sign Concordance ---
      # Set 1: All significant sc-eQTLs
      sign_df <- df_comb_all %>%
        rowwise() %>%
        mutate(sign_concordance = (sign(slope_std) == sign(beta_std))) %>%
        ungroup() %>%
        group_by(cell_type_short) %>%
        summarise(Concordance = mean(sign_concordance, na.rm = TRUE), Disconcordance = 1 - mean(sign_concordance, na.rm = TRUE), n = n())
      sign_df$tissue <- tissue
      sign_df$ancestry <- anc
      sign_df_list[[length(sign_df_list) + 1]] <- sign_df
      
      # Set 2: eQTLs significant in both sc-eQTL and bulk
      sign_df1 <- df_comb_all %>%
        filter(sig == "Significant") %>%
        rowwise() %>%
        mutate(sign_concordance = (sign(slope_std) == sign(beta_std))) %>%
        ungroup() %>%
        group_by(cell_type_short) %>%
        summarise(Concordance = mean(sign_concordance, na.rm = TRUE), Disconcordance = 1 - mean(sign_concordance, na.rm = TRUE), n = n())
      sign_df1$tissue <- tissue
      sign_df1$ancestry <- anc
      sign_df1_list[[length(sign_df1_list) + 1]] <- sign_df1
    }
  }
}


## 3. AGGREGATE RESULTS AND GENERATE SUMMARY PLOT
# ------------------------------------------------------------------------------
print("--- Section 3: Aggregating Results and Generating Summary Plot ---")

# Aggregate and save concordance results
sign_df_all <- rbindlist(sign_df_list)
write.table(sign_df_all, file.path(output_dir, "summary_sign_concordance_sc_sig.txt"), sep = '\t', row.names = F, col.names = T, quote = F)

sign_df_all1 <- rbindlist(sign_df1_list)
write.table(sign_df_all1, file.path(output_dir, "summary_sign_concordance_sc_and_bulk_sig.txt"), sep = '\t', row.names = F, col.names = T, quote = F)

# Print summary statistics to console
print("Average concordance for sc-significant eQTLs:")
print(sign_df_all %>% group_by(tissue) %>% summarise(avg_conc = mean(Concordance)))
print("Average concordance for shared significant eQTLs:")
print(sign_df_all1 %>% group_by(tissue) %>% summarise(avg_conc = mean(Concordance)))

# --- Create Diverging Bar Plot ---
y_labels <- c("Blood-EUR", "Blood-EAS", "Blood-AFR", "Blood-AMR", "Lung-EUR", "Lung-EAS", "Skin-EUR", "Colon-EUR", "Liver-EUR", "Liver-EAS")

# Plot for sc-significant eQTLs (left side)
p1 <- sign_df_all %>%
  mutate(ts_anc = factor(paste0(tissue, "-", anc), levels = y_labels)) %>%
  group_by(ts_anc) %>%
  summarise(med_conc = median(Concordance), med_disconc = 1 - med_conc, .groups = 'drop') %>%
  pivot_longer(!ts_anc, values_to = "prop", names_to = "conc") %>%
  ggplot(aes(y = fct_rev(ts_anc), x = prop * 100, fill = factor(conc, levels = c("med_disconc", "med_conc")))) +
  geom_bar(stat = 'identity', width = 0.8, position = "stack") +
  theme_classic() +
  theme(
    text = element_text(size = 16, family = "Helvetica"),
    axis.text = element_text(size = 16, color = "black"),
    axis.title = element_text(size = 16, color = "black"),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    plot.margin = margin(l = 0.5, r = 0.5, unit = "cm"),
    legend.position = "none",
    axis.title.y = element_blank()
  ) +
  scale_fill_manual(values = c("lightgrey", "#f1bd69")) +
  scale_x_reverse(limits = c(100, 0), breaks = c(100, 75, 50, 25, 0), expand = expansion(add = c(0, 1.5)), position = "top") +
  labs(x = "sc-significant eQTLs", fill = "")

# Plot for shared significant eQTLs (right side)
p2 <- sign_df_all1 %>%
  mutate(ts_anc = factor(paste0(tissue, "-", anc), levels = y_labels)) %>%
  group_by(ts_anc) %>%
  summarise(med_conc = median(Concordance), med_disconc = 1 - med_conc, .groups = 'drop') %>%
  pivot_longer(!ts_anc, values_to = "prop", names_to = "conc") %>%
  ggplot(aes(y = fct_rev(ts_anc), x = prop * 100, fill = factor(conc, levels = c("med_disconc", "med_conc")))) +
  geom_bar(stat = 'identity', width = 0.8, position = "stack") +
  theme_classic() +
  theme(
    text = element_text(size = 16, family = "Helvetica"),
    axis.text = element_text(size = 16, color = "black"),
    axis.title = element_text(size = 16, color = "black"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(r = 0.5, l = 0.5, unit = "cm"),
    legend.position = "none",
    axis.title.y = element_blank()
  ) +
  scale_fill_manual(values = c("lightgrey", "#f1bd69")) +
  scale_x_continuous(limits = c(0, 100), position = "top", expand = expansion(add = c(0, 2))) +
  labs(x = "sc-bulk shared eQTLs", fill = "")

# Create a central axis with labels
label_df <- data.frame(ts_anc = factor(y_labels, levels = y_labels), y = seq_along(y_labels))
p_y <- ggplot(label_df, aes(y = ts_anc, x = 1)) +
  geom_blank() +
  scale_y_discrete(position = 'right') +
  theme_void() +
  theme(
    axis.text.y = element_text(size = 14, family = "Helvetica", hjust = 0.5, color = "black"),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = margin(l = 0.5, r = 0.5, 0, 0, "cm")
  )

# Combine the three plots into a single figure
final_plot <- plot_grid(p1, p_y, p2, ncol = 3, align = "h", axis = "tb", rel_widths = c(1, 0.3, 1))

# Save the final plot
ggsave(file.path(output_dir, "summary_sign_concordance_diverging_bar.png"), plot = final_plot, width = 8, height = 6)
