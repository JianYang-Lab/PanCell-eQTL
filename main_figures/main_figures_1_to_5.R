################################################################################
#
# Script: Generation of Main Figures 1-5
#
# Description:
# This script contains the R code to generate the main figures for the study.
# It covers a wide range of analyses and visualizations, including:
# - Figure 1: Genotype quality control, PCA, and cell type UMAPs.
# - Figure 2: Summary of eQTL discovery rates across tissues and ancestries.
# - Figure 3: Functional enrichment analysis of eVariants (chromHMM, SnpEff).
# - Figure 4: Replication and concordance with external datasets (OneK1K, GTEx).
# - Figure 5: Visualization of mashR results for eQTL sharing.
#
# This is a complete, runnable script organized from the original Rmd file.
#
################################################################################


## 1. SETUP & GLOBAL PARAMETERS
# ------------------------------------------------------------------------------
print("--- Section 1: Setup & Global Parameters ---")

# --- 1.1: Load Libraries ---
library(tidyverse)
library(data.table)
library(cowplot)
library(RColorBrewer)
library(ggtext)
library(ComplexHeatmap)
library(ggnewscale)
library(latex2exp)
library(arrow)
library(circlize)
library(grid)
library(ggvenn)
library(scatterpie)
library(parallel)
library(harmony)

# --- 1.2: Global Theme and Working Directory ---
theme_set(
  theme_classic() +
    theme(text = element_text(family = "Helvetica"), axis.text = element_text(color = "black"))
)

# Set the base directory for the project.
base_dir <- "/path/to/your/scRNA_project_folder/"
setwd(base_dir)

# Define output directory for figures and intermediate data
output_dir <- "final_figures_and_data/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- 1.3: Global Variables & Metadata ---
tissue_list <- c("Blood", "Lung", "Skin", "Colon", "Liver")
ancestry_list <- c("EUR", "EAS", "AFR", "AMR")
ancestry_color_list <- c(EUR = "#66C2A5", EAS = "#FC8D62", AFR = "#8DA0CB", AMR = "#E78AC3")

ancestries_all_list <- list(
  Blood = c("EUR", "EAS", "AFR", "AMR"), Lung = c("EUR", "EAS"), Skin = c("EUR"),
  Colon = c("EUR"), Liver = c("EUR", "EAS")
)

celltypes_all_list <- list(
  Blood = c('CD4TNC','CD4TEM','Treg','CD8TNC','CD8TEM','CD8TEMRA','MAIT','NKp','NKn','BIN','BMem','Plasma','Plasmablasts','MonoC','MonoNC','Nph','DC','pDC'),
  Liver = c('CD4T',"CD8T","Treg","MAIT","Th","gdT",'Circulating_NK','Resident_NK','B','Plasma','Monocytes','Macrophages','DC','pDC','Neutrophils','Basophils','Hepatocytes','Cholangiocytes','Endothelium','Fibroblasts'),
  Skin = c("Th","Tc", "Treg","NK","DC1","DC2","MigDC","Macro1","Macro2","MonoMac","Mast","KCdiff","KCundiff","Melanocyte","VE1","VE2","VE3","LE1","Pericyte1","Pericyte2","F2"),
  Lung = c("CD4T","CD8T","NK","B","Monocytes","Mast","Macrophages","DC","AT1","AT2","Airway_Epi_Multiciliated","Airway_Epi_Basal","Airway_Epi_Secretory","EC_arterial","EC_capillary","EC_venous","LEC_mature","LEC_diff","Fibroblasts","Smooth_muscle","Mesothelium","Rare"),
  Colon = c("CD4T","CD8T","Treg","Th","gdT","NKT","NK","ILC3","BIN","BMem","Bcyc","Plasma","Mono","Macro","cDC2","Mast","Colonocyte","GLoblet","Tuft","TA", "EEC" ,"ECcap","ECven","Stromal1","Stromal2","Myofibroblast","Glia")
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
  tmp[celltypes_all_list[[x]]]
})
names(cell_type_short_dict) <- tissue_list
cell_type_color_dict <- cell_type_cat_colors$cell_type_colors
names(cell_type_color_dict) <- cell_type_cat_colors$cell_type_short
lineage_colors <- c(`Immune Lymphoid` = "#8AB9D8", `Immune Myeloid` = "#e5a03d", Epithelial = "#8d4c28", Endothelial = "#b84e3c", Mesenchymal = "#c76b85", Other = "#c3e0e5")

# Load and process summary statistics to define high-quality conditions
smr_list <- list()
for (t in tissue_list) {
  smr_file <- file.path(t, "sc-eQTL/results/sig/eQTL_summary_5en8_maf01.tsv")
  if(file.exists(smr_file)){
    smr <- read_table(smr_file) %>% distinct()
    smr$tissue <- t
    smr_list[[t]] <- smr
  }
}

ts_ct_anc <- read_delim("metadata/tissue_celltype_ancestry_ss_new.txt", delim = '\t') %>% distinct()

smr_list_all <- list_rbind(smr_list) %>%
  filter(!(tissue == "Blood" & cell_type %in% c("T", "B", "Eryth"))) %>%
  left_join(ts_ct_anc, by = c("tissue", "cell_type", "ancestry")) %>%
  filter(sample_size >= 50) %>%
  inner_join(cell_type_cat_colors, by = c("tissue", "cell_type")) %>%
  mutate(
    ancestry = factor(ancestry, levels = ancestry_list),
    ct_anc = paste0(cell_type_short, "-", ancestry),
    total_count = round(sample_size * mean_count)
  ) %>%
  distinct() %>%
  filter(total_count > 10000) %>%
  rename(neGenes = neGenes_5en8, neQTLs = neQTLs_5en8) %>%
  filter(!(tissue == "Skin" & ancestry == "EAS"))


## 2. HELPER FUNCTIONS
# ------------------------------------------------------------------------------
print("--- Section 2: Defining Helper Functions ---")

# For Figure 1
get_pc_plot_batch <- function(ts, ref_pop_table) {
  pop_table <- read_delim(file.path(ts, "demographics", paste0(ts, "_sample_anc_pred_topmed.txt")), delim = '\t') %>% rename(ancestry = ancestry_pred)
  pop <- rbind(pop_table[, c("sample", "ancestry")], ref_pop_table)
  
  eig <- sqrt(read.table(file.path(ts, "VCF_TopMedimputed/filtered/rsq09/merged_ref_isec_pca300.eigenval"))[, 1])
  vec <- read.table(file.path(ts, "VCF_TopMedimputed/filtered/rsq09/merged_ref_isec_pca300.eigenvec"))
  
  nPCS <- ncol(vec) - 2
  for (k in 1:nPCS) {
    vec[, 2 + k] <- vec[, 2 + k] * eig[k]
  }
  vec <- data.frame(vec[, -1])
  colnames(vec) <- c("sample", paste0("PC", 1:300))
  
  vec <- vec %>%
    filter(sample %in% pop$sample) %>%
    left_join(pop, by = "sample") %>%
    filter(ancestry %in% c(ancestry_list, "Other"))
  vec$REF <- ifelse(vec$sample %in% ref_pop_table$sample, 1, 0)
  
  if (ts == "Blood") {
    vec$Category <- "scRNA"
    vec$Category[vec$REF == 1] <- "WGS"
    wgs_patterns <- "1k1k|popCell|PRJNA728702|LLDeep|TB.*_TB.*|HMN|EGAD00001008197|PRJNA671316"
    vec$Category[str_detect(vec$sample, wgs_patterns)] <- "WGS"
    PCs_corrected <- HarmonyMatrix(vec[, 2:301], meta_data = vec, vars_use = "Category")
  } else {
    PCs_corrected <- HarmonyMatrix(vec[, 2:301], meta_data = vec, vars_use = "REF")
  }
  
  PCs_corrected_df <- as.data.frame(PCs_corrected)
  PCs_corrected_df$ancestry <- vec$ancestry
  PCs_corrected_df$REF <- vec$REF
  
  if (ts == "Liver") {
    PCs_corrected_df$PC2 <- -PCs_corrected_df$PC2
  }
  
  p <- ggplot(PCs_corrected_df, aes(x = PC1, y = PC2, color = factor(ancestry, levels = c(ancestry_list, "Other")), shape = factor(REF, labels = c("This study", "1KGP")), size = factor(REF, labels = c("This study", "1KGP")))) +
    geom_point(alpha = 0.8) +
    scale_color_manual(values = c(ancestry_color_list, Other = "grey")) +
    scale_shape_manual(values = c(17, 1)) +
    scale_size_manual(values = c(2, 1)) +
    theme_nothing() + labs(x = "PC1", y = "PC2") +
    theme(legend.text = element_text(size = 12), legend.title = element_blank(), plot.title = element_text(hjust = 0.5)) +
    guides(shape = guide_legend(override.aes = list(size = 3)), color = guide_legend(override.aes = list(size = 3)))
  return(p)
}

# For Figure 2
eqtl_rank_count_plot <- function(ind_res_all, t, cell_type_short_dict, celltypes_all_list, cell_type_color_dict) {
  df <- ind_res_all %>% filter(tissue == t)
  df$cell_type_short <- factor(df$cell_type_short, levels = cell_type_short_dict[[t]])
  
  df <- df %>%
    mutate(
      ancestry = factor(ancestry, levels = ancestry_list),
      rank = factor(ifelse(as.numeric(Var1) >= 5, ">=5", Var1), levels = c(">=5", 4:1))
    ) %>%
    group_by(rank, ancestry, cell_type_short, ct_anc) %>%
    summarise(count = sum(Freq), .groups = 'drop')
  
  alpha_level <- c(">=5" = 1, "4" = 0.8, "3" = 0.6, "2" = 0.4, "1" = 0.2)
  
  p <- df %>%
    ggplot(aes(fill = cell_type_short)) +
    geom_bar(aes(y = ct_anc, x = count / 1000, alpha = rank), stat = "identity", position = "stack") +
    facet_grid(vars(cell_type_short), scales = "free", space = "free") +
    theme_classic() +
    scale_alpha_manual(values = alpha_level) +
    scale_fill_manual(values = cell_type_color_dict) +
    theme(
      text = element_text(size = 12), axis.title = element_text(size = 12),
      axis.text = element_text(size = 12, color = "black"),
      strip.background.y = element_blank(), strip.text.y = element_blank(),
      plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
      axis.title.y = element_blank(), legend.position = c(0.75, 0.75)
    ) +
    labs(x = "Number of eGenes (k)", alpha = "# ind. eQTLs", title = t) +
    guides(fill = "none")
  return(p)
}

# For Figure 5
get_significant_top_results_new <- function(lfsr, all_pair, thresh = 0.05, condition) {
  lfsr_sub <- lfsr[, condition, with = FALSE]
  lfsr_sub <- cbind(lfsr_sub, all_pair)
  setnames(lfsr_sub, colnames(lfsr)[condition], "lfsr")
  top <- lfsr_sub[lfsr < thresh, .(row_num = .I[which.min(lfsr)]), by = c("phenotype_id")]$row_num
  return(top)
}

get_pairwise_sharing_scatter_plot <- function(pm, lfsr, pair, pairs, output_dir, factor = 0.5, lfsr_thresh = 0.05, FUN = identity) {
  col1 <- pairs[1]
  col2 <- pairs[2]
  i <- which(colnames(pm) == col1)
  j <- which(colnames(pm) == col2)
  
  sig_list <- list()
  sig_list[[i]] <- get_significant_top_results_new(lfsr, pair, thresh = lfsr_thresh, condition = i)
  sig_list[[j]] <- get_significant_top_results_new(lfsr, pair, thresh = lfsr_thresh, condition = j)
  
  a <- union(sig_list[[i]], sig_list[[j]])
  pm_tmp <- pm[a, c(i, j), with = FALSE]
  pm_tmp$share <- ifelse(pm_tmp[, get(col1)] / pm_tmp[, get(col2)] > 0.5 & pm_tmp[, get(col1)] / pm_tmp[, get(col2)] < 2, "shared", "not shared")
  prop <- mean(pm_tmp$share == "shared")

  p <- ggplot(pm_tmp, aes(x = get(col1), y = get(col2), color = share)) +
    geom_point() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
    geom_abline(slope = 2, intercept = 0, linetype = "dashed") +
    geom_abline(slope = 0.5, intercept = 0, linetype = "dashed") +
    scale_color_manual(values = c("grey", "royalblue4")) +
    annotate("text", size = 6, family = "Helvetica", x = -0.6, y = 0.8, label = paste0("Prop = ", format(prop * 100, digits = 4), "%")) +
    theme_classic() +
    theme(axis.title = element_text(size = 18), axis.text = element_text(size = 18, color = "black"), legend.position = "none", plot.title = element_text(size = 18, hjust = 0.5)) +
    scale_y_continuous(limits = c(-1, 1)) + scale_x_continuous(limits = c(-1, 1)) +
    labs(title = paste0("Blood-", str_replace(col1, "_", "-")), y = paste0("Blood-", str_replace(col2, "_", "-")), x = "")
  
  ggsave(file.path(output_dir, paste0("pairwise_scatter_", col1, "_", col2, ".png")), plot = p, width = 4.5, height = 4.5)
  return(p)
}


## 3. FIGURE 1: GENOTYPE QC & POPULATION STRUCTURE
# ------------------------------------------------------------------------------
print("--- Section 3: Generating Figure 1 ---")

# --- 3.1: Imputation Accuracy Plots ---
plot_grid_list1 <- list()
for (ts in tissue_list) {
  stats <- read_delim(file.path(ts, "VCF_TopMedimputed/stats/vcf_maf_r2_stats.txt"), delim = "\t")
  stats$bin <- factor(stats$bin, levels = stats$bin)
  p1 <- ggplot(stats, aes(x = bin, y = n / 1000000)) +
    geom_bar(width = 0.5, stat = "identity", fill = "darkorange") +
    theme_classic() + scale_y_continuous(limits = c(0, 6)) +
    theme(axis.title = element_text(size = 14), axis.ticks.x = element_blank(), axis.text.x = element_blank(),
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5)) +
    labs(x = "", y = "Variant count (million)", title = ts)
  
  p2 <- stats %>%
    pivot_longer(cols = c(mean, median), values_to = "value", names_to = "metric") %>%
    mutate(metric = factor(metric, levels = c("median", "mean"), labels = c("Median", "Mean"))) %>%
    ggplot(aes(x = bin, y = value, color = metric, group = metric, shape = metric)) +
    geom_point(size = 3) + scale_shape_manual(values = c(17, 16)) +
    scale_y_continuous(limits = c(0, 1)) + theme_classic() + scale_color_manual(values = c("orange", "royalblue2")) +
    theme(axis.title = element_text(size = 14), axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank()) +
    labs(x = "MAF bin", y = "Imputation accuracy")

  if (ts != "Blood") {
    p1 <- p1 + labs(y = "") + theme(axis.text.y = element_blank())
    p2 <- p2 + labs(y = "") + theme(axis.text.y = element_blank())
  }
  plot_grid_list1[[ts]] <- plot_grid(p1, p2, nrow = 2, rel_heights = c(0.8, 1), align = "v", axis = "lr")
}
p_fig1a <- plot_grid(plotlist = plot_grid_list1, nrow = 1, rel_widths = c(1.16, 1, 1, 1, 1.38))
ggsave(file.path(output_dir, "fig1_imputation_accuracy.png"), plot = p_fig1a, width = 17, height = 4.8)

# --- 3.2: Genotype PCA Plots ---
ref_pop_table <- read.table("metadata/1KGP_3202_sample_ancestry.txt", col.names = c("sample", "ancestry"))
plot_grid_list2 <- list()
for (ts in tissue_list) {
  p3 <- get_pc_plot_batch(ts, ref_pop_table)
  if (ts == "Liver") p3 <- p3 + theme(legend.position = "right")
  plot_grid_list2[[ts]] <- p3
}
p_fig1b <- plot_grid(plotlist = lapply(plot_grid_list2, function(x) x + theme(plot.margin = margin(r = 25))), nrow = 1, align = "h", axis = "b", rel_widths = c(1, 1, 1, 1, 1.42))
ggsave(file.path(output_dir, "fig1_genotype_pcs.png"), plot = p_fig1b, width = 16, height = 2.5)

# --- 3.3: UMAP Visualization ---
plot_grid_list3 <- list()
for (ts in tissue_list) {
  umap_df <- fread(file.path("metadata/umap_coords", paste0(ts, "_umap.txt.gz")))
  # Custom label adjustments from original script
  if (ts == "Liver") {
    umap_df$cell_type_short[umap_df$cell_type_short == "Neutrophil"] <- "Nph"
    umap_df$cell_type_short[umap_df$cell_type_short == "Basophil"] <- "Bph"
  }
  p4 <- umap_df %>%
    mutate(cell_type_short = factor(cell_type_short, levels = cell_type_short_dict[[ts]])) %>%
    ggplot(aes(x = umap1, y = umap2, color = cell_type_short)) +
    geom_point(size = 0.013) +
    scale_color_manual(values = cell_type_color_dict) +
    theme_nothing() + theme(legend.position = "none")
  plot_grid_list3[[ts]] <- p4
}
p_fig1c <- plot_grid(plotlist = lapply(plot_grid_list3, function(x) x + theme(plot.margin = margin(r = 30))), nrow = 1, align = "h")
ggsave(file.path(output_dir, "fig1_umaps.png"), plot = p_fig1c, width = 16, height = 2.5)


## 4. FIGURE 2: EQTL DISCOVERY SUMMARY
# ------------------------------------------------------------------------------
print("--- Section 4: Generating Figure 2 ---")

ind_res_all1 <- readRDS("metadata/ind_res_all1_5en8.rds")
ind_res_all1 <- ind_res_all1 %>%
  inner_join(smr_list_all %>% select(tissue, ancestry, cell_type_short, neGenes, neQTLs, sample_size, total_count), by = c("tissue", "ancestry", "cell_type_short"))

p_list <- lapply(tissue_list, function(t) eqtl_rank_count_plot(ind_res_all1, t, cell_type_short_dict, celltypes_all_list, cell_type_color_dict))

p_fig2 <- plot_grid(
  p_list[[1]],
  plot_grid(p_list[[2]], p_list[[3]], nrow = 2, axis = "l", align = "v", rel_heights = c(1, 0.5)),
  plot_grid(p_list[[5]], p_list[[4]], nrow = 2, axis = "l", align = "v", rel_heights = c(1, 0.6)),
  ncol = 3, align = "h", axis = "l", rel_widths = c(1, 1, 1)
)
ggsave(file.path(output_dir, "fig2_independent_eqtl_counts.png"), plot = p_fig2, height = 12, width = 15, limitsize = FALSE)


## 5. FIGURE 3: FUNCTIONAL ENRICHMENT ANALYSIS
# ------------------------------------------------------------------------------
print("--- Section 5: Generating Figure 3 ---")

# --- 5.1: chromHMM Enrichment ---
chrom_grid <- read.table("metadata/chromHMM_grid.txt")
pattern <- chrom_grid %>% mutate(pattern = paste0(V1, "_5en8_", V2, "_", V3)) %>% pull(pattern)
fs <- paste0("functional_enrichment_chromHMM_", pattern, ".txt")

func_enrich_chrom <- lapply(fs, function(x) {
  df <- read_table(file.path("functional_enrichment/5en8/", x))
  ts <- str_split(str_split(x, "_")[[1]][6], "[.]")[[1]][1]
  anc <- str_split(x, "_")[[1]][4]
  df$tissue <- ts
  df$ancestry <- anc
  return(df)
}) %>% list_rbind()

chrom_state <- c("TssA", "TssAFlnk", "TxFlnk", "EnhG", "Tx", "TssBiv", "BivFlnk", "ZNF", "Enh", "TxWk", "EnhBiv", "Het", "Quies", "ReprPCWk", "ReprPC")
anno_cat_chrom <- list(Promoter = c("TssA", "TssAFlnk", "TssBiv", "BivFlnk"), Enhancer = c("EnhG", "EnhBiv", "Enh"), Transcription = c("Tx", "TxFlnk", "TxWk"), Repressed = c("ReprPC", "ReprPCWk"), `ZNF & Repeates` = c("ZNF"), Heterochromatin = c("Het"), Quiescent = c("Quies"))
anno_cat_df_chrom <- data.frame(chromatin_state = unlist(anno_cat_chrom), group = rep(names(anno_cat_chrom), lengths(anno_cat_chrom)))

func_enrich_chrom <- func_enrich_chrom %>% left_join(anno_cat_df_chrom, by = "chromatin_state") %>% filter(observed >= 0.001)

p1_chrom <- func_enrich_chrom %>%
  mutate(ts_anc = factor(paste0(tissue, "-", ancestry), levels = c("Blood-EUR", "Blood-EAS", "Blood-AFR", "Blood-AMR", 'Lung-EUR', 'Lung-EAS', 'Skin-EUR', "Colon-EUR", 'Liver-EUR', "Liver-EAS"))) %>%
  ggplot(aes(y = fct_rev(factor(chromatin_state, levels = chrom_state)), x = observed, fill = group)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_brewer(palette = "Set2") +
  facet_wrap(~ts_anc, nrow = 1) + theme_classic() +
  labs(x = "Proportion of eVariants", y = "", fill = "")

p2_chrom <- func_enrich_chrom %>%
  mutate(ts_anc = factor(paste0(tissue, "-", ancestry), levels = c("Blood-EUR", "Blood-EAS", "Blood-AFR", "Blood-AMR", 'Lung-EUR', 'Lung-EAS', 'Skin-EUR', "Colon-EUR", 'Liver-EUR', "Liver-EAS"))) %>%
  ggplot(aes(y = fct_rev(factor(chromatin_state, levels = chrom_state)), x = fold_enrichment, color = group)) +
  geom_point(size = 4, shape = 18) + scale_color_brewer(palette = "Set2") +
  geom_errorbar(aes(xmin = fold_enrichment - 1.96 * enrichment_se, xmax = fold_enrichment + 1.96 * enrichment_se), width = 0.2, color = "black") +
  facet_wrap(~ts_anc, nrow = 1) + theme_classic() + geom_vline(xintercept = 1, linetype = "dashed") +
  labs(x = 'Fold enrichment', y = "", color = "Chromatin state")

p_fig3a <- plot_grid(p1_chrom + theme(legend.position = "none"), p2_chrom, nrow = 2, align = "v", rel_heights = c(1, 1.04))
ggsave(file.path(output_dir, "fig3_chromhmm_enrichment.png"), plot = p_fig3a, width = 22, height = 10)

# --- 5.2: SnpEff Enrichment ---
fs_snpeff <- list.files("functional_enrichment/5en8", pattern = "snpEff")
fs_snpeff <- fs_snpeff[str_detect(fs_snpeff, "txt") & !str_detect(fs_snpeff, "level3")]

func_enrich_snpeff <- lapply(fs_snpeff, function(x) {
    df <- read_table(file.path("functional_enrichment/5en8/", x))
    ts <- str_split(str_split(x, "_")[[1]][6], "[.]")[[1]][1]
    anc <- str_split(x, "_")[[1]][4]
    df$tissue <- ts
    df$ancestry <- anc
    return(df)
}) %>% list_rbind()

anno_cat_snpeff <- list(`UTR Variants` = c("5-prime UTR", "3-prime UTR"), `Coding Sequence Variants` = c("synonymous", "missense", "stop lost", "stop retained", "start lost", "disruptive inframe deletion", "stop gained", "frameshift", "conservative inframe insertion", "conservative inframe deletion"), `Splicing Variants` = c("splice region", "splice acceptor", "splice donor"), `Intergenic Intronic Variants` = c("intron", "intergenic region"), `Transcript Level Variants` = c("nc-transcript"), `Regulatory Region Variants` = c("downstream gene", "upstream gene"))
anno_cat_df_snpeff <- data.frame(chromatin_state = unlist(anno_cat_snpeff), group = rep(names(anno_cat_snpeff), lengths(anno_cat_snpeff)))

func_enrich_snpeff <- func_enrich_snpeff %>%
  mutate(chromatin_state = str_replace(str_replace_all(chromatin_state, "_", " "), " variant", ""))
func_enrich_snpeff$chromatin_state[func_enrich_snpeff$chromatin_state == "5 prime UTR"] <- "5-prime UTR"
func_enrich_snpeff$chromatin_state[func_enrich_snpeff$chromatin_state == "3 prime UTR"] <- "3-prime UTR"
func_enrich_snpeff$chromatin_state[func_enrich_snpeff$chromatin_state == "non coding transcript"] <- "nc-transcript"

func_enrich_snpeff <- func_enrich_snpeff %>%
  left_join(anno_cat_df_snpeff, by = "chromatin_state") %>%
  mutate(chromatin_state = factor(chromatin_state, levels = c("5-prime UTR", "splice region", "3-prime UTR", "nc-transcript", "synonymous", "missense", "intron", "upstream gene", "downstream gene", "intergenic region"))) %>%
  filter(!is.na(chromatin_state), fold_enrichment != 0, observed > 0.001)

p1_snpeff <- func_enrich_snpeff %>%
  mutate(ts_anc = factor(paste0(tissue, "-", ancestry), levels = c("Blood-EUR", "Blood-EAS", "Blood-AFR", "Blood-AMR", 'Lung-EUR', 'Lung-EAS', 'Skin-EUR', "Colon-EUR", 'Liver-EUR', "Liver-EAS"))) %>%
  ggplot(aes(y = fct_rev(chromatin_state), x = observed, fill = group)) +
  geom_bar(width = 0.6, stat = "identity") + scale_fill_brewer(palette = "Set2") +
  facet_wrap(~ts_anc, nrow = 1) + theme_classic() +
  labs(x = "Proportion of eVariants", y = "", fill = "")

p2_snpeff <- func_enrich_snpeff %>%
  mutate(ts_anc = factor(paste0(tissue, "-", ancestry), levels = c("Blood-EUR", "Blood-EAS", "Blood-AFR", "Blood-AMR", 'Lung-EUR', 'Lung-EAS', 'Skin-EUR', "Colon-EUR", 'Liver-EUR', "Liver-EAS"))) %>%
  ggplot(aes(y = fct_rev(chromatin_state), x = fold_enrichment, color = group)) +
  geom_point(size = 3, shape = 18) + scale_color_brewer(palette = "Set2") +
  geom_errorbar(aes(xmin = fold_enrichment - 1.96 * enrichment_se, xmax = fold_enrichment + 1.96 * enrichment_se), width = 0.2, color = "black") +
  facet_wrap(~ts_anc, nrow = 1) + theme_classic() + geom_vline(xintercept = 1, linetype = "dashed") + scale_x_continuous(limits = c(0.75, 2)) +
  labs(x = 'Fold enrichment', y = "", color = "snpEff variant annotation")

p_fig3b <- plot_grid(p1_snpeff + theme(legend.position = "none"), p2_snpeff, nrow = 2, align = "v", rel_heights = c(1, 1.05))
ggsave(file.path(output_dir, "fig3_snpeff_enrichment.png"), plot = p_fig3b, width = 22, height = 8)


## 6. FIGURE 4: REPLICATION AND CONCORDANCE
# ------------------------------------------------------------------------------
print("--- Section 6: Generating Figure 4 ---")

# --- 6.1: Indep1K vs. OneK1K Comparison ---
print("Step 6.1: Plotting Indep1K vs. OneK1K replication.")

# Define path to pre-computed replication data
replication_data_dir <- "replication_data/"

# Load matched data
df_comb_all_list_1k1k <- readRDS(file.path(replication_data_dir, "Blood_sub_1k1k_effect_size.rds"))
df_celltype_all_1k1k <- rbindlist(df_comb_all_list_1k1k) %>%
  filter(af_df >= 0.01 & af_df <= 0.99) %>%
  mutate(sig = ifelse(fdr < 0.05, "Significant", "Non-significant"))

# Prepare for plotting
ct_order_1k1k <- cell_type_short_dict$Blood[celltypes_all_list$Blood[celltypes_all_list$Blood %in% unique(df_celltype_all_1k1k$celltype)]]
df_celltype_all_1k1k$cell_type_short <- factor(cell_type_short_dict$Blood[df_celltype_all_1k1k$celltype], levels = ct_order_1k1k)
color_pal_1k1k <- c(cell_type_color_dict[levels(df_celltype_all_1k1k$cell_type_short)], grey = "#bebebe83")

# Plot 1: Scatter plot of effect sizes
p_scatter_1k1k <- df_celltype_all_1k1k %>%
  mutate(cell_type_color = ifelse(sig == "Significant", as.character(cell_type_short), "grey")) %>%
  ggplot(aes(x = slope_std, y = betae_std, color = cell_type_color)) +
  geom_point(size = 1) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  scale_color_manual(values = color_pal_1k1k) +
  facet_wrap(~cell_type_short, ncol = 4) +
  labs(x = "Standardized effect size of sc-eQTLs in Indep1K", y = "Standardized effect size in OneK1K", color = "") +
  theme_classic() +
  theme(text = element_text(size = 12, family = "Helvetica"), axis.text = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 12, color = 'black'), panel.grid = element_blank(),
        strip.background = element_rect(colour = "white"), legend.position = 'none', strip.text = element_text(size = 14))
ggsave(file.path(output_dir, "fig4_indep1k_vs_onek1k_scatter.png"), plot = p_scatter_1k1k, width = 7, height = 6)

# Plot 2: Sign concordance bar plot
sign_df_1k1k <- df_celltype_all_1k1k %>%
  rowwise() %>%
  mutate(sign_concordance = (sign(slope_std) == sign(beta))) %>%
  ungroup() %>%
  group_by(cell_type_short) %>%
  summarise(Concordance = mean(sign_concordance, na.rm = TRUE), Disconcordance = 1 - mean(sign_concordance, na.rm = TRUE), n = n())

p_concord_1k1k <- sign_df_1k1k %>%
  pivot_longer(cols = c(Concordance, Disconcordance), names_to = "Concordance_type", values_to = "value") %>%
  mutate(Concordance_type = factor(Concordance_type, levels = c("Disconcordance", "Concordance"))) %>%
  ggplot(aes(y = value * 100, x = cell_type_short, fill = interaction(cell_type_short, Concordance_type))) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = c(
    setNames(rep("#d1d1d1", length(unique(sign_df_1k1k$cell_type_short))), paste0(unique(sign_df_1k1k$cell_type_short), ".Disconcordance")),
    setNames(cell_type_color_dict[unique(sign_df_1k1k$cell_type_short)], paste0(unique(sign_df_1k1k$cell_type_short), ".Concordance"))
  )) +
  labs(fill = "", x = "Cell type", y = "Sign concordance (%)") +
  theme_classic() +
  theme(text = element_text(size = 12, family = "Helvetica"), legend.position = "none",
        axis.text = element_text(size = 12, color = "black"), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_text(size = 12, color = "black"))
ggsave(file.path(output_dir, "fig4_indep1k_vs_onek1k_concordance.png"), plot = p_concord_1k1k, width = 5, height = 3)

# Plot 3: eGene Venn diagram
egene_list_1k1k <- readRDS(file.path(replication_data_dir, "egene_Blood_sub_tk1k_1k1k.rds"))
p_venn_1k1k <- ggvenn(
  list(Indep1K = unique(unlist(egene_list_1k1k$egene_df)), OneK1K = unique(unlist(egene_list_1k1k$egene_1k1k))),
  fill_color = c("#b0ccff", "#2646af"), fill_alpha = 0.8, auto_scale = TRUE, text_size = 4, set_name_size = 0
)
ggsave(file.path(output_dir, "fig4_indep1k_vs_onek1k_egene_venn.png"), plot = p_venn_1k1k, width = 5, height = 3)

# Plot 4: eGene count comparison
celltype_df_list <- c("CD4TNC", "CD4TEM", "CD8TNC", "CD8TEM", "NKp", "BIN", "BMem", "Plasma", "MonoC", "MonoNC", "DC")
egene_counts_1k1k <- lapply(1:length(egene_list_1k1k$egene_tk1k), function(x) {
  data.frame(neGene_1k1k = length(egene_list_1k1k$egene_1k1k[[x]]), neGene_df = length(egene_list_1k1k$egene_df[[x]]))
}) %>% rbindlist()
egene_counts_1k1k$celltype <- celltype_df_list

p_counts_1k1k <- egene_counts_1k1k %>%
  pivot_longer(cols = !celltype, names_to = 'group', values_to = 'neGene') %>%
  mutate(
    group = factor(group, levels = c("neGene_df", "neGene_1k1k"), labels = c("Indep1K", "OneK1K")),
    celltype = factor(celltype, levels = celltypes_all_list[["Blood"]][celltypes_all_list[["Blood"]] %in% celltype_df_list])
  ) %>%
  ggplot(aes(x = celltype, y = neGene, group = group, fill = group)) +
  geom_bar(stat = 'identity', position = 'dodge', width = 0.6) +
  theme_classic() + scale_fill_manual(values = c("#b0ccff", "#2646af")) +
  labs(x = 'Cell type', y = 'No. eGenes (FDR<0.05)', fill = "") +
  scale_y_continuous(labels = scales::label_comma()) +
  theme(text = element_text(family = "Helvetica"), axis.title = element_text(size = 12),
        axis.text = element_text(size = 12, color = 'black'), axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90),
        legend.text = element_text(size = 12))
ggsave(file.path(output_dir, "fig4_indep1k_vs_onek1k_egene_counts.png"), plot = p_counts_1k1k, width = 6, height = 3.5)

# --- 6.2: Bulk Tissue Concordance Summary ---
print("Step 6.2: Plotting bulk tissue concordance summary.")

sign_df_all <- fread(file.path(replication_data_dir, "matched_sceqtl_bulk_sign.txt"))
sign_df_all1 <- fread(file.path(replication_data_dir, "matched_sceqtl_bulk_sign_allsig.txt"))
y_labels <- rev(c("Blood-EUR", "Blood-EAS", "Blood-AFR", "Blood-AMR", "Lung-EUR", "Lung-EAS", "Skin-EUR", "Colon-EUR", "Liver-EUR", "Liver-EAS"))

p1 <- sign_df_all %>%
  mutate(ts_anc = factor(paste0(tissue, "-", ancestry), levels = rev(y_labels))) %>%
  group_by(ts_anc) %>%
  reframe(med_conc = median(Concordance), med_disconc = 1 - med_conc) %>%
  pivot_longer(!ts_anc, values_to = "prop", names_to = "conc") %>%
  ggplot(aes(y = ts_anc, x = prop * 100, fill = factor(conc, levels = c("med_disconc", "med_conc")))) +
  geom_bar(stat = 'identity', width = 0.8, position = "stack") +
  theme_classic() +
  theme(text = element_text(size = 16), axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none", axis.title.y = element_blank()) +
  scale_fill_manual(values = c("lightgrey", "#f1bd69")) +
  scale_x_reverse(limits = c(100, 0), breaks = c(100, 75, 50, 25, 0), expand = expansion(add = c(1.5, 0)), position = "top") +
  labs(x = "sc-significant eQTLs", fill = "")

p2 <- sign_df_all1 %>%
  mutate(ts_anc = factor(paste0(tissue, "-", ancestry), levels = rev(y_labels))) %>%
  group_by(ts_anc) %>%
  reframe(med_conc = median(Concordance), med_disconc = 1 - med_conc) %>%
  pivot_longer(!ts_anc, values_to = "prop", names_to = "conc") %>%
  ggplot(aes(y = ts_anc, x = prop * 100, fill = factor(conc, levels = c("med_disconc", "med_conc")))) +
  geom_bar(stat = 'identity', width = 0.8, position = "stack") +
  theme_classic() +
  theme(text = element_text(size = 16), axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none", axis.title.y = element_blank()) +
  scale_fill_manual(values = c("lightgrey", "#f1bd69")) +
  scale_x_continuous(limits = c(0, 100), position = "top", expand = expansion(add = c(2, 0))) +
  labs(x = "sc-bulk shared eQTLs", fill = "")

label_df <- data.frame(ts_anc = factor(y_labels, levels = y_labels))
p_y <- ggplot(label_df, aes(y = ts_anc, x = 1)) +
  geom_blank() + scale_y_discrete(position = 'right') + theme_void() +
  theme(axis.text.y = element_text(size = 12, family = "Helvetica", hjust = 0.5), plot.margin = margin(l = -10, r = -10))

p_bulk_conc <- plot_grid(p1, p_y, p2, ncol = 3, align = "h", axis = "tb", rel_widths = c(1, 0.4, 1))
ggsave(file.path(output_dir, "fig4_bulk_concordance_summary.png"), plot = p_bulk_conc, width = 7.2, height = 4)

# --- 6.3: eQTLGen eGene Discovery Analysis ---
print("Step 6.3: Plotting eQTLGen eGene discovery.")

out_eqtlgen <- readRDS(file.path(replication_data_dir, "Blood_eQTLGen_eGene.rds"))

p1_eqtlgen <- out_eqtlgen %>%
  group_by(egene_df_all_sig_in_bulk, expr_decile) %>% reframe(count = n()) %>%
  mutate(group = ifelse(egene_df_all_sig_in_bulk, "eGene in eQTLGen", "Not an eGene in eQTLGen"), expr_decile = factor(expr_decile, labels = paste0(seq(10, 100, 10), "%"))) %>%
  ggplot(aes(x = expr_decile, y = count, fill = group)) +
  geom_bar(stat = 'identity', position = 'dodge', width = 0.8) + theme_classic() +
  scale_x_discrete(breaks = c("50%", "100%")) + scale_y_continuous(labels = scales::label_comma()) +
  scale_fill_manual(values = c("#f1bd69", "#2A629A")) +
  theme(text = element_text(size = 16), legend.position = "none", axis.text = element_text(color = "black", size = 14), axis.title = element_text(size = 14)) +
  labs(fill = "", x = "Max expression decile", y = "No. eGenes")

p2_eqtlgen <- out_eqtlgen %>%
  group_by(egene_df_all_sig_in_bulk, ct_count) %>% reframe(count = n()) %>%
  mutate(group = ifelse(egene_df_all_sig_in_bulk, "eGene in eQTLGen", "Not an eGene in eQTLGen")) %>%
  ggplot(aes(x = ct_count, y = count, fill = group)) +
  geom_bar(stat = 'identity', width = 0.8, position = 'dodge') +
  scale_fill_manual(values = c("#f1bd69", "#2A629A")) + theme_classic() +
  scale_y_continuous(labels = scales::label_comma()) +
  theme(text = element_text(size = 16), axis.text = element_text(color = "black", size = 12), axis.title = element_text(size = 14), legend.position = c(0.4, 0.8), legend.text = element_text(size = 14)) +
  labs(fill = "", x = "Number of cell types gene is tested in", y = "No. eGenes")

p_eqtlgen <- plot_grid(p1_eqtlgen, p2_eqtlgen, nrow = 2)
ggsave(file.path(output_dir, "fig4_eQTLGen_eGene_stats.png"), plot = p_eqtlgen, width = 4.8, height = 5.4)

# --- 6.4: Natri Lung eGene/eQTL Count Comparison ---
print("Step 6.4: Plotting Natri Lung eGene/eQTL counts.")

stats_natri <- readRDS(file.path(replication_data_dir, "matched_sceqtl_natri_Lung_EAS_egene_count_new.rds"))
stats_natri$cell_type_short <- factor(cell_type_short_dict[["Lung"]][stats_natri$celltype_df], levels = cell_type_short_dict[["Lung"]])

p1_natri <- stats_natri %>%
  ggplot(aes(x = negene_natri, y = negene_df, color = cell_type_short)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_color_manual(values = cell_type_color_dict) +
  labs(x = "No. eGenes in Natri et al. (EUR)", y = "No. eGenes in Lung-EAS", color = "Cell type") +
  geom_point(size = 2.5) + theme_classic() +
  theme(axis.text = element_text(size = 14, color = "black"), axis.title = element_text(size = 14, color = "black"), text = element_text(size = 14))

p2_natri <- stats_natri %>%
  ggplot(aes(x = neqtl_natri, y = neqtl_df, color = cell_type_short)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_color_manual(values = cell_type_color_dict) +
  labs(x = "No. eQTLs in Natri et al. (EUR)", y = "No. eQTLs in Lung-EAS", color = "Cell type") +
  scale_y_continuous(labels = scales::label_comma()) +
  geom_point(size = 2.5) + theme_classic() +
  theme(axis.text = element_text(size = 14, color = "black"), axis.title = element_text(size = 14, color = "black"), text = element_text(size = 14))

p_natri <- plot_grid(
  p1_natri + theme(legend.position = "none", plot.margin = margin(10, 10, 5, 10)),
  p2_natri + theme(legend.position = "none", plot.margin = margin(10, 10, 5, 10)),
  nrow = 2, axis = "lr", align = "v"
)
ggsave(file.path(output_dir, "fig4_natri_lung_counts.png"), plot = p_natri, width = 3.6, height = 5.4)


## 7. FIGURE 5: MASH RESULTS VISUALIZATION
# ------------------------------------------------------------------------------
print("--- Section 7: Generating Figure 5 ---")


# Define paths for mashR and heterogeneity results
mashr_intra_tissue_dir <- "mashR/intra_tissue/"
mashr_cross_tissue_dir <- "mashR/cross_tissue/"
heterogeneity_dir <- "heterogeneity_results/cross_tissue/"

# --- 7.1: Intra-Tissue Sharing Heatmap (Blood) ---
print("Step 7.1: Generating intra-tissue mashR heatmap for Blood.")
ht_opt$HEATMAP_LEGEND_PADDING <- unit(4, "cm")
ht_opt$ANNOTATION_LEGEND_PADDING <- unit(4, "cm")

for (tissue in tissue_list[1]){
    print(tissue)
    m.pairwise_PM <- readRDS(paste0("crossTissue_analysis/mashR/v9/step6_pairwise_PM_all_pairwisetop_",tissue,".rds"))
    #m.pairwise_PM <- readRDS(paste0("crossTissue_analysis/mashR/v9/step5_pairwise_PM_",tissue,".rds"))
    ct_anc_grid=expand.grid(celltypes_all_list[[tissue]],ancestries_all_list[[tissue]])
    colnames(ct_anc_grid) = c("cell_type","ancestry")
    ct_anc_grid = ct_anc_grid %>% 
        dplyr::left_join(smr_list_all[smr_list_all$tissue == tissue,] %>% dplyr::select(cell_type,ancestry,lineage,cell_type_short,cell_type_colors,sample_size, total_count, neGenes, neQTLs)%>%distinct(),by=c("ancestry","cell_type")) %>% 
        dplyr::mutate(ct_anc = paste0(cell_type,"_",ancestry), ct_anc1 = paste0(cell_type_short,"_",ancestry)) %>% 
        dplyr::filter(ct_anc %in% rownames(m.pairwise_PM) & !is.na(cell_type_short)) %>% 
        dplyr::filter(!(tissue=="Skin" & ancestry=="EAS"))

    m.pairwise_PM = m.pairwise_PM[rownames(m.pairwise_PM) %in% ct_anc_grid$ct_anc, rownames(m.pairwise_PM) %in% ct_anc_grid$ct_anc] 
    m.pairwise_PM = m.pairwise_PM[ct_anc_grid$ct_anc,ct_anc_grid$ct_anc]
    rownames(m.pairwise_PM) = colnames(m.pairwise_PM) = ct_anc_grid$ct_anc1
   

    row_annotations <- data.frame(
        celltype = factor(ct_anc_grid$cell_type_short, levels=cell_type_short_dict[[tissue]][celltypes_all_list[[tissue]]]),
        ancestry = factor(ct_anc_grid$ancestry, levels=ancestries_all_list[[tissue]]),
        sample_size = ct_anc_grid$sample_size,
        total_count = ct_anc_grid$total_count,
        lineage = ct_anc_grid$lineage,
        neQTLs = ct_anc_grid$neQTLs,
        neGenes = ct_anc_grid$neGenes
    ) %>% arrange(ancestry,celltype)

    # Set colors
    celltype_colors <- setNames(ct_anc_grid$cell_type_colors, nm=ct_anc_grid$ct_anc1)
    ancestry_colors <- ancestry_color_list[ancestries_all_list[[tissue]]] 

    col_fun_sample_size <- colorRamp2(
        c(min(row_annotations$sample_size), max(row_annotations$sample_size)),
        c("#fdfdfd", "#d02e7d")  # Light gray to dark purple
    )

    col_fun_total_count <- colorRamp2(
        c(min(log10(row_annotations$total_count)), max(log10(row_annotations$total_count))),
        c("#fdfdfd", "#7136c2")  # Light gray to dark purple
    )

    col_fun_eGenes <- colorRamp2(
        c(min(log10(row_annotations$neGenes)), max(log10(row_annotations$neGenes))),
        c("#fdfdfd", "#35b8bb")  # Light orange to red
    )

    col_fun_eQTLs <- colorRamp2(
        c(min(log10(row_annotations$neQTLs)), max(log10(row_annotations$neQTLs))),
        c("#fdfdfd", "#ecbf1d")  # Light cyan to dark green
    )
    # Create row annotation with blocks
    col_anno <- columnAnnotation(
        `No. Indiv.` = anno_simple(row_annotations$sample_size,col=col_fun_sample_size, simple_anno_size = unit(sqrt(nrow(m.pairwise_PM))/3.5, "cm")),
        `log10(nCells)` = anno_simple(log10(row_annotations$total_count),col=col_fun_total_count, simple_anno_size = unit(sqrt(nrow(m.pairwise_PM))/3.5, "cm")),
        `log10(neGenes)` = anno_simple(log10(row_annotations$neGenes),col=col_fun_eGenes, simple_anno_size = unit(sqrt(nrow(m.pairwise_PM))/3.5, "cm")), 
        `log10(neQTLs)` =  anno_simple(log10(row_annotations$neQTLs),col=col_fun_eQTLs, simple_anno_size = unit(sqrt(nrow(m.pairwise_PM))/3.5, "cm")),
        `Lineage` = anno_simple(row_annotations$lineage,col=lineage_colors, simple_anno_size = unit(sqrt(nrow(m.pairwise_PM))/3.5, "cm")),
        ancestry=anno_block(
            gp=gpar(fill=ancestry_colors,col="white"), 
            labels=ancestries_all_list[[tissue]], 
            labels_gp=gpar(fontsize=38, fontface="bold",col="white"), 
            height=unit(sqrt(nrow(m.pairwise_PM))/3.5, "cm")),
        #`ancestry` = row_annotations$ancestry,
        annotation_name_gp = gpar(fontsize=38, fontface="bold"),
        annotation_name_side = "left",
        annotation_label = list(
        `No. Indiv.` = expression(bold("No. Indiv.")),
        `log10(nCells)` = expression(bold(log[10]("nCells"))),
        `log10(neGenes)` = expression(bold(log[10]("neGenes"))),
        `log10(neQTLs)` = expression(bold(log[10]("neQTLs")))
    ),
        simple_anno_size = unit(sqrt(nrow(m.pairwise_PM))/3.5, "cm")
    ) 
    row_anno <- rowAnnotation(
        celltype = anno_text(
            row_annotations$celltype, 
            just = "right",  # Right-align the text
            location = unit(0.95, "npc"),  # Set the anchor point to the right edge
            gp = gpar(col=celltype_colors, fontsize=38, fontface="bold")),
        show_annotation_name = FALSE
    )

    right_anno <- rowAnnotation(
        ancestry=anno_block(gp=gpar(fill=ancestry_colors,col="white"), 
                            labels=ancestries_all_list[[tissue]], 
                            labels_gp=gpar(fontsize=38, fontface="bold",col="white"),
                            width=unit(sqrt(nrow(m.pairwise_PM))/3.5, "cm"))
    )

    bottom_anno <- columnAnnotation(
        celltype = anno_text(
            row_annotations$celltype, 
            gp = gpar(col=celltype_colors, fontsize=38, fontface="bold")),
        show_annotation_name = FALSE
    )
    col_fun = colorRamp2(c(quantile(as.vector(m.pairwise_PM ))*100), c("#4aa7ff", "#fff2d2","#ffd0c2","#ff8e6c","#ff3c00"))

    create_legend <- function(title,col_fun,title_fontsize=40,labels_fontsize=36,legend_width=unit(nrow(m.pairwise_PM)/3.6, "cm"),legend_height=unit(nrow(m.pairwise_PM)/25, "cm"),grid_width=unit(nrow(m.pairwise_PM)/4.2, "cm"),grid_height=unit(nrow(m.pairwise_PM)/25, "cm")){
        return(Legend(
            title = title, 
            col_fun = col_fun, 
            title_gp = gpar(fontsize = title_fontsize, fontface="bold"),
            labels_gp = gpar(fontsize = labels_fontsize, fontface="bold"),
            legend_width = legend_width,
            legend_height = legend_height,
            grid_width = grid_width,
            grid_height = grid_height,
            direction = "horizontal",
            title_position = "topleft"
        ) )
    }
    # Draw the heatmap

    png(file.path(output_dir, paste0("fig5_mashr_heatmap_", tissue, ".png")), width = sqrt(nrow(m.pairwise_PM))*450, height = sqrt(nrow(m.pairwise_PM))*450)

    h <- Heatmap(
        m.pairwise_PM*100, name = "Proportion",
        na_col = "white", 
        top_annotation = col_anno,
        bottom_annotation = bottom_anno,
        right_annotation = right_anno, 
        left_annotation = row_anno,
        show_row_names = FALSE,
        show_column_names = FALSE,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_split = row_annotations[,c("ancestry")],   # Split rows by ancestry
        column_split = row_annotations[,c("ancestry")], # Split columns by ancestry
        width = unit(nrow(m.pairwise_PM)*1.35, "cm"),
        height = unit(nrow(m.pairwise_PM)*1.35, "cm"),
        show_heatmap_legend = FALSE,
        row_title = NULL,
        column_title = NULL,
        column_title_side = "bottom", 
        row_gap = unit(6, "mm"),  # Adjust the gap between row blocks
        column_gap = unit(6, "mm"),
        col = col_fun,
        rect_gp = gpar(col = "white") # Set grid border color to white
    )
    legend_lineage = Legend(
        title = "Lineage", 
        labels = unique(row_annotations$lineage),
        legend_gp = gpar(fill=lineage_colors[unique(row_annotations$lineage)]), 
        title_gp = gpar(fontsize=40, fontface="bold"),  # Slightly smaller font size
        labels_gp = gpar(fontsize=36, fontface="bold"),  # Slightly smaller font size
        grid_height = unit(nrow(m.pairwise_PM)/20, "cm"),
        grid_width = unit(nrow(m.pairwise_PM)/20, "cm"),
        direction = "horizontal",
        title_position = "topleft"
    ) 

    legend_sample_size = create_legend(title=expression(bold("No. Indiv.")), col_fun = col_fun_sample_size)

    legend_cell_count = create_legend(title=expression(bold(log[10]("nCells"))), col_fun = col_fun_total_count)

    legend_eGenes = create_legend(title = expression(bold(log[10]("neGenes"))), col_fun = col_fun_eGenes)

    legend_eQTLs = create_legend(title = expression(bold(log[10]("neQTLs"))), col_fun = col_fun_eQTLs)

    legend_heatmap = create_legend(title = "Pairwise sharing (%)",col_fun = col_fun)

    draw(h, heatmap_legend_side = "right", annotation_legend_side = "right" , merge_legend = TRUE, legend_gap = unit(nrow(m.pairwise_PM)/10, "cm"), annotation_legend_list = list(legend_sample_size, legend_cell_count, legend_eGenes,legend_eQTLs, legend_lineage, legend_heatmap)) # optionally combine them into one list 

    dev.off()
}


# --- 7.2: Pairwise Sharing Scatter Plots ---
print("Step 7.2: Generating example pairwise sharing scatter plots.")
out <- readRDS(file.path(mashr_intra_tissue_dir, paste0("step6_fit_all_results_", tissue, ".rds")))
all_pair <- lapply(1:22, function(x) readRDS(file.path(mashr_intra_tissue_dir, paste0("step5_all_pairs_", tissue, "_chr", x, ".rds")))) %>% rbindlist()
pm_all <- out$pm
lfsr_all <- out$lfsr
rm(out); gc()

sig_rows <- readRDS(file.path(mashr_intra_tissue_dir, paste0("step6_top_eQTL_rows_", tissue, ".rds")))
pm <- pm_all[sig_rows,]
lfsr <- lfsr_all[sig_rows,]
pair <- all_pair[sig_rows,]

get_pairwise_sharing_scatter_plot(pm, lfsr, pair, pairs = c("CD8TEM_EUR", "NKp_EUR"), output_dir)
get_pairwise_sharing_scatter_plot(pm, lfsr, pair, pairs = c("MAIT_EUR", "MAIT_AFR"), output_dir)
get_pairwise_sharing_scatter_plot(pm, lfsr, pair, pairs = c("CD8TEMRA_AFR", "MonoC_AFR"), output_dir)
get_pairwise_sharing_scatter_plot(pm, lfsr, pair, pairs = c("MonoNC_EAS", "CD8TNC_AFR"), output_dir)

# --- 7.3: Pairwise Sharing Boxplots ---
print("Step 7.3: Generating pairwise sharing summary boxplots.")
m.pairwise_PM_up <- as.matrix(m.pairwise_PM)
m.pairwise_PM_up[row(m.pairwise_PM_up) >= col(m.pairwise_PM_up)] <- NA

rb <- as.data.frame(m.pairwise_PM_up)
rb$anc_ct1 <- rownames(rb)
rb <- rb %>%
  pivot_longer(cols = !anc_ct1, values_to = "prop", names_to = "anc_ct2") %>%
  filter(!is.na(prop)) %>%
  rowwise() %>%
  mutate(
    celltype1 = str_split(anc_ct1, "_")[[1]][1], ancestry1 = str_split(anc_ct1, "_")[[1]][2],
    celltype2 = str_split(anc_ct2, "_")[[1]][1], ancestry2 = str_split(anc_ct2, "_")[[1]][2]
  ) %>%
  ungroup() %>%
  left_join(distinct(data.frame(celltype1 = ct_anc_grid$cell_type_short, celltype_meta1 = ct_anc_grid$lineage)), by = "celltype1") %>%
  left_join(distinct(data.frame(celltype2 = ct_anc_grid$cell_type_short, celltype_meta2 = ct_anc_grid$lineage)), by = "celltype2") %>%
  mutate(Case = case_when(
    (ancestry1 == ancestry2 & celltype_meta1 == celltype_meta2) ~ "Same both",
    (ancestry1 == ancestry2 & celltype_meta1 != celltype_meta2) ~ "Same ancestry",
    (ancestry1 != ancestry2 & celltype_meta1 == celltype_meta2) ~ "Same lineage",
    (ancestry1 != ancestry2 & celltype_meta1 != celltype_meta2) ~ "Different both"
  ))

p_boxplot_sharing <- rb %>%
  mutate(Case = factor(Case, levels = c("Same both", "Same lineage", "Same ancestry", "Different both"))) %>%
  ggplot(aes(x = Case, y = prop * 100)) +
  geom_violin(fill = "#3558c3") +
  geom_boxplot(color = "black", fill = "white", width = 0.1, linewidth = 0.5) +
  theme_classic() +
  labs(x = "", y = "Pairwise sharing (%)") +
  scale_y_continuous(limits = c(0, 100)) +
  theme(axis.text = element_text(size = 18), axis.title = element_text(size = 18),
        plot.title = element_text(size = 24, hjust = 0.5, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave(file.path(output_dir, "fig5_mashr_boxplot_sharing.png"), plot = p_boxplot_sharing, width = 1 * length(unique(rb$Case)), height = 4.5)

# --- 7.4: Cross-Tissue Sharing & Heterogeneity ---
print("Step 7.4: Summarizing cross-tissue sharing and heterogeneity.")

ct_cross_tissue <- read_delim("metadata/cross_tissue_celltype.txt", delim = "\t")
celltype_cross_tissue_list <- unique(ct_cross_tissue$celltype)

median_df_list <- list()
for (ct in celltype_cross_tissue_list) {
  rds_file <- file.path(mashr_cross_tissue_dir, paste0("step6_pairwise_PM_all_pairwisetop_", ct, ".rds"))
  if(file.exists(rds_file)){
    m.pairwise_PM <- readRDS(rds_file)
    median_df_list[[ct]] <- data.frame(celltype = ct, prop = m.pairwise_PM[lower.tri(m.pairwise_PM)])
  }
}
median_df <- rbindlist(median_df_list)

p_cross_tissue_sharing <- median_df %>%
  mutate(celltype = factor(celltype, levels = c("CD4T", "CD8T", "NK", "B", "Mono", "EC", "Epi"))) %>%
  ggplot(aes(x = celltype, y = prop * 100, color = celltype)) +
  geom_boxplot(width = 0.6) + geom_point() +
  scale_color_manual(values = c(`Epi` = "#8d4c28", cell_type_color_dict[celltype_cross_tissue_list])) +
  theme_classic() +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 16), legend.position = "none") +
  scale_y_continuous(limits = c(0, 100)) +
  labs(x = "Cell type", y = "Median pairwise sharing (%)", title = "Sign & Magnitude")
ggsave(file.path(output_dir, "fig5_cross_tissue_sharing.png"), plot = p_cross_tissue_sharing, width = 6, height = 4)

# Heterogeneity Test Summary
i2_all <- list()
for (ct in celltype_cross_tissue_list) {
  rds_file <- file.path(heterogeneity_dir, paste0("I2_REML_", ct, ".rds"))
  if (file.exists(rds_file)) {
    i2_list <- readRDS(rds_file)
    i2_list[, celltype := ct]
    i2_list$QEp.adj <- p.adjust(i2_list$QEp, method = "BH")
    i2_all[[ct]] <- i2_list
  }
}
i2_all <- rbindlist(i2_all)
i2_all[, het_cat := cut(I2, breaks = c(-1, 25, 50, 75, 100), labels = c("[0,25]", "(25,50]", "(50,75]", "(75,100]"))]

p_heterogeneity <- i2_all %>%
  mutate(celltype = factor(celltype, levels = c("CD4T", "CD8T", "NK", "B", "Mono", "EC", "Epi"))) %>%
  group_by(celltype) %>%
  reframe(n = n(), het_cat = het_cat) %>%
  group_by(celltype, het_cat) %>%
  reframe(het_cat_prop = n() / unique(n)) %>%
  ggplot(aes(y = het_cat_prop * 100, x = celltype, fill = het_cat)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6) +
  scale_fill_brewer(palette = 1) +
  scale_y_continuous(labels = scales::comma_format()) +
  theme_classic() +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.text = element_text(size = 16), legend.title = element_text(size = 14)) +
  labs(x = "Cell type", y = "Proportion of eQTLs (%)", fill = expression(paste(I^2, " (%)")))
ggsave(file.path(output_dir, "fig5_heterogeneity_summary.png"), plot = p_heterogeneity, width = 6, height = 3.6)

