# Usage: Rscript cell_count_sensitivity_public.R [project_root] [output_dir]
# Generating figure for Liver–EUR cell-count sensitivity analysis.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) args[1] else "."
results <- file.path(root, "Liver", "sc-eQTL_R", "results")
summary_dir <- file.path(results, "sensitivity_analysis")
out_dir <- if (length(args) >= 2) args[2] else summary_dir
thresholds <- c(5, 10, 20)

cell_order <- c(
  "CD4TNC", "CD4TEM", "Treg", "CD8TEMRA", "CD8TRM", "MAIT", "NKp", "NKn",
  "BIN", "BMem", "Plasma", "MonoC", "MonoNC", "Macro", "KC", "Mast", "Nph",
  "DC1", "DC2", "HepPP", "HepIZ", "HepPC", "CGC", "Endo_portal", "LSEC_PP",
  "LSEC_CV", "Fibro", "HSC"
)
metadata <- fread(file.path(root, "crossTissue_analysis_R", "tissue_ancestry_celltype_ss_cc_colors.txt"))
cells <- metadata[tissue == "Liver" & ancestry == "EUR" & sample_size >= 50 & total_count >= 10000,
                  .(cell_type, cell_type_short)]
labels <- cells$cell_type_short[order(match(cells$cell_type, cell_order))]
cells[, cell_type_short := factor(cell_type_short, levels = rev(labels))]

# A. Donors retained at each minimum cell count.
donors <- fread(file.path(summary_dir, "cell_count_ndonor.txt"))
donors <- tidyr::pivot_longer(cells[donors, on = "cell_type", nomatch = NULL], c(ct5, ct10, ct20),
                             names_to = "ct", names_prefix = "ct", values_to = "ndonor")
p_donors <- ggplot(donors, aes(factor(ct, levels = thresholds), cell_type_short, fill = ndonor)) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_text(aes(label = ndonor), size = 3, family = "Helvetica") +
  scale_fill_gradientn(colors = c("#f6ffeb", "#aed0a9", "#7DA78C", "#35858E")) +
  labs(x = "Cell count threshold", y = "Cell type", fill = "No. donors", tag = "A") +
  theme_minimal(base_size = 12) +
  theme(text = element_text(family = "Helvetica"), axis.text = element_text(color = "black"),
        panel.grid = element_blank(), axis.ticks = element_blank())

# B. eGenes with nominal P < 5e-8.
genes <- rbindlist(lapply(thresholds, function(ct) {
  fread(file.path(summary_dir, sprintf("eQTL_summary_5en8_ct%d.txt", ct)),
        select = c("celltype", "ct", "neGene"))
}))
setnames(genes, "celltype", "cell_type")
genes <- cells[genes, on = "cell_type", nomatch = NULL]
p_genes <- ggplot(genes, aes(factor(ct), cell_type_short, fill = log10(neGene))) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradientn(colors = c("#e1edf5", "#91d5ff", "#2f99db", "#2C5EAD")) +
  labs(x = "Cell Count Threshold", y = NULL, fill = "log10 neGene", tag = "B") +
  theme_minimal(base_size = 12) +
  theme(text = element_text(family = "Helvetica"), axis.text = element_text(color = "black"),
        panel.grid = element_blank())

# C. Correlate slopes for variant–gene pairs shared by each pair of thresholds.
comparisons <- c("5 vs 10", "5 vs 20", "10 vs 20")
slope_correlations <- function(cell_type) {
  slopes <- lapply(thresholds, function(ct) {
    suffix <- if (ct == 5) "" else paste0("c", ct)
    path <- file.path(results, "sig", paste0(cell_type, suffix, "_EUR_pc5_sig.tsv.gz"))
    x <- fread(path, select = c("phenotype_id", "variant_id", "slope"))
    setnames(x, "slope", paste0("ct", ct))
    x
  })
  shared <- Reduce(function(x, y) merge(x, y, by = c("phenotype_id", "variant_id"), all = TRUE), slopes)
  r <- combn(paste0("ct", thresholds), 2, function(pair) {
    cor(shared[[pair[1]]], shared[[pair[2]]], method = "pearson", use = "complete.obs")
  })
  data.table(cell_type, comparison = factor(comparisons, levels = comparisons), r = r)
}
eligible <- cell_order[cell_order %in% genes[ct == 20 & neGene > 0, cell_type]]
correlations <- rbindlist(lapply(eligible, slope_correlations))
p_corr <- ggplot(correlations, aes(comparison, r)) +
  geom_violin(fill = "royalblue") +
  geom_boxplot(fill = "white", color = "black", width = 0.1) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Cell count threshold comparison", y = "Pearson's correlation of shared eQTLs", tag = "C") +
  theme_classic(base_size = 12) +
  theme(text = element_text(family = "Helvetica"), axis.text = element_text(color = "black"))

figure <- (p_donors | p_genes | p_corr) + plot_layout(widths = c(1, 1, 2)) &
  theme(plot.tag.position = c(0, 1), plot.tag = element_text(family = "Helvetica", size = 18, face = "bold"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(out_dir, "figure_R.png"), figure, width = 9.8, height = 4.8, units = "in", dpi = 300)
