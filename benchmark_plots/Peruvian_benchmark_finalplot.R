#!/usr/bin/env Rscript
# Render Figure S23 (A–J). Requires ggplot2, ggvenn, cowplot, patchwork, hexbin, data.table and PLINK 1.9.
# Usage: Rscript Peruvian_benchmark_finalplot.R [project_root] [output_dir] [eqtl_table.rds]

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggvenn)
  library(cowplot)
  library(patchwork)
})
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
source(file.path(dirname(normalizePath(script)), "Peruvian_benchmark_prepare.R"))
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[1] else "."
benchmark <- file.path(root, "Peruvian_benchmark")
out_dir <- if (length(args) >= 2L) args[2] else file.path(benchmark, "plots")
eqtl_file <- if (length(args) >= 3L) args[3] else {
  file.path(root, "PeruvianscRNA", "sc-eQTL_R", "figures", "sig_list_all_new_new.rds")
}
eqtls <- if (length(args) >= 3L || file.exists(eqtl_file)) as.data.table(readRDS(eqtl_file)) else prepare_eqtls(root)
eqtls <- eqtls[ancestry == "AMR" & tissue %in% platforms & cell_type %in% names(cell_labels)]
cell_factor <- function(x) factor(x, levels = names(cell_labels), labels = unname(cell_labels))
cell_colors <- c(CD4TNC = "#dce6f2", CD4TEM = "#b9cde6", CD8TRM = "#3b6ebb",
                 CD8TEMRA = "#2c61b4", Treg = "#94b1da")
platform_colors <- c("SNP-array" = "royalblue2", "scRNA-seq" = "orange1")
distance_colors <- c("<=50kb" = "#ffe0a8", ">50kb" = "#ffa200")
classic <- theme_classic(base_size = 8) +
  theme(text = element_text(family = "Helvetica"), axis.text = element_text(color = "black"))

# A. Per-sample genotype concordance before and after imputation.
concordance_dirs <- c("arraytopmed_vs_scrnamonopogen", "arraytopmed_vs_scrnatopmed_new")
concordance <- rbindlist(lapply(concordance_dirs, function(d) {
  fread(file.path(benchmark, d, "accuracy.sample_accuracy.tsv"), select = "overall_concordance")
}), idcol = "stage")
concordance[, stage := factor(stage, levels = 1:2, labels = c("Before\nimputation", "After\nimputation"))]
p_a <- ggplot(concordance, aes(stage, overall_concordance)) +
  geom_violin(fill = "#00b0dc") + geom_boxplot(width = 0.05, fill = "white") +
  scale_y_continuous(limits = c(0, 1)) + classic +
  labs(x = "scRNA-seq-derived genotype", y = "Per-sample concordance") +
  theme(plot.margin = margin(15, 5, 15, 5))

# B. Match 1KGP imputation INFO to dosage r² by chromosome and position.
info_dir <- file.path(benchmark, "arraytopmed_vs_scrna1kgp")
info <- fread(file.path(info_dir, "scrna1kgp_imputation_summary.tsv"),
              select = c("CHROM", "POS", "INFO_SCORE"))[INFO_SCORE > 0.3]
info[, CHROM := sub("^chr", "", CHROM)]
accuracy <- fread(file.path(info_dir, "accuracy.variant_accuracy.tsv.gz"), select = c("chr", "bp", "dosage_r2"))
accuracy[, chr := as.character(chr)]
dosage <- info[accuracy, on = .(CHROM = chr, POS = bp), nomatch = NULL, .(INFO_SCORE, dosage_r2)]
rm(info, accuracy)
p_b <- ggplot(dosage, aes(INFO_SCORE, dosage_r2)) +
  geom_hex(bins = 500, size = 0.3) + classic +
  scale_fill_viridis_c(trans = "log10", name = "Variant\ncount") +
  coord_cartesian(ylim = c(0, 1)) + labs(x = "Imputation INFO score", y = expression("Dosage " * r^2)) +
  theme(legend.position = "top", legend.direction = "horizontal", legend.box.margin = margin(0, 0, 0, 0)) +
  guides(fill = guide_colorbar(barwidth = grid::unit(18, "mm"), barheight = grid::unit(2, "mm"),
                               ticks.linewidth = 0.3, frame.linewidth = 0.3))

# C. Concordance of pairwise LD within loci centered on scRNA-seq lead variants.
ld <- fread(file.path(benchmark, "arraytopmed_vs_scrnatopmed_new", "ld_compare_scRNA.tsv"))
ld <- ld[status == "ok" & is.finite(pearson_r2_matrix_correlation) & is.finite(mean_abs_pairwise_r2_difference)]
median_r <- median(ld$pearson_r2_matrix_correlation)
p_c <- ggplot(ld, aes(pearson_r2_matrix_correlation)) +
  geom_histogram(bins = 50, fill = "royalblue4", color = "white") +
  geom_vline(xintercept = median_r, linetype = "dashed", color = "darkred", linewidth = 0.4) +
  annotate("text", x = -Inf, y = Inf, vjust = 1.8, hjust = -0.2, size = 2.5, family = "Helvetica",
           label = paste0("#Loci = ", scales::comma(nrow(ld)))) +
  annotate("text", x = -Inf, y = Inf, vjust = 1.7, hjust = -0.05, size = 2.5, family = "Helvetica",
           label = sprintf("\nMedian Pearson r = %.2f", median_r)) +
  classic + labs(x = "Pearson correlation of pairwise LD", y = "Count")

# D–E. Count distinct significant genes or genomic positions across the five cell types.
sig <- eqtls[pval_nominal < cutoff]
sig[, variant := paste(Chr, BP, sep = "-")]
venn <- function(column, title) {
  sets <- lapply(platforms, function(platform) unique(sig[tissue == platform][[column]]))
  ggvenn(setNames(sets, names(platform_colors)), fill_color = unname(platform_colors), auto_scale = TRUE,
         text_size = 2, set_name_color = "white", set_name_size = 0, fill_alpha = 0.8,
         show_percentage = FALSE, stroke_size = 0.2) +
    labs(title = title) + coord_cartesian(clip = "off") +
    theme(plot.margin = margin(l = 20, r = 20), text = element_text(family = "Helvetica"),
          plot.title = element_text(face = "bold", hjust = 0.5, size = 8, margin = margin(b = -10)))
}
legend_panel <- ggplot(data.frame(x = 1:2, group = names(platform_colors)), aes(x, 1)) +
  geom_tile(aes(fill = group), width = 0.18, height = 0.18) +
  geom_text(aes(label = group), nudge_x = 0.12, hjust = 0, size = 2, family = "Helvetica") +
  scale_fill_manual(values = platform_colors) + xlim(0.75, 2.8) + theme_void() + theme(legend.position = "none")
p_de <- ggdraw() +
  draw_label("D", x = 0, y = 0.97, hjust = 0, vjust = 1, fontfamily = "Helvetica", fontface = "bold", size = 12) +
  draw_plot(venn("phenotype_id", "eGene"), x = 0.15, y = 0.6, width = 0.8, height = 0.35) +
  draw_label("E", x = 0, y = 0.50, hjust = 0, vjust = 1, fontfamily = "Helvetica", fontface = "bold", size = 12) +
  draw_plot(venn("variant", "eVariant"), x = 0.15, y = 0.15, width = 0.8, height = 0.35) +
  draw_plot(legend_panel, x = 0.15, y = 0.1, width = 0.8, height = 0.05)

# F–G. Calculate Pearson r before applying the display limits. Shared pairs pass P < 5e-8 in both analyses.
effect_row <- function(lead_only) {
  plots <- lapply(names(cell_labels), function(ct) {
    x <- matched_effects(eqtls, ct, lead_only)
    x <- x[!is.na(slope_array) & !is.na(slope_scrna) & pval_nominal_scrna < cutoff]
    x[, color := fifelse(pval_nominal_array < cutoff, cell_colors[[ct]], "grey")]
    correlations <- function(d) if (nrow(d) >= 2L) cor(d$slope_array, d$slope_scrna) else NA_real_
    label <- sprintf("scRNA r = %.3f\nShared r = %.3f", correlations(x), correlations(x[pval_nominal_array < cutoff]))
    ggplot(x, aes(slope_array, slope_scrna, color = color)) +
      geom_point(alpha = 0.8, size = 1.2) + scale_color_identity() +
      geom_hline(yintercept = 0, linewidth = 0.5) + geom_vline(xintercept = 0, linewidth = 0.5) +
      scale_x_continuous(limits = c(-2, 2)) + scale_y_continuous(limits = c(-1.5, 1.5)) +
      annotate("text", x = -Inf, y = Inf, label = label, hjust = -0.05, vjust = 1.1, size = 2) +
      labs(x = "Effect Size (SNP-array)", y = "Effect Size (scRNA)", title = cell_labels[[ct]]) +
      theme_classic() +
      theme(legend.position = "none", plot.margin = margin(3, 3, 3, 3), axis.line = element_blank(),
            axis.title = element_text(size = 8), axis.text = element_text(size = 8, color = "black"),
            axis.ticks = element_blank(), panel.grid.major = element_line(),
            plot.title = element_text(size = 8, hjust = 0.5, face = "bold"))
  })
  wrap_plots(plots, nrow = 1) + plot_layout(axis_titles = "collect")
}
p_f <- wrap_elements(full = effect_row(FALSE)) + labs(tag = "F")
p_g <- wrap_elements(full = effect_row(TRUE)) + labs(tag = "G")

# H. Distance to TSS of each platform's significant lead variant per cell type and gene.
leads <- sig[order(pval_nominal), .SD[1L], by = .(tissue, cell_type, phenotype_id)]
leads[, platform := factor(tissue, levels = platforms, labels = names(platform_colors))]
p_h <- ggplot(leads, aes(start_distance / 1e6, color = platform)) +
  geom_density(linewidth = 0.8, alpha = 0.2) + facet_wrap(~cell_factor(cell_type), nrow = 1) +
  classic + scale_color_manual(values = platform_colors) + labs(x = "Distance to TSS (Mb)", y = "Density") +
  theme(legend.position = "none", strip.background = element_rect(color = "white"),
        strip.text = element_text(size = 8, face = "bold"), plot.margin = margin(t = 5, b = 10))

# I. Fraction of array eGenes recovered in scRNA-seq, by distance of the array lead to the TSS.
array_leads <- leads[tissue == "PeruvianWGS",
                     .(cell_type, phenotype_id, variant_array = variant_id, start_distance)]
array_leads[, distance := factor(fifelse(abs(start_distance) > 50000, ">50kb", "<=50kb"),
                                 levels = names(distance_colors))]
array_leads[, recovered := FALSE]
array_leads[leads[tissue == "PeruvianscRNA"], on = .(cell_type, phenotype_id), recovered := TRUE]
recovery <- array_leads[, .(rate = mean(recovered)), by = .(cell_type, distance)]
bottom_theme <- theme(plot.margin = margin(b = 10, l = 10), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
p_i <- ggplot(recovery, aes(cell_factor(cell_type), rate, fill = distance)) +
  geom_col(position = "dodge", width = 0.5) + scale_fill_manual(values = distance_colors) +
  scale_y_continuous(labels = scales::label_percent(), limits = c(0, 1)) + classic + bottom_theme +
  theme(legend.position = c(0.5, 0.9), legend.direction = "horizontal") +
  labs(x = "Cell type", y = "eGene recovery rate", fill = "Distance to TSS")

# J. Compare the array lead to the best available scRNA-seq variant, including subthreshold hits.
scrna_leads <- eqtls[tissue == "PeruvianscRNA"][order(pval_nominal), .(variant_scrna = variant_id[1L]),
                                                by = .(cell_type, phenotype_id)]
bfile <- file.path(root, "PeruvianWGS", "sc-eQTL_R", "genotype", "PeruvianWGS_AMR_indivpruned")
ld_pairs <- lead_ld(merge(array_leads, scrna_leads, by = c("cell_type", "phenotype_id")), bfile)
p_j <- ggplot(ld_pairs[!is.na(r2)], aes(cell_factor(cell_type), r2, fill = distance)) +
  geom_violin(position = position_dodge(width = 0.9), color = "black", alpha = 0.8, scale = "width") +
  geom_boxplot(position = position_dodge(width = 0.9), width = 0.15, color = "black", alpha = 0.7, outlier.size = 0.8) +
  scale_fill_manual(values = distance_colors) + classic + bottom_theme + theme(legend.position = "none") +
  labs(x = "Cell type", y = expression(r^2 ~ "between lead variants"))

tag_theme <- theme(plot.tag = element_text(family = "Helvetica", size = 12, face = "bold"), plot.tag.position = c(0, 1))
top <- ((p_a + labs(tag = "A")) | (p_b + labs(tag = "B")) | (p_c + labs(tag = "C")) | wrap_elements(full = p_de)) +
  plot_layout(widths = c(0.8, 1, 1.15, 1.2)) & tag_theme
middle <- p_f / p_g
bottom <- (p_i + labs(tag = "I")) | (p_j + labs(tag = "J"))
figure <- top / middle / (p_h + labs(tag = "H")) / bottom +
  plot_layout(heights = c(0.8, 2, 0.5, 0.8)) & tag_theme

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(out_dir, "figure_S23_reduce.png"), plot = figure, width = 8, height = 9, units = "in", dpi = 300)
