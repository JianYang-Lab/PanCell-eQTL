#!/usr/bin/env Rscript
# CIMA Figure R4, panels A–E. Keep CIMA_egene_recovery.R beside this script.
# Usage: Rscript CIMA_benchmark_finalplot_compact.R [project_root] [output_dir]

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
source(file.path(dirname(normalizePath(script_file)), "CIMA_egene_recovery.R"))
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggvenn)
  library(patchwork)
})
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[1] else "."
benchmark_dir <- file.path(root, "CIMA_benchmark")
input <- file.path(benchmark_dir, "eqtl")
out <- if (length(args) >= 2L) args[2] else file.path(benchmark_dir, "plots", "compact")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
group_colors <- c(wgs = "#009E73", array = "#4169E1", scrna = "#FFA500")
comparisons <- c("array", "scrna")
read_table <- function(name) fread(file.path(input, name))
if (file.exists(file.path(input, "PREPARATION_INCOMPLETE")))
  stop("Benchmark preparation is still incomplete. Run this script after it finishes.")
cell_types <- read_table("cell_types.tsv")[order(order), cell_type]
availability <- read_table("input_availability.tsv")[cell_type %in% cell_types]
settings <- read_table("benchmark_settings.tsv")
genotypes_skipped <- identical(settings$genotype_calculation, "skipped")
coverage <- availability[, .(n_available = sum(status == "available")), by = group]
legend_labels <- c(wgs = "WGS", array = "SNP-array", scrna = "scRNA-seq")
font <- "Arial"
page_width <- 7.2
theme_set(theme_classic(base_size = 8, base_family = font) + theme(
  text = element_text(color = "black"), axis.text = element_text(size = 7, color = "black"),
  axis.title = element_text(size = 8), axis.line = element_line(linewidth = .25),
  axis.ticks = element_line(linewidth = .25), axis.ticks.length = grid::unit(2, "pt"),
  strip.background = element_blank(), strip.text = element_text(size = 7, face = "bold"),
  plot.title = element_text(size = 8, face = "bold", hjust = .5),
  plot.subtitle = element_text(size = 6.5, hjust = .5),
  plot.tag = element_text(size = 11, face = "bold"), plot.tag.position = "topleft",
  plot.caption = element_text(size = 6.5, hjust = 0),
  plot.margin = margin(4, 4, 4, 4), panel.spacing = grid::unit(5, "pt"),
  legend.position = "top", legend.title = element_blank(), legend.text = element_text(size = 7),
  legend.key.size = grid::unit(9, "pt"), legend.margin = margin(0, 0, 1, 0)))

observed_limits <- function(x, lower_bound, step, minimum_span) {
  x <- x[is.finite(x)]
  if (!length(x)) return(c(lower_bound, 1))
  lower <- max(lower_bound, floor((min(x) - .02 * minimum_span) / step) * step)
  lower <- min(lower, 1 - minimum_span)
  c(lower, 1)
}
blank_panel <- function(label, title, tag) ggplot() +
  annotate("text", x = 0, y = 0, label = label, family = font, size = 2.5) +
  theme_void(base_family = font) + labs(title = title, tag = tag)

samples <- read_table("sample_concordance.tsv")[group %in% comparisons & is.finite(overall_concordance)]
geno_summary <- read_table("genotype_concordance_summary.tsv")
samples[, group := factor(group, levels = comparisons)]
geno_stats <- samples[, .(n = .N, varying = .N >= 3L && sd(overall_concordance) > 0), by = group]
geno_missing <- data.table(group = factor(setdiff(comparisons, as.character(geno_stats$group)),
                                         levels = comparisons))
geno_missing[, label := if (genotypes_skipped) "Not calculated" else
  fifelse(as.character(group) %in% geno_summary$group, "No shared calls", "Pending")]
geno_limits <- c(0, 1)
p_a <- ggplot(samples, aes(group, overall_concordance, fill = group)) +
  geom_violin(data = samples[group %in% geno_stats[varying == TRUE, group]],
              width = .75, linewidth = .25, trim = TRUE) +
  geom_boxplot(width = .15, linewidth = .25, outlier.size = .3) +
  geom_text(data = geno_missing, aes(group, mean(geno_limits), label = label),
             inherit.aes = FALSE, angle = 90, size = 2.1) +
  scale_fill_manual(values = group_colors[comparisons], limits = comparisons, guide = "none") +
  scale_x_discrete(drop = FALSE, labels = legend_labels[comparisons]) +
  scale_y_continuous(breaks = seq(0, 1, .25), expand = expansion(mult = 0)) +
  coord_cartesian(ylim = geno_limits) +
  labs(x = NULL, y = "Concordance with WGS", title = "Genotype concordance", tag = "A") +
  theme(axis.text.x = element_text(size = 6.5))

egenes <- read_table("significant_eGenes.tsv")
evariants <- read_table("significant_eVariants.tsv.gz")
complete_groups <- names(group_colors)[names(group_colors) %in% coverage[n_available == length(cell_types), group]]
venns <- list()
for (entity in c("eGene", "eVariant")) {
  x <- if (entity == "eGene") egenes else evariants
  id <- if (entity == "eGene") "phenotype_id" else "variant_id"
  sets <- setNames(lapply(complete_groups, function(g) unique(x[cell_type %in% cell_types & group == g, get(id)])),
                   complete_groups)
  tag <- if (entity == "eGene") "B" else "C"
  if (length(sets) < 2L || sum(lengths(sets)) == 0L) {
    venns[[entity]] <- blank_panel("No complete overlap estimate", paste0(entity, " overlap"), tag)
  } else {
    venns[[entity]] <- ggvenn(sets, fill_color = unname(group_colors[complete_groups]),
      auto_scale = FALSE, text_size = 2.05, set_name_size = 0,
      show_percentage = FALSE, fill_alpha = .8, stroke_size = .25) +
      labs(title = paste0(entity, " overlap"), tag = tag,
        subtitle = if (length(sets) == 3L) NULL else "scRNA-seq overlap pending") +
      theme(text = element_text(family = font, color = "black"),
        plot.title = element_text(size = 8, face = "bold", hjust = .5),
        plot.subtitle = element_text(size = 6.5, hjust = .5),
        plot.tag = element_text(size = 11, face = "bold"), plot.margin = margin(4, 1, 4, 1))
  }
}

# Every observation is a saved, ALT-aligned Pearson r for one cell type.
# Never pool gene-variant pairs or weight cell types by their pair counts.
correlations <- read_table("effect_correlations.tsv")[cell_type %in% cell_types &
  comparison %in% comparisons & panel %in% c("all_pairs", "lead_union") & subset %in% c("Query", "Shared")]
if (anyDuplicated(correlations[, .(cell_type, comparison, panel, subset)]))
  stop("Duplicate cell-type correlation records in effect_correlations.tsv")
if (any(is.finite(correlations$pearson_r) & abs(correlations$pearson_r) > 1 + 1e-10))
  stop("Pearson correlations must lie in [-1, 1].")
correlations[, valid := comparison_available & is.finite(pearson_r)]
correlations[valid == FALSE, pearson_r := NA_real_]
correlations[, `:=`(comparison = factor(comparison, levels = comparisons),
  panel = factor(panel, levels = c("all_pairs", "lead_union")),
  subset = factor(subset, levels = c("Query", "Shared")))]
corr_summary <- correlations[, {
  r <- pearson_r[valid]
  .(n_cell_types = length(r),
    varying = length(r) >= 3L && sd(r) > 0)
}, by = .(panel, subset, comparison)]
corr_limits <- observed_limits(correlations$pearson_r, -1, .1, .2)
correlation_plot <- function(panel_name, tag, title) {
  x <- correlations[panel == panel_name & valid == TRUE]
  summary <- corr_summary[panel == panel_name]
  density_data <- merge(x, summary[varying == TRUE, .(subset, comparison)], by = c("subset", "comparison"))
  ggplot(x, aes(comparison, pearson_r, fill = comparison)) +
    geom_violin(data = density_data, width = .8, trim = TRUE, linewidth = .25, alpha = .65) +
    geom_boxplot(width = .17, linewidth = .25, outlier.shape = NA, fill = "white") +
    geom_point(position = position_jitter(width = .07, height = 0, seed = 20260928),
               size = .65, alpha = .75, shape = 21, stroke = .12) +
    geom_text(data = summary, aes(comparison, -Inf, label = paste0("n = ", n_cell_types)),
               inherit.aes = FALSE, vjust = -.5, size = 2.0, family = font) +
    scale_fill_manual(values = group_colors[comparisons], limits = comparisons, guide = "none") +
    scale_x_discrete(drop = FALSE, labels = c(array = "SNP-array", scrna  = "scRNA-seq")) +
    scale_y_continuous(breaks = pretty(corr_limits, n = 3)) +
    coord_cartesian(ylim = corr_limits) +
    facet_grid(cols = vars(subset), drop = FALSE,
      labeller = labeller(subset = c(Query = "Comparator-significant", Shared = "Shared-significant"))) +
    labs(x = NULL, y = "Pearson r with WGS", title = title, tag = tag) +
    theme(axis.text.x = element_text(size = 6.5), strip.text = element_text(size = 6.3),
          panel.grid.major.y = element_line(color = "#EEEEEE", linewidth = .2))
}
p_d <- correlation_plot("all_pairs", "D", "Effect correlation: all pairs")
# Both recovery estimates use the same 12 memory T-cell types and WGS reference.
# The x-axis retains total processed expression cells per final scRNA donor.
recovery_data <- cima_egene_recovery_data(root, 5e-8)
recovery_points <- copy(recovery_data$summary)
recovery_points[, label := paste0(
  fifelse(cohort == "CIMAscRNAsub", "Downsampled to Peruvian\n", ""),
  sprintf("%.2f%%", egene_recovery_percent))]
recovery_points[, `:=`(label_hjust = fifelse(cohort == "CIMAscRNA", 0, 1),
                   label_vjust = fifelse(cohort == "CIMAscRNA", 1.6, 1.2))]
p_e <- ggplot(recovery_points, aes(mean_cells_per_donor, egene_recovery_percent)) +
  geom_hline(yintercept = 100, color = "#BBBBBB", linewidth = .25, linetype = "dashed") +
  geom_point(aes(shape = cohort), color = group_colors[["scrna"]], size = 2.5) +
  geom_text(aes(label = label, hjust = label_hjust, vjust = label_vjust),
            size = 2.3, family = font, lineheight = 1.1) +
  scale_shape_manual(values = c(CIMAscRNA = 16, CIMAscRNAsub = 17), guide = "none") +
  scale_x_reverse(limits = c(max(recovery_points$mean_cells_per_donor) * 1.12, 0),
                     labels = function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE),
                     breaks = function(x) pretty(x, n = 4), expand = expansion(mult = 0)) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25),
                     expand = expansion(mult = c(0, .06))) +
  labs(x = "Mean total cells per donor", y = "WGS eGene recovery (%)",
       title = "eGene recovery in 12 memory T cells", tag = "E") +
  theme(panel.grid.major.y = element_line(color = "#EEEEEE", linewidth = .2))

# Measure the shared source legend on the same device used in the original layout.
legend_plot <- ggplot(data.table(group = names(group_colors)), aes(group, 1, fill = group)) +
  geom_col() + scale_fill_manual(values = group_colors, limits = names(group_colors), labels = legend_labels)
top_row <- (p_a | venns$eGene | venns$eVariant) + plot_layout(widths = c(.65, 1, 1))
make_source_legend <- function() {
  cairo_pdf(tempfile(fileext = ".pdf"), width = page_width, height = 1)
  legend_table <- ggplotGrob(legend_plot)
  invisible(dev.off())
  guide_indices <- grep("^guide-box", legend_table$layout$name)
  guide_indices <- guide_indices[vapply(legend_table$grobs[guide_indices],
    function(g) inherits(g, "gtable"), logical(1))]
  stopifnot(length(guide_indices) == 1L)
  wrap_elements(full = legend_table$grobs[[guide_indices]])
}

figure_ae <- (make_source_legend() / top_row / (p_d | p_e)) +
  plot_layout(heights = c(.2, 1.65, 1.9))
figure_ae <- figure_ae & theme(legend.position = "none")
ggsave(file.path(out, "figure_R4_compact_A-E.png"), figure_ae, width = page_width, height = 4.6,
       units = "in", dpi = 300, bg = "white")
