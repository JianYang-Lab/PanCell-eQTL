#!/usr/bin/env Rscript
# OneK1K figure, from the prepared OneK1K_benchmark/eqtl tables.
# Usage: Rscript OneK1K_benchmark_finalplot_compact.R [project_root] [output_dir]

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})
setDTthreads(4)
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[1] else "."
input <- file.path(root, "OneK1K_benchmark/eqtl")
out <- if (length(args) >= 2L) args[2] else file.path(root, "OneK1K_benchmark/plots/compact")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
cell_types <- c("BIN", "BMemNS", "BMemS", "Plasma", "CD4TNC", "CD4TEM",
                "CD4TEMRA", "Treg", "CD8TNC", "CD8TRM", "CD8TEMRA", "MAIT",
                "gdT", "NKp", "MonoC", "MonoNC", "DC2", "pDC")
group_colors <- c(array = "#4169E1", scrna = "#FFA500")
group_labels <- c(array = "SNP-array", scrna = "scRNA-seq")
font <- "Arial"
page_width <- 8.6
page_height <- 2.65
read_table <- function(name) fread(file.path(input, name))
required <- c("sample_concordance.tsv",
  "significant_eGenes.tsv", "significant_eVariants.tsv.gz", "effect_correlations.tsv")
missing <- required[!file.exists(file.path(input, required))]
if (length(missing)) stop("Missing prepared benchmark tables: ", paste(missing, collapse = ", "))

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

# Panel A retains the before/after comparison and its stage-specific SNP sets.
samples <- read_table("sample_concordance.tsv")[stage %in% c("before", "after")]
if (anyDuplicated(samples[, .(stage, sample_id)])) stop("Duplicate donor concordance records")
samples[, stage := factor(stage, c("before", "after"))]
p_a <- ggplot(samples[is.finite(overall_concordance)], aes(stage, overall_concordance)) +
  geom_violin(fill = group_colors[["scrna"]], width = .75, linewidth = .25, trim = TRUE) +
  geom_boxplot(width = .15, fill = "white", linewidth = .25, outlier.size = .3) +
  scale_x_discrete(labels = c(before = "Before\nimputation",
                              after = "After imputation\nand filtering"), drop = FALSE) +
  scale_y_continuous(breaks = seq(0, 1, .25), expand = expansion(mult = 0)) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = NULL, y = "Concordance with SNP-array",
       title = "Genotype concordance", tag = "A") +
  theme(axis.text.x = element_text(size = 6.5))

# Two-circle geometry with areas proportional to the set and intersection counts.
# The center separation is solved from the exact circle-intersection area.
proportional_venn <- function(sets, text_size) {
  n <- lengths(sets)
  shared <- length(intersect(sets[[1]], sets[[2]]))
  if (any(n <= 0)) stop("Proportional Venn plots require two nonempty sets")
  radius <- sqrt(n / max(n))
  a <- unname(radius[1])
  b <- unname(radius[2])
  intersection_area <- function(d) {
    if (d >= a + b) return(0)
    if (d <= abs(a - b)) return(pi * min(a, b)^2)
    a^2 * acos(pmax(-1, pmin(1, (d^2 + a^2 - b^2) / (2 * d * a)))) +
      b^2 * acos(pmax(-1, pmin(1, (d^2 + b^2 - a^2) / (2 * d * b)))) -
      .5 * sqrt(pmax(0, (-d + a + b) * (d + a - b) * (d - a + b) * (d + a + b)))
  }
  target <- pi * shared / max(n)
  d <- if (shared == 0) a + b + .1 else if (shared == min(n)) abs(a - b) else
    uniroot(function(z) intersection_area(z) - target, c(abs(a - b), a + b), tol = 1e-12)$root
  theta <- seq(0, 2 * pi, length.out = 721)
  shapes <- rbindlist(list(
    data.table(group = names(sets)[1], x = a * cos(theta), y = a * sin(theta)),
    data.table(group = names(sets)[2], x = d + b * cos(theta), y = b * sin(theta))))
  counts <- data.table(region = c("array_only", "shared", "scrna_only"),
    n = c(n[1] - shared, shared, n[2] - shared),
    x = c((-a + d - b) / 2, (max(-a, d - b) + min(a, d + b)) / 2, (a + d + b) / 2),
    y = 0, hjust = .5)
  segments <- data.table(x = numeric(), y = numeric(), xend = numeric(), yend = numeric())
  # Move a narrow exclusive-region label outside, retaining a leader to the region.
  for (side in c("array_only", "scrna_only")) {
    region_n <- counts[region == side, n]
    width <- if (side == "array_only") d - b + a else d + b - a
    if (region_n > 0 && width < .35) {
      anchor <- counts[region == side, x]
      direction <- if (side == "array_only") -1 else 1
      label_x <- if (side == "array_only") -a - .18 else d + b + .18
      counts[region == side, `:=`(x = label_x, y = .55, hjust = if (direction < 0) 1 else 0)]
      segments <- rbind(segments, data.table(x = label_x, y = .48, xend = anchor, yend = 0))
    }
  }
  counts <- counts[n > 0]
  x_limits <- range(c(shapes$x, counts$x)) + c(-.15, .15)
  if (any(counts$hjust == 0)) x_limits[2] <- x_limits[2] + .2
  if (any(counts$hjust == 1)) x_limits[1] <- x_limits[1] - .2
  p <- ggplot(shapes, aes(x, y, group = group)) +
    geom_polygon(aes(fill = group), alpha = .8) +
    geom_polygon(fill = NA, color = "black", linewidth = .25) +
    geom_segment(data = segments, aes(x, y, xend = xend, yend = yend),
                 inherit.aes = FALSE, linewidth = .25) +
    geom_text(data = counts, aes(x, y, label = n, hjust = hjust), inherit.aes = FALSE,
              size = text_size, family = font) +
    scale_fill_manual(values = group_colors, guide = "none") +
    coord_fixed(xlim = x_limits, ylim = c(-max(radius), max(radius)) + c(-.1, .1),
                expand = FALSE, clip = "off") + theme_void(base_family = font)
  p
}

# Panels B/C aggregate exact identifiers across all 18 cell types.
entities <- list(eGene = read_table("significant_eGenes.tsv"),
                 eVariant = read_table("significant_eVariants.tsv.gz"))
venns <- list()
for (entity in names(entities)) {
  x <- entities[[entity]]
  id <- if (entity == "eGene") "phenotype_id" else "variant_id"
  sets <- setNames(lapply(names(group_colors), function(g)
    unique(x[cell_type %in% cell_types & group == g, get(id)])), names(group_colors))
  venns[[entity]] <- proportional_venn(sets, text_size = if (entity == "eVariant") 1.85 else 2.15) +
    labs(title = paste0(entity, " overlap"), subtitle = "All 18 cell types",
         tag = if (entity == "eGene") "B" else "C") +
    theme(text = element_text(family = font, color = "black"),
      plot.title = element_text(size = 8, face = "bold", hjust = .5),
      plot.subtitle = element_text(size = 6.5, hjust = .5),
      plot.tag = element_text(size = 11, face = "bold"), plot.tag.position = "topleft",
      plot.margin = margin(4, 1, 4, 1))
  venns[[entity]]$coordinates$clip <- "off"
}
rm(entities)

# Panel D uses one saved all-pair ALT-aligned Pearson r per cell type and subset.
correlations <- read_table("effect_correlations.tsv")[cell_type %in% cell_types &
  panel == "all_pairs" & subset %in% c("scRNA", "Shared")]
if (anyDuplicated(correlations[, .(cell_type, panel, subset)]))
  stop("Duplicate cell-type correlation records")
correlations <- merge(CJ(cell_type = cell_types, panel = "all_pairs",
  subset = c("scRNA", "Shared")), correlations,
  by = c("cell_type", "panel", "subset"), all.x = TRUE)
if (any(abs(correlations$pearson_r[is.finite(correlations$pearson_r)]) > 1 + 1e-10))
  stop("Pearson correlations must lie in [-1, 1]")
correlations[, `:=`(valid = is.finite(pearson_r),
  subset = factor(subset, c("scRNA", "Shared")))]
correlations[valid == FALSE, pearson_r := NA_real_]
corr_summary <- correlations[, {
  r <- pearson_r[valid]
  .(n_cell_types = length(r),
    varying = length(r) >= 3L && sd(r) > 0)
}, by = .(panel, subset)]
finite_r <- correlations[valid == TRUE, pearson_r]
corr_limits <- c(if (length(finite_r)) max(-1, min(.8, floor((min(finite_r) - .004) / .1) * .1)) else -1, 1)
correlation_plot <- function(panel_name, tag, title) {
  x <- correlations[panel == panel_name & valid == TRUE]
  summary <- corr_summary[panel == panel_name]
  density_data <- x[subset %in% summary[varying == TRUE, subset]]
  ggplot(x, aes(subset, pearson_r)) +
    geom_violin(data = density_data, fill = group_colors[["scrna"]],
                width = .8, trim = TRUE, linewidth = .25, alpha = .65) +
    geom_boxplot(width = .17, linewidth = .25, outlier.shape = NA, fill = "white") +
    geom_point(position = position_jitter(width = .07, height = 0, seed = 20260929),
               size = .8, alpha = .85, shape = 21, stroke = .15, fill = group_colors[["scrna"]]) +
    geom_text(data = summary, aes(x = subset, y = -Inf, label = paste0("n = ", n_cell_types)),
               inherit.aes = FALSE, vjust = -.5, size = 2, family = font) +
    scale_x_discrete(limits = c("scRNA", "Shared"),
                     labels = c(scRNA = "scRNA-\nsignificant", Shared = "Shared-\nsignificant"), drop = FALSE) +
    scale_y_continuous(breaks = pretty(corr_limits, n = 3)) +
    coord_cartesian(ylim = corr_limits) +
    labs(x = NULL, y = "Pearson r with SNP-array", title = title, tag = tag) +
    theme(axis.text.x = element_text(size = 6.5), strip.text = element_text(size = 6.5),
          panel.grid.major.y = element_line(color = "#EEEEEE", linewidth = .2))
}
p_d <- correlation_plot("all_pairs", "D", "Effect correlation")

# Use the PDF font device when measuring the shared source legend.
legend_plot <- ggplot(data.table(group = names(group_colors)), aes(group, 1, fill = group)) +
  geom_col() + scale_fill_manual(values = group_colors, breaks = names(group_colors), labels = group_labels)
cairo_pdf(tempfile(fileext = ".pdf"), width = page_width, height = 1)
legend_table <- ggplotGrob(legend_plot)
invisible(dev.off())
guide_indices <- grep("^guide-box", legend_table$layout$name)
guide_indices <- guide_indices[vapply(legend_table$grobs[guide_indices],
  function(g) inherits(g, "gtable"), logical(1))]
stopifnot(length(guide_indices) == 1L)
source_legend <- wrap_elements(full = legend_table$grobs[[guide_indices]])
panel_row <- (p_a | venns$eGene | venns$eVariant | p_d) +
  plot_layout(widths = c(.9, 1.05, 1.15, .95))
figure_four <- (source_legend / panel_row) + plot_layout(heights = c(.16, 1))
figure_four <- figure_four & theme(legend.position = "none")
ggsave(file.path(out, "figure_R4_compact.png"), figure_four, width = page_width, height = page_height,
       units = "in", dpi = 300, bg = "white")
