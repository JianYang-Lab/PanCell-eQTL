# Usage: Rscript healthy_disease_public.R [project_root] [output_dir]
# Plot Lung–EUR eQTL effect sizes in mixed and healthy-only samples.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) args[1] else "."
analysis_dir <- file.path(root, "crossTissue_analysis_R")
figure_dir <- file.path(analysis_dir, "healthy_disease")
out_dir <- if (length(args) >= 2) args[2] else figure_dir

# Each list entry contains mixed-sample eQTLs (P < 5e-8) matched to healthy-only
# results by gene and variant. Preserve the cell type and point order in the list.
pairs <- readRDS(file.path(figure_dir, "merged_df_list.rds"))
metadata <- fread(file.path(analysis_dir, "tissue_ancestry_celltype_ss_cc_colors.txt"))
eqtl_summary <- fread(file.path(root, "Lung", "sc-eQTL_R", "results", "sig", "eQTL_summary_5en8.tsv"))
cells <- metadata[tissue == "Lung" & ancestry == "EUR" & sample_size >= 50 & total_count >= 10000 &
                  cell_type %in% eqtl_summary[ancestry == "EUR", cell_type]]
cell_types <- intersect(names(pairs), cells$cell_type_short)
colors <- setNames(cells$color, cells$cell_type_short)

plot_cell_type <- function(cell_type) {
  tab <- pairs[[cell_type]]
  healthy <- paste0(cell_type, "Healthy")
  p_healthy <- tab[[paste0("pval_nominal_", healthy)]]
  df <- data.frame(mixed = tab[[paste0("slope_", cell_type)]],
                   healthy = tab[[paste0("slope_", healthy)]],
                   significant = !is.na(p_healthy) & p_healthy < 5e-8)
  df <- df[complete.cases(df$mixed, df$healthy), ]
  if (!nrow(df)) return(NULL)

  r <- if (nrow(df) >= 2 && sd(df$mixed) > 0 && sd(df$healthy) > 0) {
    cor(df$mixed, df$healthy, method = "pearson")
  } else NA_real_
  xr <- range(df$mixed)
  yr <- range(df$healthy)
  xpad <- if (diff(xr) == 0) 0.05 else 0.05 * diff(xr)
  ypad <- if (diff(yr) == 0) 0.05 else 0.05 * diff(yr)

  ggplot(df, aes(mixed, healthy, color = significant)) +
    geom_hline(yintercept = 0, color = "grey85", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "grey85", linewidth = 0.3) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50", linewidth = 0.5) +
    geom_point(alpha = 0.5, size = 0.5) +
    scale_x_continuous(breaks = scales::breaks_width(1), labels = scales::label_number(accuracy = 1)) +
    scale_y_continuous(breaks = scales::breaks_width(1), labels = scales::label_number(accuracy = 1)) +
    scale_color_manual(values = c("FALSE" = "grey75", "TRUE" = colors[[cell_type]])) +
    annotate("text", x = xr[1] + xpad, y = yr[2] - ypad, label = sprintf("r = %.3f", r),
             hjust = 0, vjust = 1, size = 2, fontface = "bold", family = "Helvetica") +
    coord_cartesian(xlim = xr + c(-xpad, xpad), ylim = yr + c(-ypad, ypad)) +
    labs(x = "Effect Size (Mixed samples)", y = "Effect Size (Healthy samples)", title = cell_type) +
    theme_classic(base_size = 7, base_family = "Helvetica") +
    theme(panel.grid.major = element_line(color = "gray"), panel.grid.minor = element_blank(),
          axis.title = element_text(size = 7, color = "black", face = "bold"),
          axis.text = element_text(size = 5.5, color = "black"),
          plot.title = element_text(size = 6, hjust = 0.5, face = "bold"),
          axis.line = element_blank(), axis.ticks = element_blank(),
          legend.position = "none", plot.margin = margin(5, 5, 5, 5))
}

plots <- Filter(Negate(is.null), lapply(cell_types, plot_cell_type))
figure <- wrap_plots(plots, nrow = 4, axis_titles = "collect")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(out_dir, "effect_size_compare_all_reduce.png"), figure,
       width = 5.8, height = 3.8, units = "in", dpi = 300)
