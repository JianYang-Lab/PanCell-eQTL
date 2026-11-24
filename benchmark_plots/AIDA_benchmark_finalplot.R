#!/usr/bin/env Rscript
# AIDA and MixedEAS eGene discovery and shared significant eQTL effects
# Usage: Rscript AIDA_benchmark_finalplot.R [project_root] [output_dir]

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggvenn)
  library(cowplot)
  library(patchwork)
})
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[1] else "."
out <- if (length(args) >= 2L) args[2] else file.path(root, "crossTissue_analysis_R", "AIDA_benchmark")
results <- file.path(root, "Blood", "sc-eQTL_R", "results")
ct_list <- c("CD4TNC","CD4TEM","CD4TEMRA","CD8TEMRA","NKp","Plasma","pDC")
common_theme <- theme_classic(base_size = 12) +
    theme(axis.text=element_text(color="black"), text=element_text(family="Helvetica"),plot.title=element_text(hjust=0.5))

metadata <- fread(file.path(root, "crossTissue_analysis_R", "tissue_ancestry_celltype_ss_cc_colors.txt"))
metadata <- unique(metadata[tissue == "Blood" & sample_size >= 50 & total_count >= 10000 & cell_type %in% ct_list,
                            .(cell_type, cell_type_short, color)])
cell_labels <- setNames(metadata$cell_type_short, metadata$cell_type)
cell_type_color_dict <- setNames(metadata$color, metadata$cell_type_short)

# Keep pairs with nominal P < 5e-8 in each cohort. A stricter export threshold requires the nominal files.
cohorts <- c("AIDA", "sub5")
sig_list_all <- rbindlist(lapply(ct_list, function(ct) {
  rbindlist(lapply(cohorts, function(grp) {
    anc <- if (grp == "AIDA") "EASAIDA" else paste0("EAS", ct, "sub5")
    prefix <- paste0(ct, grp, "_", anc, "_pc5")
    sig <- fread(file.path(results, "sig", paste0(prefix, "_sig.tsv.gz")),
                 select = c("phenotype_id", "variant_id", "pval_nominal", "pval_nominal_threshold", "slope"))
    if (!nrow(sig) || any(sig$pval_nominal_threshold < 5e-8, na.rm = TRUE)) {
      files <- file.path(results, "raw", paste0(prefix, ".cis_qtl_pairs.", 1:22, ".parquet"))
      sig <- arrow::open_dataset(files, format = "parquet") |>
        dplyr::filter(pval_nominal < 5e-8) |>
        dplyr::select(phenotype_id, variant_id, slope) |> dplyr::collect() |> as.data.table()
    } else {
      sig <- sig[pval_nominal < 5e-8, .(phenotype_id, variant_id, slope)]
    }
    sig[, `:=`(group = grp, cell_type = ct, cell_type_short = cell_labels[[ct]])]
    sig
  }))
}))
smr_list_ts <- sig_list_all[, .(neGenes = uniqueN(phenotype_id)), by = .(cell_type_short, group)]
smr_list_ts <- smr_list_ts[CJ(cell_type_short = unname(cell_labels[ct_list]), group = cohorts),
                           on = .(cell_type_short, group)]
smr_list_ts[is.na(neGenes), neGenes := 0L]

# eGene count ----------------------------------------
p1 <- smr_list_ts %>% dplyr::mutate(cell_type=factor(cell_type_short, levels=unname(cell_labels[ct_list])), group=factor(group,levels=c("AIDA","sub5"), labels=c("AIDA","MixedEAS"))) %>%
 ggplot(aes(x=cell_type,y=neGenes, group=group, fill=group)) +
    geom_bar(stat="identity", position="dodge", width=0.7) +
    scale_fill_manual(values = c(AIDA="#2646af",MixedEAS="#b0ccff")) +
    common_theme +
    theme(legend.position="none",
        axis.text=element_text(size=12,color="black"),
        axis.text.x=element_text(angle=45,hjust=1,vjust=1),
        axis.title=element_text(size=12,color="black"), plot.margin=margin(20,20,20,20)) +
    labs(x="Cell type", y="Number of eGenes", fill="")
p2 <- ggvenn(
        list("AIDA"=unique(sig_list_all[group=="AIDA"]$phenotype_id),
            "MixedEAS"=unique(sig_list_all[group=="sub5"]$phenotype_id)
           ),  auto_scale=TRUE,
        set_name_size=0,fill_color=c("#2646af","#b0ccff"),fill_alpha=0.8, show_percentage=FALSE
        ) + theme(plot.title = element_text(hjust=0.5,size=18))

set_names <- c("AIDA", "MixedEAS")
legend_plot <- ggplot(data.frame(x=1:2, y=c(1,1), Group=factor(set_names, levels=c("AIDA","MixedEAS"))), aes(x, y, fill=Group)) +
    geom_tile() +
    scale_fill_manual(values = c("AIDA" = "#2646af", "MixedEAS" = "#b0ccff")) +
    theme_minimal() +
    theme(text=element_text(family="Helvetica"),legend.position = "right",
          legend.title = element_blank(),
          legend.text = element_text(size = 12))

shared_legend <- get_legend(legend_plot)

# eQTL correlation ------------------------------
shared_eqtl <- function(sig_list_all, celltype, group1, group2){
    df1 <- sig_list_all[cell_type==celltype & group==group1, .(phenotype_id,variant_id,slope,cell_type_short)]
    df2 <- sig_list_all[cell_type==celltype & group==group2, .(phenotype_id,variant_id,slope,cell_type_short)]
    merged_df <- merge(df1, df2, by = c("phenotype_id","variant_id","cell_type_short"), all.x=FALSE, suffixes=c(paste0("_",group1),paste0("_",group2)))
    return(merged_df)
}

plot_slope_comparison <- function(data, x_col, y_col, x_label, y_label,title, celltype_color, xmax=2, xmin=-2, ymax=1.5, ymin=-1.5) {
  plot_data <- data[!is.na(data[[x_col]]) & !is.na(data[[y_col]]), ]

  if(nrow(plot_data)>0){
  # Calculate Pearson correlation coefficient
  r_val <- cor(plot_data[[x_col]], plot_data[[y_col]], method = "pearson")

  x_pos <- -1.9; y_pos <- 1.5

  p <- ggplot(plot_data, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_point(alpha = 0.8, color = celltype_color, size = 1.5) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
    scale_x_continuous(limits=c(xmin,xmax)) +
    scale_y_continuous(limits=c(ymin,ymax))+
    annotate("text", x = x_pos, y = y_pos,
             label = sprintf("r = %.3f", r_val),
             hjust = 0, vjust = 1, size = 3, fontface = "bold") +
    labs(x = x_label, y = y_label, title = title) +
    theme_classic() +
    theme(
        text=element_text(family="Helvetica"),
        axis.line = element_blank(),
        axis.title = element_text(size = 10,face="bold"),
        axis.text = element_text(size = 10,color="black"),
        axis.ticks = element_blank(),
        panel.grid.major = element_line(color="#bababa5a"),
      plot.title = element_text(size=12, hjust=0.5, face = "bold")
    )
    return(p)
  }else{
    # Return empty plot
   empty_p <- ggplot() +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50", linewidth = 0.8) +
      annotate("text", x = 0, y = 0, label = "No shared eQTLs",
               size = 5, fontface = "italic", color = "gray40") +
      coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1)) +
      labs(x = x_label, y = y_label, title = title) +
      theme_minimal() +
      theme(
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10,color="black"),
      plot.title = element_text(size=12, hjust=0.5, face = "bold")
      )
    return(empty_p)
  }

}

merged_df_list <- list()
group1="AIDA"
group2="sub5"
for (ct in ct_list){
    merged_df_list[[ct]] <- shared_eqtl(sig_list_all, ct, group1, group2)
}

p_list <- lapply(ct_list, function(ct) {
    tmp <- merged_df_list[[ct]]
    plot_slope_comparison(tmp,
                                x_col = "slope_AIDA", y_col = "slope_sub5",
                                x_label = "Effect Size (AIDA)", y_label = "Effect Size (MixedEAS)", title = cell_labels[[ct]], celltype_color=cell_type_color_dict[cell_labels[[ct]]])
    })

p3 <-  wrap_plots(p_list,nrow=1) + plot_layout(axes="collect")

legend_panel <- ggdraw() +
  draw_grob(
    shared_legend,
    x = 0.5,
    y = 0.5,
    hjust = 0.5,
    vjust = 0.5
  )

p_top <- (
  (p1 + labs(tag = "A")) |
  (p2 + labs(tag = "B")) |
  wrap_elements(full = legend_panel)
) +
  plot_layout(widths = c(3, 1.5, 0.7))

figure <- p_top / (wrap_elements(full=p3) + labs(tag="C")) + plot_layout(height=c(1,1)) &
  theme(
    plot.tag = element_text(
      family = "Helvetica",
      size = 18,
      face = "bold"
    ),
    plot.tag.position = c(0,1)
  )

dir.create(out, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(out, "FigureR_AIDA.png"), plot = figure, width = 12, height = 6, dpi = 300)
