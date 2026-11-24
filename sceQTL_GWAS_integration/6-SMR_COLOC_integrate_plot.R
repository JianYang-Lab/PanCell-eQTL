#!/usr/bin/env Rscript
# Usage: Rscript 5-SMR_COLOC_plot_public.R [project_root] [output_dir] [auto|raw|snapshot]
# Auto uses the archived SMR and COLOC lists for each analysis when both exist.
# Raw reads per-context SMR/output/*.smr and COLOC/output/*.coloc from each tissue.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scatterpie)
})

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[1] else "."
analysis_dir <- file.path(root, "crossTissue_analysis_R")
out <- if (length(args) >= 2) args[2] else file.path(analysis_dir, "SMR_COLOC", "public")
mode <- if (length(args) >= 3) match.arg(args[3], c("auto", "raw", "snapshot")) else "auto"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

blood_traits <- c("Asthma", "AS", "AR", "COPD", "CD", "IBD", "RA", "SLE")
resp_traits <- c("IPF", "COPD", "CPD", "LuC", paste0("COVID19-HGI-", c("A2","B1","B2","C2")),
                 "CAD", "SBP", "DBP", "AFib", "Pleurisy", "Bronchiectasis", "Chronic_bronchitis",
                 "Bronchitis", "Chronic_Sinusitis", "Ped_asthma", "Asthma", "Pneumoconiosis",
                 "Pneumonia", "Pneumothorax", "PTB", "SAS", "Tonsillitis")
cell_order <- list(
  Blood = c("CD4TNC","CD4TNCInf","CD4TEM","Treg","CD4TEMRA","CD8TNC","CD8TRM","CD8TEMRA",
            "MAIT","gdT","NKp","NKn","BIN","BMemNS","BMemS","Plasma","Plasmablasts","MonoC",
            "MonoInf","MonoNC","DC1","DC2","Macro","Nph","pDC","Eryth","HSPC","Platelet"),
  Lung = c("CD4TNC","CD4TEM","Treg","CD8TEMRA","CD8TEMTRM","CD8TRM","Tfh","NKp","NKn","NKCyc",
           "BIN","BMem","Plasma","MonoC","MonoNC","Macro","MacroCyc","MacroInt","MacroAlv","Mast",
           "DC1","DC2","MigDC","pDC","Eryth","Basal","Suprabasal","Ciliated","Club","Goblet","Ion",
           "AT1","AT2","AeroCap","ECArt","ECCap","ECVenPulm","ECVenSys","LECdiff","LECmat",
           "FibroAdv","FibroAlv","Myofibro","Pericyte","Meso")
)
metadata <- fread(file.path(analysis_dir, "tissue_ancestry_celltype_ss_cc_colors.txt"))
tested <- rbindlist(lapply(names(cell_order), function(ts) {
  fread(file.path(root, ts, "sc-eQTL_R/results/sig/eQTL_summary_5en8.tsv"),
        select = c("cell_type", "ancestry"))[, tissue := ts]
}))
labels <- unique(merge(metadata[sample_size >= 50 & total_count >= 10000], tested,
                       by = c("tissue", "ancestry", "cell_type"))[, .(tissue, cell_type, cell_type_short, color)])
stopifnot(!anyDuplicated(labels, by = c("tissue", "cell_type")))
display_order <- function(ts) {
  mapping <- labels[tissue == ts]
  na.omit(mapping$cell_type_short[match(cell_order[[ts]], mapping$cell_type)])
}

read_raw <- function(contexts, method) {
  prefix <- ifelse(contexts$celltype == "JCTF_Japanese", "JCTF_Japanese",
                   paste(contexts$tissue, contexts$celltype, sep = "_"))
  files <- file.path(root, contexts$tissue, "sceqtl_gwas_integration_R", method, "output",
                     paste0(method, "_", contexts$trait, "_", prefix, "_", contexts$ancestry,
                            if (method == "SMR") ".smr" else ".coloc"))
  available <- which(file.size(files) > 1)
  if (!length(available)) stop("No ", method, " input files found for the requested contexts.")
  message(method, ": reading ", length(available), " nonempty files")
  lapply(available, function(i) {
    x <- fread(files[i], sep = "\t", nThread = 1)
    if (!nrow(x)) return(NULL)
    for (field in names(contexts)) set(x, j = field, value = contexts[[field]][i])
    x
  })
}

read_method <- function(contexts, method, snapshot) {
  if (is.null(snapshot)) {
    blocks <- read_raw(contexts, method)
  } else {
    message(method, ": reading ", basename(snapshot))
    blocks <- readRDS(snapshot)
  }
  if (method == "COLOC") {
    # Older headerless COLOC files use the same nine-column layout.
    schema <- c("probe_name","gene_name","pp0","pp1","pp2","pp3","pp4","topsnp","topsnp_pp")
    for (b in blocks) if (all(paste0("V", 1:9) %in% names(b))) setnames(b, paste0("V", 1:9), schema)
  }
  x <- rbindlist(blocks, use.names = TRUE)
  if (!"tissue" %in% names(x)) x[, tissue := "Blood"]
  x <- x[contexts, on = names(contexts), nomatch = NULL]
  if (method == "SMR") {
    x[, `:=`(p_SMR_bhadj = p.adjust(p_SMR, "BH"), p_SMR_bfradj = p.adjust(p_SMR, "bonferroni")),
      by = .(tissue, trait, celltype, ancestry)]
  } else setnames(x, "gene_name", "Gene")
  x <- merge(x, labels[, .(tissue, celltype = cell_type, cell_type_short)],
             by = c("tissue", "celltype"), all.x = TRUE)
  x[celltype == "JCTF_Japanese", cell_type_short := "JCTF Japanese"]
  x[trait == "Chronic_bronchitis", trait := "Chronic Bronchitis"]
  x[trait == "Chronic_Sinusitis", trait := "Chronic Sinusitis"]
  x
}

integrate <- function(contexts, smr_file, coloc_file) {
  paths <- file.path(analysis_dir, "SMR_COLOC", c(smr_file, coloc_file))
  use_snapshot <- mode == "snapshot" || (mode == "auto" && all(file.exists(paths)))
  smr <- read_method(contexts, "SMR", if (use_snapshot) paths[1])
  coloc <- read_method(contexts, "COLOC", if (use_snapshot) paths[2])
  keys <- c("tissue", "trait", "cell_type_short", "ancestry", "Gene", "probe")
  s <- smr[!is.na(cell_type_short), .(tissue, trait, cell_type_short, ancestry, Gene,
            probe = sub("[.][0-9]+$", "", probeID), p_SMR, p_SMR_bhadj, p_HEIDI)]
  pp <- coloc[!is.na(cell_type_short), .(tissue, trait, cell_type_short, ancestry, Gene,
              probe = sub("[.][0-9]+$", "", probe_name), pp3, pp4)]
  # Retain both methods' measurements whenever at least one method passes.
  hits <- unique(rbind(s[p_SMR_bhadj < .05 & p_HEIDI > .05, ..keys], pp[pp4 > .8, ..keys]))
  evidence <- merge(s[hits, on = keys], pp[hits, on = keys], by = keys)
  evidence[, `:=`(SMR_support = fifelse(p_SMR_bhadj < .05, "yes", "no"),
                 HEIDI_pass = fifelse(p_HEIDI > .05, "yes", "no"),
                 COLOC_PP4 = pp4, COLOC_PP3_PP4 = pp3 + pp4,
                 supported_by_COLOC = fifelse(pp4 > .8, "yes", "no"))]
  evidence[, supported_by_SMR_HEIDI := fcase(SMR_support == "yes" & HEIDI_pass == "yes", "yes",
                                             SMR_support == "yes" & HEIDI_pass == "no", "partial", default = "no")]
  evidence[, evidence := fcase(supported_by_SMR_HEIDI == "yes" & supported_by_COLOC == "yes", "yes",
                               default = "partial")]
  evidence[, c("pp3", "pp4", "probe") := NULL]
  setcolorder(evidence, c("Gene","trait","tissue","cell_type_short","ancestry","p_SMR","p_SMR_bhadj","p_HEIDI",
                          "SMR_support","HEIDI_pass","COLOC_PP4","COLOC_PP3_PP4",
                          "supported_by_SMR_HEIDI","supported_by_COLOC","evidence"))
  list(smr = smr[p_SMR_bhadj < .05 & p_HEIDI > .05], coloc = coloc[pp4 > .8], evidence = evidence)
}

write_tables <- function(result, blood_only) {
  x <- copy(result$evidence[cell_type_short != "JCTF Japanese"])
  prefix <- if (blood_only) "autoimmune_Blood" else "respiratory_Blood_Lung"
  if (blood_only) x[, tissue := NULL]
  setorderv(x, if (blood_only) c("ancestry","trait","Gene","cell_type_short") else
              c("trait","Gene","tissue","cell_type_short"))
  fwrite(x, file.path(out, paste0(prefix, if (blood_only) "_evidence.tsv" else "_evidence_R2.tsv")),
         sep = "\t", na = "NA")
  for (method in c("smr", "coloc")) {
    y <- copy(result[[method]][!is.na(cell_type_short) & celltype != "JCTF_Japanese"])
    y[, celltype := cell_type_short][, cell_type_short := NULL]
    if (blood_only && method == "smr") y[, tissue := NULL]
    # Keep the native method columns first, followed by the analysis context.
    context_cols <- if (method == "smr") {
      c("tissue", "ancestry", "celltype", "trait", "p_SMR_bhadj", "p_SMR_bfradj")
    } else if (blood_only) c("trait", "celltype", "ancestry", "tissue") else c("tissue", "trait", "celltype", "ancestry")
    setcolorder(y, c(setdiff(names(y), context_cols), intersect(context_cols, names(y))))
    suffix <- if (blood_only) "_results.tsv" else "_results_EUR_R2.tsv"
    fwrite(y, file.path(out, paste0(prefix, "_", method, suffix)), sep = "\t", na = "NA")
    if (!blood_only) for (ts in c("Blood", "Lung")) {
      fwrite(y[tissue == ts], file.path(out, paste0("respiratory_", ts, "_", method, "_results_EUR.tsv")),
             sep = "\t", na = "NA")
    }
  }
}

blood <- integrate(CJ(tissue = "Blood", trait = blood_traits, celltype = c(cell_order$Blood, "JCTF_Japanese"),
                      ancestry = c("EAS_EAS", "EUR_EAS")), "Blood_smr_results.rds", "Blood_coloc_results.rds")
resp <- integrate(rbindlist(lapply(names(cell_order), function(ts) {
  CJ(tissue = ts, trait = resp_traits, celltype = cell_order[[ts]], ancestry = "EUR_EUR")
})), "Lung_Blood_smr_results_EUR_EUR_R2.rds", "Blood_Lung_coloc_results_EUR_EUR_R2.rds")
write_tables(blood, TRUE)
write_tables(resp, FALSE)

cluster_genes <- function(x, traits) {
  presence <- unique(x[, .(Gene, trait = factor(trait, levels = traits))])[order(trait, Gene)]
  genes <- unique(presence$Gene)
  if (length(genes) < 2) return(genes)
  presence[, Gene := factor(Gene, levels = genes)]
  tree <- hclust(dist(xtabs(~ Gene + trait, presence), method = "binary"))
  tree$labels[tree$order]
}

pie_panel <- function(x, cells, genes, traits, xlabel, blood_comparator = FALSE) {
  counts <- unique(x[, .(cell_type_short, Gene, trait)])[, present := 1L]
  pies <- dcast(counts, cell_type_short + Gene ~ trait, value.var = "present", fill = 0L)
  shown <- intersect(traits, names(pies))
  xpos <- seq_along(cells) + ifelse(cells == "JCTF Japanese", 2L, 0L)
  pies[, `:=`(x = xpos[match(cell_type_short, cells)], y = match(Gene, genes))]
  grid <- CJ(x = xpos, y = seq_along(genes))
  p <- ggplot() +
    geom_tile(data = grid, aes(x, y), width = if (blood_comparator) 0.95 else 1,
              height = if (blood_comparator) 0.95 else 1, fill = NA, color = "lightgrey") +
    geom_scatterpie(data = as.data.frame(pies), aes(x, y), cols = shown, color = NA,
                    pie_scale = if (blood_comparator) 1 else 1.4) +
    scale_x_continuous(breaks = xpos, labels = cells, expand = c(.03,.03)) +
    scale_y_continuous(breaks = seq_along(genes), labels = genes, expand = c(.01,.01)) +
    coord_fixed() + labs(x = xlabel, y = "Gene", fill = "Trait") + theme_classic() +
    theme(text = element_text(family = "Helvetica"), axis.text = element_text(size = 10, color = "black"),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
          axis.text.y = element_text(face = "italic"), axis.title = element_text(size = if (blood_comparator) 14 else 12))
  if (blood_comparator) p + scale_fill_brewer(palette = "Set2") else
    p + scale_fill_manual(values = c(IPF = "dodgerblue", COPD = "firebrick"))
}

save_plot <- function(p, filename, width, height) {
  ggsave(file.path(out, filename), p, width = width, height = height, units = "in", dpi = 300, bg = "white")
}

dual <- blood$evidence[ancestry == "EAS_EAS" & evidence == "yes"]
retain <- dual[cell_type_short != "JCTF Japanese" & !grepl("^RPS|[.]", Gene), unique(Gene)]
dual <- dual[Gene %in% retain]
cells <- intersect(c(display_order("Blood"), "JCTF Japanese"), dual$cell_type_short)
save_plot(pie_panel(dual, cells, cluster_genes(dual, blood_traits), blood_traits, "", TRUE),
          "SMR_COLOC_Blood_EAS_celltype_by_gene.png", 6, 15)

blood_count_plot <- function(hits, title) {
  traits <- c("CD","IBD","RA","SLE","Asthma","COPD","AR")
  genes <- unique(hits[celltype != "JCTF_Japanese" & trait %in% traits, .(trait, Gene, ancestry)])
  eas <- genes[ancestry == "EAS_EAS", .(trait, Gene)]
  eur <- genes[ancestry == "EUR_EAS", .(trait, Gene)]
  eas[, shared := FALSE][eur, on = .(trait, Gene), shared := TRUE]
  counts <- eas[, .N, by = .(trait, shared)]
  counts[, group := factor(shared, levels = c(FALSE, TRUE), labels = c("EAS-specific","EUR & EAS prioritized"))]
  ggplot(counts, aes(factor(trait, levels = traits), N, fill = group)) + geom_col(width = .6) +
    scale_x_discrete(drop = FALSE) + scale_fill_manual(values = c("#ffa23f","#6384e9"), drop = FALSE) +
    labs(x = "Trait", y = "Number of prioritized genes", fill = NULL, title = title) + theme_classic() +
    theme(text = element_text(size = 12, family = "Helvetica"), axis.text = element_text(color = "black", size = 12),
          axis.title = element_text(size = 14, face = "bold"), plot.title = element_text(hjust = .5, face = "bold"),
          legend.text = element_text(size = 12), plot.margin = margin(b = 10),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
}
p <- (blood_count_plot(blood$coloc, "COLOC") / blood_count_plot(blood$smr, "SMR")) +
  plot_layout(axis_titles = "collect", guides = "collect") & theme(legend.position = "bottom")
save_plot(p, "Blood_SMR_COLOC_EAS_ngenes_new.png", 5, 7)

dual <- resp$evidence[evidence == "yes" & trait %in% c("IPF","COPD")]
genes <- cluster_genes(dual, c("IPF","COPD"))
panels <- lapply(c("Lung","Blood"), function(ts) {
  x <- dual[tissue == ts]
  pie_panel(x, intersect(display_order(ts), x$cell_type_short), genes, c("IPF","COPD"), paste0(ts, "-EUR"))
})
panels[[2]] <- panels[[2]] + labs(y = NULL) + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
p <- wrap_plots(panels, nrow = 1, guides = "collect") & theme(legend.position = "top")
save_plot(p, "SMR_COLOC_Lung_Blood_EUR_celltype_by_gene_COPD_IPF.png", 5, 6)

counts <- rbindlist(lapply(resp[c("smr","coloc")], function(x) {
  x[tissue == "Lung" & !is.na(cell_type_short) & !trait %in% c("SBP","DBP"),
    .(n_genes = uniqueN(Gene)), by = cell_type_short]
}), idcol = "method")
counts[, method := toupper(method)]
counts[, cell_type_short := factor(cell_type_short, levels = rev(display_order("Lung")))]
lung_colors <- labels[tissue == "Lung", setNames(color, cell_type_short)]
p <- ggplot(counts, aes(n_genes, cell_type_short, color = cell_type_short)) + geom_point(size = 2.5) +
  scale_color_manual(values = lung_colors) + facet_wrap(~method, axes = "all", axis.labels = "margins", scales = "free_x") +
  labs(x = "Number of prioritized genes", y = "Cell type") + theme_classic() +
  theme(text = element_text(size = 14, family = "Helvetica"), axis.text = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 14), strip.text = element_text(size = 14), axis.ticks.y = element_blank(),
        strip.background = element_blank(), legend.position = "none", panel.grid.major.y = element_line(linetype = "dashed", color = "grey"))
save_plot(p, "SMR_COLOC_Lung_EUR_celltype_neGenes.png", 6, 7)
