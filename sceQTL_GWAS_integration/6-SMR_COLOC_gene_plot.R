#!/usr/bin/env Rscript
# Source this file for effect_plot(), locus_plot(), and expression_plot(),
# or run Rscript 5-SMR_COLOC_gene_plots_public.R --help for file-based usage.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

# SMR --plot output supplies harmonized GWAS/eQTL effects and signed LD.
effect_plot <- function(plot_file, title = NULL, r2_min = 0.1, cis_kb = 2000) {
  stopifnot(r2_min >= 0, r2_min < 1, cis_kb > 0)
  lines <- readLines(plot_file)
  fields <- function(line) strsplit(trimws(line), "[[:space:]]+")[[1]]
  section <- function(name) {
    i <- grep(paste0("^\\$", name, "[[:space:]]"), lines)
    if (length(i) != 1) stop("Missing or ambiguous $", name, " section in ", plot_file)
    i
  }
  records <- function(i, n, columns) {
    if (!is.finite(n) || n < 1 || i + n > length(lines)) stop("Invalid record count in ", plot_file)
    x <- fread(text = paste0(paste(lines[i + seq_len(n)], collapse = "\n"), "\n"), header = FALSE)
    if (ncol(x) < length(columns)) stop("Incomplete record in ", plot_file)
    setnames(x, seq_along(columns), columns)
    x[, ..columns]
  }
  read_section <- function(name, columns) {
    i <- section(name)
    records(i, as.integer(fields(lines[i])[2]), columns)
  }
  target <- fields(lines[section("probe")])[3]
  probes <- read_section("probe", c("probe", "chr", "bp", "gene"))
  probe <- probes[probe == target]
  if (nrow(probe) != 1L) stop("Expected one target probe: ", target)
  snps <- read_section("SNP", c("SNP", "chr", "bp"))
  gwas <- read_section("GWAS", c("SNP", "b_gwas", "se_gwas"))
  i <- section("eQTL")
  n_probes <- as.integer(fields(lines[i])[2])
  eqtl <- NULL
  for (j in seq_len(n_probes)) {
    i <- i + 1L
    header <- fields(lines[i])
    n <- as.integer(header[2])
    if (header[1] == target) {
      eqtl <- records(i, n, c("SNP", "b_eqtl", "se_eqtl", "ld"))
      break
    }
    i <- i + n
  }
  if (is.null(eqtl)) stop("No eQTL section for ", target)
  x <- Reduce(function(a, b) merge(a, b, by = "SNP"), list(eqtl, gwas, snps))
  x <- x[chr == probe$chr & abs(bp - probe$bp) < 1000 * cis_kb & is.finite(ld) &
         is.finite(b_eqtl) & se_eqtl > 0 & is.finite(se_eqtl) & is.finite(b_gwas) &
         is.finite(se_gwas) & se_gwas > 0]
  if (!nrow(x)) stop("No matched cis variants with effect sizes and LD for ", target)
  lead <- x[which.max(abs(b_eqtl / se_eqtl))]
  x <- x[ld^2 > r2_min]
  if (!nrow(x) || lead$b_eqtl == 0) stop("No usable effect estimates at the requested LD threshold.")
  x[, type := factor(SNP == lead$SNP, levels = c(TRUE, FALSE), labels = c("top cis-eQTL","cis-eQTL"))]
  if (is.null(title)) title <- probe$gene
  ggplot(x, aes(b_eqtl, b_gwas)) +
    geom_segment(aes(x = b_eqtl - se_eqtl, xend = b_eqtl + se_eqtl, yend = b_gwas), alpha = .5) +
    geom_segment(aes(y = b_gwas - se_gwas, yend = b_gwas + se_gwas, xend = b_eqtl), alpha = .5) +
    geom_abline(intercept = 0, slope = lead$b_gwas / lead$b_eqtl, linetype = "dashed", color = "orange") +
    geom_point(aes(color = type, shape = type, size = ld^2), alpha = .8) +
    scale_color_manual(values = c("top cis-eQTL" = "red", "cis-eQTL" = "black"), name = "SNP") +
    scale_shape_manual(values = c("top cis-eQTL" = 17, "cis-eQTL" = 16), name = "SNP") +
    scale_size_continuous(name = expression(italic(r)^2), range = c(1, 5)) +
    labs(x = "eQTL effect sizes", y = "GWAS effect sizes", title = title) + theme_classic(base_family = "Helvetica")
}

# One manifest row per GWAS/eQTL track, in display order. Paths are relative to
# the manifest. Supply summary statistics on the same genome build as the gene track.
locus_plot <- function(manifest, chr, start, end, genes_file = NULL) {
  tracks <- fread(manifest, na.strings = c("", "NA"))
  stopifnot(all(c("track", "file") %in% names(tracks)), nrow(tracks) > 0,
            !anyDuplicated(tracks$track), is.finite(start), is.finite(end), start < end)
  resolve <- function(path) {
    if (grepl("^(/|[A-Za-z]:[/\\\\])", path)) path else file.path(dirname(manifest), path)
  }
  chrom <- function(x) sub("^chr", "", as.character(x), ignore.case = TRUE)
  region_chr <- chrom(chr)
  left_bp <- start
  right_bp <- end
  panels <- lapply(seq_len(nrow(tracks)), function(i) {
    x <- fread(resolve(tracks$file[i]))
    if (!"chr" %in% names(x)) setnames(x, "Chr", "chr")
    if (!"pos" %in% names(x)) setnames(x, "BP", "pos")
    x <- x[chrom(chr) == region_chr & pos >= left_bp & pos <= right_bp & is.finite(p) & p >= 0 & p <= 1]
    if ("probe" %in% names(tracks) && !is.na(tracks$probe[i])) {
      if (!"Probe" %in% names(x)) stop("Missing Probe column in ", tracks$file[i])
      x <- x[Probe == tracks$probe[i]]
    }
    if (!nrow(x) || anyDuplicated(x$SNP)) stop("Expected one row per SNP in region for ", tracks$track[i])
    lead <- if ("lead" %in% names(tracks) && !is.na(tracks$lead[i])) tracks$lead[i] else x$SNP[which.min(x$p)]
    if (!lead %in% x$SNP) stop("Lead SNP is absent from track ", tracks$track[i], ": ", lead)
    if ("ld_file" %in% names(tracks) && !is.na(tracks$ld_file[i])) {
      ld <- fread(resolve(tracks$ld_file[i]), select = c("SNP_A", "SNP_B", "R2"))
      ld <- ld[SNP_A == lead]
      if (!nrow(ld) || anyDuplicated(ld$SNP_B)) stop("Expected LD relative to ", lead, " for ", tracks$track[i])
      x[, r2 := ld$R2[match(SNP, ld$SNP_B)]]
    }
    if (!"r2" %in% names(x)) x[, r2 := NA_real_]
    if (any(x$r2 < 0 | x$r2 > 1, na.rm = TRUE)) stop("Expected squared LD (0 to 1) for ", tracks$track[i])
    x[SNP == lead, r2 := 1]
    x[, logp := -log10(pmax(p, .Machine$double.xmin))]
    ggplot(x, aes(pos / 1e6, logp)) +
      geom_hline(yintercept = -log10(5e-8), linetype = "dashed", color = "grey75") +
      geom_point(aes(color = r2), size = 1.5) +
      geom_point(data = x[SNP == lead], shape = 23, fill = "purple", size = 2) +
      geom_text(data = x[SNP == lead], aes(label = SNP), vjust = -0.7, size = 3) +
      scale_color_gradientn(colors = c("navy","skyblue","gold","red"), limits = c(0,1),
                            na.value = "grey70", name = expression(italic(r)^2)) +
      scale_y_continuous(expand = expansion(mult = c(.03,.18))) +
      coord_cartesian(xlim = c(start, end) / 1e6) +
      labs(title = tracks$track[i], x = paste0("Chromosome ", region_chr, " (Mb)"), y = expression(-log[10](P))) +
      theme_classic(base_family = "Helvetica")
  })
  heights <- rep(1, length(panels))
  if (!is.null(genes_file)) {
    g <- fread(genes_file)
    g <- g[chrom(chr) == region_chr & end >= left_bp & start <= right_bp]
    if (nrow(g)) {
      setorder(g, start)
      g[, lane := seq_len(.N)]
      g[, mid := (pmax(start, left_bp) + pmin(end, right_bp)) / 2e6]
      panels[[length(panels) + 1L]] <- ggplot(g, aes(y = lane)) +
        geom_segment(aes(x = start / 1e6, xend = end / 1e6, yend = lane), color = "grey40") +
        geom_text(aes(x = mid, label = gene_name), vjust = -0.35, fontface = "italic", size = 3) +
        scale_y_continuous(expand = expansion(add = 1)) + coord_cartesian(xlim = c(start,end) / 1e6) +
        labs(x = paste0("Chromosome ", chrom(chr), " (Mb)"), y = NULL) + theme_classic(base_family = "Helvetica") +
        theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank())
      heights <- c(heights, max(.55, .13 * nrow(g)))
    }
  }
  wrap_plots(panels, ncol = 1, heights = heights, guides = "collect", axis_titles = "collect_x")
}

# A tensorQTL expression BED, donor covariates, PLINK additive .raw file, and
# SMR --query table for the same cell type and ancestry. A1 is the eQTL effect allele.
expression_plot <- function(expr_bed, cov_file, genotype_raw, eqtl_file, probe, snp,
                            title = NULL, covariates = NULL, dosage_column = NULL) {
  eqtl <- fread(eqtl_file)
  v <- eqtl[SNP == snp & Probe == probe]
  if (nrow(v) != 1L) stop("Expected one eQTL row for ", probe, " and ", snp)
  ids <- fread(expr_bed, select = "TargetID")[[1]]
  row <- which(ids == probe)
  if (length(row) != 1L) stop("Expected one BED row for ", probe)
  samples <- names(fread(expr_bed, nrows = 0))[-(1:4)]
  values <- fread(expr_bed, skip = row, nrows = 1, header = FALSE)
  expr <- data.table(IID = samples, expression = as.numeric(unlist(values[, -(1:4)], use.names = FALSE)))
  covs <- fread(cov_file, colClasses = list(character = "sample"))
  if (anyDuplicated(covs$sample) || anyDuplicated(samples)) stop("Duplicate expression or covariate sample IDs.")
  if (is.null(covariates)) {
    covariates <- intersect(c(grep("PF|Condition_group", names(covs), value = TRUE),
                              "Gender_pred","Age_pred","Smoking_status_pred",paste0("PC",1:5)), names(covs))
  }
  selected <- c("sample", covariates)
  x <- na.omit(merge(expr, covs[, ..selected], by.x = "IID", by.y = "sample"))
  if (!nrow(x)) stop("No complete expression and covariate records.")
  varying <- covariates[vapply(x[, ..covariates], uniqueN, integer(1)) > 1L]
  model_data <- x[, c("expression", varying), with = FALSE]
  fit <- if (length(varying)) lm(expression ~ ., data = model_data) else lm(expression ~ 1, data = model_data)
  x[, residual := residuals(fit)]

  columns <- names(fread(genotype_raw, nrows = 0))
  if (is.null(dosage_column)) {
    ids <- c(snp, paste(v$Chr, v$BP, v$A2, v$A1, "b38", sep = "_"),
             paste(v$Chr, v$BP, v$A1, v$A2, "b38", sep = "_"))
    dosage_column <- intersect(columns, c(paste0(ids, "_", v$A1), paste0(ids, "_", v$A2)))
  }
  if (length(dosage_column) != 1L || !dosage_column %in% columns) stop("Specify the PLINK dosage column for ", snp)
  counted <- sub(".*_", "", dosage_column)
  if (!counted %in% c(v$A1, v$A2)) stop("PLINK counted allele does not match the eQTL alleles.")
  geno <- fread(genotype_raw, select = c("IID", dosage_column), colClasses = list(character = "IID"))
  if (anyDuplicated(geno$IID)) stop("Duplicate genotype sample IDs.")
  setnames(geno, dosage_column, "dosage")
  if (counted == v$A2) geno[, dosage := 2 - dosage]
  x <- merge(x, geno, by = "IID")[!is.na(dosage)]
  if (!nrow(x) || any(!x$dosage %in% 0:2)) stop("Expected matched diploid hard calls (0, 1, 2).")
  genotype <- c(paste0(v$A2,v$A2), paste0(v$A2,v$A1), paste0(v$A1,v$A1))
  if (is.null(title)) title <- paste(v$Gene, snp)
  ggplot(x, aes(dosage, residual)) +
    geom_violin(aes(group = dosage), color = "skyblue", linewidth = .6) +
    geom_boxplot(aes(group = dosage), color = "skyblue", width = .2, outlier.size = .5) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "#be0000", linewidth = .6) +
    scale_x_continuous(breaks = 0:2, labels = genotype, limits = c(-.5, 2.5)) +
    labs(title = title, subtitle = sprintf("eQTL P = %.3g", v$p), x = snp, y = "Covariate-adjusted expression") +
    theme_classic(base_family = "Helvetica") + theme(plot.title = element_text(hjust = .5, face = "bold"),
                                                    plot.subtitle = element_text(hjust = .5))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  usage <- paste(
    "Rscript 5-SMR_COLOC_gene_plots_public.R effect SMR_PLOT.txt OUT.png [TITLE]",
    "Rscript 5-SMR_COLOC_gene_plots_public.R locus TRACKS.tsv CHR START_BP END_BP OUT.png [GENES.tsv]",
    "Rscript 5-SMR_COLOC_gene_plots_public.R expression EXPR.bed COV.tsv GENOTYPE.raw EQTL.tsv PROBE SNP OUT.png [TITLE]",
    "TRACKS.tsv: track, file; optional ld_file (PLINK .ld), lead (LD index SNP), probe.",
    "Regional files: SNP, chr, pos, p (or SMR query Chr, BP); optional r2 relative to lead.",
    "Optional gene annotation: chr, start, end, gene_name.",
    "Expression: tensorQTL BED with TargetID; covariates keyed by sample; PLINK .raw keyed by IID;",
    "SMR --query output with SNP, Chr, BP, A1, A2, Probe, Gene, p. Dosages must be diploid hard calls.",
    "Source expression_plot() to specify covariates or a nonstandard PLINK dosage column.", sep = "\n")
  if (!length(args) || args[1] %in% c("-h", "--help")) {
    cat(usage, "\n")
  } else {
    mode <- args[1]
    a <- args[-1]
    if (mode == "effect" && length(a) %in% 2:3) {
      p <- effect_plot(a[1], if (length(a) == 3) a[3])
      output <- a[2]; width <- 5; height <- 3.2
    } else if (mode == "locus" && length(a) %in% 5:6) {
      p <- locus_plot(a[1], a[2], as.numeric(a[3]), as.numeric(a[4]), if (length(a) == 6) a[6])
      output <- a[5]; width <- 7; height <- 2 * length(p)
    } else if (mode == "expression" && length(a) %in% 7:8) {
      p <- expression_plot(a[1], a[2], a[3], a[4], a[5], a[6], if (length(a) == 8) a[8])
      output <- a[7]; width <- 4; height <- 3.2
    } else stop(usage, call. = FALSE)
    ggsave(output, p, width = width, height = height, units = "in", dpi = 300, bg = "white")
  }
}
