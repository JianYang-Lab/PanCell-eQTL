#!/usr/bin/env Rscript
# EUR single-cell eGenes versus eQTLGen (Blood) and GTEx v8 (other tissues).
# Run from the project root, for example:
#   Rscript egene_bulk_status_comparison.R
#   Rscript egene_bulk_status_comparison.R --workers=2 Blood Lung
#
# Supply the full eQTLGen association file in references/eQTLGen/ and the GTEx
# v8 allpairs files in references/GTEx_v8/. To use other locations, set
# EQTLGEN_FILE to the eQTLGen file and GTEX_DIR to the directory of GTEx files.
# Single-cell inputs: <tissue>/sc-eQTL_R/results/raw/<cell_type>_EUR_pc5_perm.tsv
# Expression inputs: crossTissue_analysis_R/expression_decile/<tissue>_EUR_cross_celltype.txt
# Gene annotation: <tissue>/resource/genes.gtf (tabular gene_id and gene_name).
#
# A tissue eGene has nominal p < 5e-8 in at least one EUR cell type or, for bulk,
# at least one cis association below the same threshold. Single-cell tests use
# autosomes and lead-variant AF in [0.01, 0.99]. GTEx comparisons exclude chrX.
# One-sided Wilcoxon and two-sided Fisher tests compare additional and shared
# eGenes. BH adjustment is separate for each test family across tissues,
# three contrasts and two features.
# Requires data.table, ggplot2, cowplot, ragg, bash, gzip and awk.

suppressPackageStartupMessages(library(data.table))

eqtlgen_filename <- "2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz"
bulk_names <- c(Blood = "eQTLGen", Lung = "Lung",
                Skin = "Skin_Sun_Exposed_Lower_leg", Liver = "Liver",
                Colon = "Colon_Sigmoid")
celltypes_all_list <- list(
    Blood = c("CD4TNC", "CD4TNCInf", "CD4TEM", "Treg", "CD4TEMRA", "CD8TNC", "CD8TRM", "CD8TEMRA", "MAIT", "NKp", "NKn", "gdT", "BIN", "BMemNS", "BMemS", "Plasma", "Plasmablasts", "MonoC", "MonoInf", "MonoNC", "DC1", "DC2", "Macro", "Nph", "pDC", "Eryth", "HSPC", "Platelet"),
    Lung = c("CD4TNC", "CD4TEM", "Treg", "CD8TEMRA", "CD8TEMTRM", "CD8TRM", "Tfh", "NKp", "NKn", "NKCyc", "BIN", "BMem", "Plasma", "MonoC", "MonoNC", "Macro", "MacroCyc", "MacroInt", "MacroAlv", "Mast", "DC1", "DC2", "MigDC", "pDC", "Eryth", "Basal", "Suprabasal", "Ciliated", "Club", "Goblet", "Ion", "AT1", "AT2", "AeroCap", "ECArt", "ECCap", "ECVenPulm", "ECVenSys", "LECdiff", "LECmat", "FibroAdv", "FibroAlv", "Myofibro", "Pericyte", "Meso"),
    Skin = c("CD4TNC", "CD4TEM", "Treg", "CD8TEMRA", "CD8TEMTRM", "ILC3", "NKp", "NKn", "BMem", "Plasma", "MonoMac", "Macro1", "Macro2", "MacroInf", "Mast", "DC1", "DC2", "moDC", "pDC", "LC", "MigLC", "MigDC", "KCundiff", "KCdiff", "VEC1", "VEC2", "VEC3", "LE1", "LE2", "Fibro1", "Fibro2", "Fibro3", "Pericyte1", "Pericyte2", "Melano", "Schwann1", "Schwann2"),
    Liver = c("CD4TNC", "CD4TEM", "Treg", "CD8TEMRA", "CD8TRM", "MAIT", "NKp", "NKn", "BIN", "BMem", "Plasma", "MonoC", "MonoNC", "Macro", "KC", "Mast", "Nph", "DC1", "DC2", "HepPP", "HepIZ", "HepPC", "CGC", "Endo_portal", "LSEC_PP", "LSEC_CV", "Fibro", "HSC"),
    Colon = c("CD4TNC", "CD4TEM", "Tfh", "Th1", "Th17", "Treg", "CD8TEMTRM", "CD8TRM", "gdT", "ILC", "NKp", "NKn", "Bf", "BGC", "BIN", "BMem", "Plasma", "Mono", "Macro", "MacroIntest", "Mast", "DC2", "Stem", "TA", "TA_CycSec", "Entero", "Entero_BEST4", "Colono", "GobletImm", "Goblet", "Paneth", "Tuft", "EEC", "ECArt", "ECCap", "ECVen", "LEC", "Stromal1", "Stromal2", "Stromal3", "Myofibro", "Pericyte", "Glia")
)

normalize_gene <- function(x) sub("\\.[0-9]+$", "", x)
safe_fraction <- function(n, d) if (d == 0L) NA_real_ else n / d

load_sc <- function(tissue, project_dir, cutoff) {
    tables <- manifest <- list()
    for (ct in celltypes_all_list[[tissue]]) {
        path <- file.path(project_dir, tissue, "sc-eQTL_R", "results", "raw",
                          paste0(ct, "_EUR_pc5_perm.tsv"))
        entry <- data.table(tissue = tissue, ancestry = "EUR", cell_type = ct,
                            path = path, available = file.exists(path),
                            n_sc_tested_genes = NA_integer_, n_sc_egenes = NA_integer_)
        if (entry$available) {
            x <- fread(path, select = c("phenotype_id", "variant_id", "pval_nominal", "af"))
            x <- x[sub("^chr", "", sub("_.*", "", variant_id)) %in% as.character(1:22)]
            if (any(!is.na(x$pval_nominal) & (!is.finite(x$pval_nominal) |
                                           x$pval_nominal < 0 | x$pval_nominal > 1))) {
                stop("Invalid sc p-value in ", path)
            }
            x[, phenotype_id := normalize_gene(phenotype_id)]
            if (anyNA(x$phenotype_id) || any(x$phenotype_id == "")) stop("Missing gene ID: ", path)
            x[, eligible_p := fifelse(is.finite(pval_nominal) & !is.na(af) &
                                      af >= 0.01 & af <= 0.99, pval_nominal, Inf)]
            x <- x[, .(sc_min_pval = min(eligible_p)), by = phenotype_id]
            x[, `:=`(cell_type = ct, sc_is_egene = sc_min_pval < cutoff)]
            entry[, `:=`(n_sc_tested_genes = nrow(x), n_sc_egenes = sum(x$sc_is_egene))]
            tables[[ct]] <- x
        }
        manifest[[ct]] <- entry
    }
    manifest <- rbindlist(manifest)
    missing <- manifest[available == FALSE, cell_type]
    if (length(missing)) warning(tissue, ": missing EUR pc5 files: ", paste(missing, collapse = ", "))
    if (!length(tables)) stop("No sc input files for ", tissue)
    list(data = rbindlist(tables), manifest = manifest,
         cell_types = manifest[available == TRUE, cell_type])
}

load_bulk_gene_summary <- function(path, source, cache_dir) {
    info <- file.info(path)
    if (is.na(info$size)) stop("Missing bulk input: ", path)
    signature <- list(version = 1L, path = normalizePath(path), size = info$size,
                      mtime = as.numeric(info$mtime), source = source)
    cache <- file.path(cache_dir, paste0(source, "_gene_min_p.rds"))
    if (file.exists(cache)) {
        cached <- readRDS(cache)
        if (identical(cached$signature, signature)) {
            message("Using bulk cache: ", source)
            return(copy(cached$data))
        }
    }
    message("Streaming full bulk associations: ", source)
    awk_file <- tempfile(fileext = ".awk")
    result_file <- tempfile(fileext = ".tsv")
    on.exit(unlink(c(awk_file, result_file)), add = TRUE)
    # Track all genes before applying a significance cutoff. Missing p-values
    # must never be coerced to zero (which would make them significant).
    writeLines(r"(
BEGIN { FS = OFS = "\t" }
NR == 1 {
    sub(/\r$/, "", $NF)
    for (i = 1; i <= NF; i++) col[$i] = i
    g = col[gene_name]; p = col[p_name]; v = col["variant_id"]
    if (!g || !p || (exclude_x && !v)) { failed = 1; exit 2 }
    next
}
{
    if (exclude_x && $v ~ /^(chr)?X_/) next
    gene = $g; sub(/[.][0-9]+$/, "", gene)
    if (gene == "" || gene == "NA") { failed = 1; exit 3 }
    count[gene]++
    value = $p; sub(/\r$/, "", value)
    if (value == "" || value == "NA" || value == "NaN" || value == ".") next
    if (value !~ /^[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][+-]?[0-9]+)?$/ ||
        value + 0 < 0 || value + 0 > 1) { failed = 1; exit 4 }
    value += 0
    valid[gene]++
    if (!(gene in best) || value < best[gene]) best[gene] = value
}
END {
    if (failed) exit 5
    print "phenotype_id", "bulk_n_pairs", "bulk_n_valid_pvalues", "bulk_min_pval"
    for (gene in count) {
        if (valid[gene]) printf "%s\t%.0f\t%.0f\t%.17g\n", gene, count[gene], valid[gene], best[gene]
        else printf "%s\t%.0f\t0\tNA\n", gene, count[gene]
    }
}
)", awk_file)
    is_eqtlgen <- source == "eQTLGen"
    command <- paste(
        "gzip -dc --", shQuote(path), "| LC_ALL=C awk",
        "-v", shQuote(paste0("gene_name=", if (is_eqtlgen) "Gene" else "gene_id")),
        "-v", shQuote(paste0("p_name=", if (is_eqtlgen) "Pvalue" else "pval_nominal")),
        "-v", shQuote(paste0("exclude_x=", as.integer(!is_eqtlgen))),
        "-f", shQuote(awk_file)
    )
    status <- system2("bash", c("-o", "pipefail", "-c", shQuote(command)), stdout = result_file)
    if (status != 0L) stop("Bulk scan failed (exit ", status, "): ", path)
    x <- fread(result_file, colClasses = c(phenotype_id = "character"))
    if (!nrow(x)) stop("No bulk genes found: ", path)
    if (any(x$bulk_n_valid_pvalues == 0L)) {
        stop("Bulk genes with no valid p-value cannot be classified: ", path)
    }
    setorder(x, phenotype_id)
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    temporary_cache <- tempfile(tmpdir = cache_dir, fileext = ".rds")
    saveRDS(list(signature = signature, data = x), temporary_cache)
    if (!file.rename(temporary_cache, cache)) stop("Could not save cache: ", cache)
    x
}

summarize_status <- function(x) {
    total <- nrow(x)
    identified <- sum(x$bulk_status == "significant")
    not_tested <- sum(x$bulk_status == "not_tested")
    nonsignificant <- sum(x$bulk_status == "tested_not_significant")
    not_identified <- not_tested + nonsignificant
    stopifnot(total == identified + not_identified)
    data.table(n_sc_egenes = total, n_bulk_identified = identified,
               n_bulk_not_identified = not_identified, n_bulk_not_tested = not_tested,
               n_bulk_tested_not_significant = nonsignificant,
               frac_bulk_not_identified = safe_fraction(not_identified, total),
               frac_bulk_not_tested = safe_fraction(not_tested, total),
               frac_bulk_tested_not_significant = safe_fraction(nonsignificant, total),
               frac_not_identified_not_tested = safe_fraction(not_tested, not_identified),
               frac_not_identified_tested_not_significant = safe_fraction(nonsignificant, not_identified))
}

add_expression <- function(genes, tissue, cell_types, project_dir) {
    genes[, `:=`(gene_name = NA_character_, expr_max = NA_real_,
                  expr_ct_count = NA_integer_, expr_decile = NA_integer_)]
    annotation_path <- file.path(project_dir, tissue, "resource", "genes.gtf")
    if (!file.exists(annotation_path)) annotation_path <- file.path(project_dir, "Blood", "resource", "genes.gtf")
    expression_path <- file.path(project_dir, "crossTissue_analysis_R", "expression_decile",
                                 paste0(tissue, "_EUR_cross_celltype.txt"))
    if (!file.exists(annotation_path) || !file.exists(expression_path)) {
        warning("Expression annotation unavailable for ", tissue, "; retaining all eGenes in counts.")
        return(genes)
    }
    annotation <- fread(annotation_path, select = c("gene_id", "gene_name"))
    annotation[, gene_id := normalize_gene(gene_id)]
    annotation <- unique(annotation, by = "gene_id")
    genes[, gene_name := annotation$gene_name[match(phenotype_id, annotation$gene_id)]]
    expr <- fread(expression_path)
    columns <- intersect(paste0("expr_mean_", cell_types, "_EUR"), names(expr))
    if (!length(columns)) stop("No matching expression columns: ", expression_path)
    if (anyDuplicated(expr$gene_name)) stop("Duplicate expression gene names: ", expression_path)
    values <- as.matrix(expr[, ..columns])
    expr[, expr_ct_count := rowSums(is.finite(values))]
    values[!is.finite(values)] <- NA_real_
    expr[, expr_max := do.call(pmax, c(as.data.frame(values), list(na.rm = TRUE)))]
    expr[expr_ct_count == 0L, expr_max := NA_real_]
    ix <- match(genes$gene_name, expr$gene_name)
    genes[, `:=`(expr_max = expr$expr_max[ix], expr_ct_count = expr$expr_ct_count[ix])]
    # Quantile bins within each tissue's sc eGenes; repeated boundaries are safe.
    good <- which(is.finite(genes$expr_max))
    if (length(good)) {
        boundaries <- quantile(genes$expr_max[good], seq(0, 1, 0.1), names = FALSE)
        genes[good, expr_decile := findInterval(expr_max, boundaries,
                                               left.open = TRUE, all.inside = TRUE)]
    }
    genes
}

compare_tissue <- function(tissue, project_dir, gtex_dir, eqtlgen_file, cache_dir, cutoff) {
    message("Comparing ", tissue, " EUR against ", bulk_names[[tissue]])
    sc <- load_sc(tissue, project_dir, cutoff)
    bulk_file <- if (tissue == "Blood") eqtlgen_file else file.path(gtex_dir, paste0(
        "GTEx_Analysis_v8_QTLs_GTEx_Analysis_v8_eQTL_all_associations_",
        bulk_names[[tissue]], ".allpairs.txt.gz"))
    bulk <- load_bulk_gene_summary(bulk_file, bulk_names[[tissue]], cache_dir)
    genes <- sc$data[, .(sc_min_pval = min(sc_min_pval), ct_count = .N,
                         sc_egene_ct_count = sum(sc_is_egene),
                         sc_egene_cell_types = paste(cell_type[sc_is_egene], collapse = ";")),
                     by = phenotype_id][sc_egene_ct_count > 0L]
    genes <- merge(genes, bulk, by = "phenotype_id", all.x = TRUE, sort = FALSE)
    genes[, bulk_tested := !is.na(bulk_n_pairs)]
    genes[, bulk_significant := bulk_tested & bulk_min_pval < cutoff]
    genes[, bulk_status := fifelse(!bulk_tested, "not_tested",
                                   fifelse(bulk_significant, "significant", "tested_not_significant"))]
    genes[, `:=`(tissue = tissue, ancestry = "EUR", bulk_source = bulk_names[[tissue]],
                  pval_cutoff = cutoff)]
    genes <- add_expression(genes, tissue, sc$cell_types, project_dir)
    summary <- summarize_status(genes)
    summary[, `:=`(tissue = tissue, ancestry = "EUR", bulk_source = bulk_names[[tissue]],
                    pval_cutoff = cutoff, n_sc_cell_types = length(sc$cell_types),
                    n_sc_tested_genes = uniqueN(sc$data$phenotype_id),
                    n_bulk_tested_genes_total = nrow(bulk),
                    n_bulk_egenes_total = sum(bulk$bulk_min_pval < cutoff),
                    n_sc_egenes_missing_expression = sum(!is.finite(genes$expr_max)))]
    per_cell <- rbindlist(lapply(sc$cell_types, function(ct) {
        ids <- sc$data[cell_type == ct & sc_is_egene == TRUE, phenotype_id]
        result <- summarize_status(genes[phenotype_id %in% ids])
        result[, `:=`(tissue = tissue, ancestry = "EUR", bulk_source = bulk_names[[tissue]],
                       pval_cutoff = cutoff, cell_type = ct,
                       n_sc_tested_genes = sc$manifest[cell_type == ct, n_sc_tested_genes])]
        result
    }))
    info <- file.info(bulk_file)
    source <- data.table(tissue = tissue, bulk_source = bulk_names[[tissue]],
                         path = bulk_file, bytes = info$size,
                         mtime = as.character(info$mtime))
    list(genes = genes, summary = summary, per_cell = per_cell,
         sc_inputs = sc$manifest, bulk_inputs = source)
}

expression_association_tests <- function(genes) {
    stopifnot(!anyDuplicated(genes[, .(tissue, phenotype_id)]))
    # All genes are sc eGenes: bulk-significant genes are shared; the others are additional.
    comparisons <- list(
        not_identified_vs_identified = list(
            case = c("not_tested", "tested_not_significant"), reference = "significant"),
        not_tested_vs_identified = list(
            case = "not_tested", reference = "significant"),
        tested_not_significant_vs_identified = list(
            case = "tested_not_significant", reference = "significant")
    )
    wilcoxon_rows <- fisher_rows <- list()
    row_id <- 0L
    for (tissue_name in unique(genes$tissue)) {
        d <- genes[tissue == tissue_name]
        for (contrast in names(comparisons)) {
            groups <- comparisons[[contrast]]
            is_case <- d$bulk_status %in% groups$case
            is_reference <- d$bulk_status %in% groups$reference
            metadata <- data.table(tissue = tissue_name, ancestry = "EUR",
                bulk_source = d$bulk_source[1], pval_cutoff = d$pval_cutoff[1],
                contrast = contrast, case_status = paste(groups$case, collapse = ";"),
                reference_status = paste(groups$reference, collapse = ";"))
            for (feature in c("expr_max", "expr_ct_count")) {
                row_id <- row_id + 1L
                values <- d[[feature]]
                x <- values[is_case & is.finite(values)]
                y <- values[is_reference & is.finite(values)]
                n_case <- length(x)
                n_reference <- length(y)
                stat <- p <- probability_lower <- NA_real_
                status <- "missing_group"
                if (n_case > 0L && n_reference > 0L) {
                    if (uniqueN(c(x, y)) == 1L) {
                        stat <- n_case * n_reference / 2
                        p <- 1
                        status <- "constant_feature"
                    } else {
                        test <- wilcox.test(x, y, alternative = "less", exact = FALSE, correct = TRUE)
                        stat <- unname(test$statistic)
                        p <- test$p.value
                        status <- "ok"
                    }
                    # P(case < reference) + 0.5 * P(tie); >0.5 means lower in cases.
                    probability_lower <- 1 - stat / (n_case * n_reference)
                }
                wilcoxon_rows[[row_id]] <- cbind(metadata, data.table(
                    feature = feature, alternative = "less", n_case = n_case,
                    n_reference = n_reference, n_case_excluded = sum(is_case) - n_case,
                    n_reference_excluded = sum(is_reference) - n_reference,
                    median_case = median(x), median_reference = median(y),
                    wilcoxon_W = stat, probability_case_lower = probability_lower,
                    p_value = p, test_status = status))

                low <- if (feature == "expr_max") d$expr_decile == 1L else d$expr_ct_count == 1L
                x_low <- low[is_case & !is.na(low)]
                y_low <- low[is_reference & !is.na(low)]
                counts <- matrix(c(sum(x_low), sum(!x_low), sum(y_low), sum(!y_low)),
                                 nrow = 2L, byrow = TRUE,
                                 dimnames = list(c("case", "reference"), c("low", "other")))
                odds <- ci_low <- ci_high <- p <- NA_real_
                status <- "missing_group"
                if (all(rowSums(counts) > 0L)) {
                    if (any(colSums(counts) == 0L)) {
                        p <- 1
                        status <- "constant_feature"
                    } else {
                        test <- fisher.test(counts, alternative = "two.sided")
                        odds <- unname(test$estimate)
                        ci_low <- test$conf.int[1]
                        ci_high <- test$conf.int[2]
                        p <- test$p.value
                        status <- "ok"
                    }
                }
                fisher_rows[[row_id]] <- cbind(metadata, data.table(
                    feature = feature, low_definition = if (feature == "expr_max")
                        "expr_decile == 1" else "expr_ct_count == 1",
                    alternative = "two.sided", n_case = length(x_low), n_reference = length(y_low),
                    n_case_excluded = sum(is_case) - length(x_low),
                    n_reference_excluded = sum(is_reference) - length(y_low),
                    n_case_low = counts[1, 1], n_case_other = counts[1, 2],
                    n_reference_low = counts[2, 1], n_reference_other = counts[2, 2],
                    fraction_case_low = safe_fraction(sum(x_low), length(x_low)),
                    fraction_reference_low = safe_fraction(sum(y_low), length(y_low)),
                    odds_ratio = odds, odds_ratio_ci_low = ci_low, odds_ratio_ci_high = ci_high,
                    p_value = p, test_status = status))
            }
        }
    }
    adjust <- function(rows, method) {
        x <- rbindlist(rows)
        x[, `:=`(p_adj_BH = p.adjust(p_value, method = "BH"),
                  test_method = method, n_tests_in_family = .N)]
        x
    }
    list(wilcoxon_tests = adjust(wilcoxon_rows, "Wilcoxon rank sum (normal approximation, continuity correction)"),
         fisher_tests = adjust(fisher_rows, "Fisher exact (conditional odds ratio, 95% CI)"))
}

make_figures <- function(summary, stratified, prefix) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Install ggplot2 to generate figures.")
    colors <- c(significant = "#E3B566", tested_not_significant = "#3671AC", not_tested = "#A5C4DD")
    labels <- c(significant = "Significant in bulk", tested_not_significant = "Tested, not significant",
                not_tested = "Not tested in bulk")
    save_plot <- function(plot, suffix, width, height) {
        ggplot2::ggsave(paste0(prefix, suffix, ".png"), plot, width = width, height = height, dpi = 300)
        ggplot2::ggsave(paste0(prefix, suffix, ".pdf"), plot, width = width, height = height)
    }
    counts <- melt(summary, id.vars = "tissue", measure.vars = c(
        "n_bulk_identified", "n_bulk_tested_not_significant", "n_bulk_not_tested"),
        variable.name = "bulk_status", value.name = "count")
    counts[, bulk_status := factor(bulk_status, levels = c("n_bulk_identified",
        "n_bulk_tested_not_significant", "n_bulk_not_tested"), labels = names(colors))]
    counts[, tissue := factor(tissue, levels = summary$tissue)]
    p <- ggplot2::ggplot(counts, ggplot2::aes(tissue, count, fill = bulk_status)) +
        ggplot2::geom_col(width = 0.7) + ggplot2::theme_classic(base_size = 12) +
        ggplot2::scale_fill_manual(values = colors, labels = labels, drop = FALSE) +
        ggplot2::labs(x = NULL, y = "No. single-cell eGenes", fill = NULL,
                       subtitle = "EUR single-cell; nominal p < 5e-8 in single-cell and bulk") +
        ggplot2::theme(legend.position = "bottom")
    save_plot(p, "_counts", 8, 5)
    for (variable in c("expr_decile", "ct_count", "expr_ct_count")) {
        d <- copy(stratified[stratification == variable & !is.na(stratum)])
        if (!nrow(d)) next
        d <- melt(d, id.vars = c("tissue", "stratum"),
                  measure.vars = c("frac_bulk_tested_not_significant", "frac_bulk_not_tested"),
                  variable.name = "bulk_status", value.name = "fraction")
        d[, bulk_status := factor(bulk_status, levels = c("frac_bulk_tested_not_significant",
            "frac_bulk_not_tested"), labels = c("tested_not_significant", "not_tested"))]
        d[, tissue := factor(tissue, levels = summary$tissue)]
        p <- ggplot2::ggplot(d, ggplot2::aes(stratum, fraction, fill = bulk_status)) +
            ggplot2::geom_col(width = 0.8) + ggplot2::facet_wrap(~ tissue, scales = "free_x", ncol = 3) +
            ggplot2::scale_y_continuous(labels = function(x) paste0(round(x * 100), "%"), limits = c(0, 1)) +
            ggplot2::scale_fill_manual(values = colors, labels = labels) +
            ggplot2::theme_classic(base_size = 11) + ggplot2::theme(legend.position = "bottom") +
            ggplot2::labs(x = switch(variable,
                expr_decile = "Max expression decile",
                expr_ct_count = "No. cell types with expression data for the gene",
                ct_count = "No. single-cell types in which the gene was tested"),
                y = "eGenes not identified in bulk (%)", fill = NULL)
        if (variable == "expr_decile") p <- p + ggplot2::scale_x_continuous(breaks = 1:10)
        save_plot(p, paste0("_", variable), 12, 6)
    }
    make_publication_figure(summary, stratified, prefix)
}

make_publication_figure <- function(summary, stratified, prefix, width = 7.1, height = 5.2) {
    for (package in c("ggplot2", "cowplot", "ragg")) {
        if (!requireNamespace(package, quietly = TRUE)) stop("Install ", package, " to generate the combined figure.")
    }
    family <- "Helvetica"
    colors <- c(significant = "#E3B566", tested_not_significant = "#3671AC", not_tested = "#A5C4DD")
    labels <- c(significant = "Significant in bulk", tested_not_significant = "Tested, not significant",
                not_tested = "Not tested in bulk")
    tissue_order <- as.character(summary$tissue)
    common_theme <- ggplot2::theme_classic(base_size = 8, base_family = family) +
        ggplot2::theme(
            text = ggplot2::element_text(family = family, color = "black"),
            axis.text = ggplot2::element_text(size = 7.5, color = "black"),
            axis.title = ggplot2::element_text(size = 8, color = "black"),
            axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 3)),
            axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 4)),
            axis.line = ggplot2::element_line(linewidth = 0.3, color = "black"),
            axis.ticks = ggplot2::element_line(linewidth = 0.3, color = "black"),
            axis.ticks.length = grid::unit(1, "mm"),
            strip.background = ggplot2::element_blank(),
            strip.text = ggplot2::element_text(size = 8, face = "bold", margin = ggplot2::margin(b = 3)),
            panel.spacing.x = grid::unit(2.2, "mm"),
            legend.title = ggplot2::element_blank(),
            legend.text = ggplot2::element_text(size = 8,
                margin = ggplot2::margin(l = 2, r = 6)),
            legend.key.size = grid::unit(3, "mm"),
            legend.key.spacing.x = grid::unit(4, "mm"),
            legend.margin = ggplot2::margin(4, 4, 6, 4),
            plot.margin = ggplot2::margin(t = 12, r = 5, b = 12, l = 7)
        )
    fill_scale <- function() ggplot2::scale_fill_manual(values = colors, breaks = names(colors),
        labels = labels, drop = FALSE)
    percentage_axis <- function() ggplot2::scale_y_continuous(breaks = c(0, 0.5, 1),
        labels = c("0", "50", "100"), limits = c(0, 1),
        expand = ggplot2::expansion(mult = c(0, 0.02)))

    proportions <- melt(summary, id.vars = c("tissue", "n_sc_egenes"), measure.vars = c(
        "n_bulk_identified", "n_bulk_tested_not_significant", "n_bulk_not_tested"),
        variable.name = "bulk_status", value.name = "count")
    proportions[, bulk_status := factor(bulk_status, levels = c("n_bulk_identified",
        "n_bulk_tested_not_significant", "n_bulk_not_tested"), labels = names(colors))]
    proportions[, `:=`(tissue = factor(tissue, levels = tissue_order), proportion = count / n_sc_egenes)]
    p_a <- ggplot2::ggplot(proportions, ggplot2::aes(1, proportion, fill = bulk_status)) +
        ggplot2::geom_col(width = 0.2) + ggplot2::facet_wrap(~ tissue, nrow = 1) +
        ggplot2::scale_x_continuous(limits = c(0.45, 1.55), breaks = NULL,
            expand = ggplot2::expansion(mult = 0)) +
        percentage_axis() + fill_scale() + common_theme +
        ggplot2::labs(x = NULL, y = "Single-cell eGenes (%)") +
        ggplot2::theme(legend.position = "top", legend.direction = "horizontal",
            axis.line.x = ggplot2::element_blank()) +
        ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1))
    stratified_panel <- function(variable) {
        d <- copy(stratified[stratification == variable & !is.na(stratum)])
        if (!nrow(d)) stop("No strata for combined figure: ", variable)
        d <- melt(d, id.vars = c("tissue", "stratum"),
            measure.vars = c("frac_bulk_tested_not_significant", "frac_bulk_not_tested"),
            variable.name = "bulk_status", value.name = "proportion")
        d[, `:=`(tissue = factor(tissue, levels = tissue_order),
            bulk_status = factor(bulk_status, levels = c("frac_bulk_tested_not_significant",
                "frac_bulk_not_tested"), labels = c("tested_not_significant", "not_tested")))]
        p <- ggplot2::ggplot(d, ggplot2::aes(stratum, proportion, fill = bulk_status)) +
            ggplot2::geom_col(width = 0.85) +
            ggplot2::facet_wrap(~ tissue, nrow = 1, scales = "free_x") +
            percentage_axis() + fill_scale() + common_theme +
            ggplot2::theme(legend.position = "none") +
            ggplot2::labs(x = if (variable == "ct_count") "No. cell types tested" else
                "Maximum expression decile", y = "sc-only eGenes (%)")
        if (variable == "expr_decile") {
            p <- p + ggplot2::scale_x_continuous(breaks = c(1, 5, 10),
                expand = ggplot2::expansion(add = 0.5))
        } else {
            p <- p + ggplot2::scale_x_continuous(breaks = function(limits) {
                maximum <- floor(limits[2])
                unique(c(1, seq(10, max(10, maximum), 10)))
            }, expand = ggplot2::expansion(add = 0.5))
        }
        p
    }
    figure_prefix <- paste0(prefix, "_combined")
    # Measure text on the final raster device, not the IDE/default PDF device:
    # device-specific Helvetica metrics otherwise make legend entries overlap.
    ragg::agg_png(paste0(figure_prefix, ".png"), width = width, height = height,
                  units = "in", res = 300, background = "white")
    output_device <- grDevices::dev.cur()
    figure <- tryCatch({
        # Recent ggplot2 versions include empty guide boxes for unused positions.
        plot_grob <- ggplot2::ggplotGrob(p_a)
        legend_candidates <- plot_grob$grobs[grepl("^guide-box", plot_grob$layout$name)]
        legend_candidates <- Filter(function(x) inherits(x, "gtable"), legend_candidates)
        if (length(legend_candidates) != 1L) stop("Could not find the shared bulk-status legend.")
        legend <- legend_candidates[[1]]
        panels <- cowplot::plot_grid(
            p_a + ggplot2::theme(legend.position = "none"),
            stratified_panel("ct_count"), stratified_panel("expr_decile"),
            ncol = 1, align = "v", axis = "lr", rel_heights = c(0.85, 1, 1),
            labels = c("A", "B", "C"), label_fontfamily = family,
            label_size = 10, label_fontface = "bold", label_x = 0.005, label_y = 0.98,
            hjust = 0, vjust = 1
        )
        figure <- cowplot::plot_grid(legend, panels, ncol = 1, rel_heights = c(0.085, 0.915))
        print(figure)
        figure
    }, finally = grDevices::dev.off(output_device))
    proportions[, tissue := as.character(tissue)]
    fwrite(proportions, paste0(figure_prefix, "_proportions.tsv"), sep = "\t", na = "NA") 
    message("Publication figure saved: ", figure_prefix)
    invisible(figure)
}

main <- function(tissues = names(bulk_names), project_dir = ".",
                 output_dir = file.path(project_dir, "crossTissue_analysis_R", "scatterplot"),
                 gtex_dir = Sys.getenv("GTEX_DIR", file.path(project_dir, "references", "GTEx_v8")),
                 eqtlgen_file = Sys.getenv("EQTLGEN_FILE", file.path(project_dir, "references", "eQTLGen", eqtlgen_filename)),
                 workers = 1L, plots = TRUE) {
    cutoff <- 5e-8
    if (!length(tissues) || any(!tissues %in% names(bulk_names)) || anyDuplicated(tissues)) {
        stop("Choose distinct tissues from: ", paste(names(bulk_names), collapse = ", "))
    }
    if (length(workers) != 1L || is.na(workers) || workers < 1 || workers != as.integer(workers)) {
        stop("workers must be a positive integer")
    }
    setDTthreads(1L)
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    cache_dir <- file.path(output_dir, "bulk_eGene_cache")
    compare <- function(tissue) compare_tissue(tissue, project_dir, gtex_dir, eqtlgen_file, cache_dir, cutoff)
    results <- if (workers == 1L) lapply(tissues, compare) else
        parallel::mclapply(tissues, compare, mc.cores = min(workers, length(tissues)), mc.preschedule = FALSE)
    if (any(vapply(results, inherits, logical(1), "try-error"))) stop("A tissue comparison failed; see worker errors.")
    names(results) <- tissues
    combined <- lapply(names(results[[1]]), function(part) rbindlist(lapply(results, `[[`, part)))
    names(combined) <- names(results[[1]])
    genes <- combined$genes
    combined$stratified <- rbindlist(lapply(c("expr_decile", "ct_count", "expr_ct_count"), function(column) {
        x <- genes[, summarize_status(.SD), by = c("tissue", "ancestry", "bulk_source", "pval_cutoff", column)]
        setnames(x, column, "stratum")
        x[, stratification := column]
        x
    }))
    combined <- c(combined, expression_association_tests(genes))
    # Keep subset runs separate from the complete five-tissue output.
    tag <- if (setequal(tissues, names(bulk_names))) "" else paste0("_", paste(tissues, collapse = "_"))
    prefix <- file.path(output_dir, paste0("bulk_eGene_EUR_p5e-8", tag))
    for (part in names(combined)) {
        x <- combined[[part]]
        first <- intersect(c("tissue", "ancestry", "bulk_source", "pval_cutoff", "cell_type", "phenotype_id"), names(x))
        setcolorder(x, c(first, setdiff(names(x), first)))
        fwrite(x, paste0(prefix, "_", part, ".tsv"), sep = "\t", na = "NA")
    }
    for (tissue_name in tissues) {
        saveRDS(genes[tissue == tissue_name], file.path(output_dir, paste0(tissue_name, "_",
            bulk_names[[tissue_name]], "_eGene_EUR_p5e-8.rds")))
    }
    saveRDS(list(parameters = list(ancestry = "EUR", pval_cutoff = cutoff, nPC = 5L,
        sc_af_range = c(0.01, 0.99), tissues = tissues, created = Sys.time()),
        tables = combined), paste0(prefix, ".rds"))
    if (plots) make_figures(combined$summary, combined$stratified, prefix)
    print(combined$summary[, .(tissue, n_sc_egenes, n_bulk_identified, n_bulk_not_identified,
                               n_bulk_not_tested, n_bulk_tested_not_significant)])
    message("Results saved: ", prefix)
    invisible(combined)
}

if (sys.nframe() == 0L) {
    args <- commandArgs(trailingOnly = TRUE)
    if ("--help" %in% args) {
        cat("Usage: Rscript script_replication_validation/egene_bulk_comparison_public.R [--workers=N] [Blood Lung Skin Liver Colon]\n",
            "Defaults: all five tissues, one worker, EUR only, nominal p < 5e-8.\n",
            "Run from the project root. Set EQTLGEN_FILE or GTEX_DIR to override the reference paths.\n",
            "Results: crossTissue_analysis_R/scatterplot/bulk_eGene_EUR_p5e-8*\n", sep = "")
    } else {
        worker_arg <- grep("^--workers=", args, value = TRUE)
        if (length(worker_arg) > 1L) stop("Specify --workers once")
        workers <- if (length(worker_arg)) as.numeric(sub("^--workers=", "", worker_arg)) else 1L
        tissues <- args[!grepl("^--workers=", args)]
        main(tissues = if (length(tissues)) tissues else names(bulk_names), workers = workers)
    }
}
