# Prepare the two cohort-level points in compact final-plot panel E.
# Uses data.table, loaded by CIMA_benchmark_common.R.
cima_egene_recovery_data <- function(root, p_cutoff = 5e-8) {
  pairs <- data.table(reference = c("CIMAWGS", "CIMAWGSsub"),
                      cohort = c("CIMAscRNA", "CIMAscRNAsub"))
  suffix <- "_EAS_pc5_perm.tsv"
  # Hold the eGene comparison scope fixed across full and downsampled cohorts.
  memory_t_cells <- sort(c("CD4_CTL", "CD4_T_Helper", "CD4_Tcm", "CD4_Tem",
                          "CD4_Tr1", "CD4_Treg", "CD8_CTL", "CD8_Tcm",
                          "CD8_Tem", "Proliferative_T", "gdT1", "gdT2"))
  inputs <- donor_tables <- summary_tables <- matched_tables <- list()
  record_input <- function(path, cohort, cell_type, kind) {
    info <- file.info(path)
    inputs[[length(inputs) + 1L]] <<- data.table(
      cohort = cohort, cell_type = cell_type, kind = kind,
      path = normalizePath(path, mustWork = TRUE), size_bytes = info$size,
      modified_at = format(info$mtime, "%Y-%m-%dT%H:%M:%S%z"))
  }
  permutation_files <- function(cohort) {
    directory <- file.path(root, cohort, "sc-eQTL_R/results/raw")
    paths <- list.files(directory, pattern = "_EAS_pc5_perm\\.tsv$", full.names = TRUE)
    if (!length(paths)) stop("No permutation tables in ", directory)
    setNames(paths, substr(basename(paths), 1L, nchar(basename(paths)) - nchar(suffix)))
  }
  significant_gene_cells <- function(paths, cohort) {
    genes <- lapply(names(paths), function(cell) {
      path <- paths[[cell]]
      x <- fread(path, select = c("phenotype_id", "pval_nominal"), na.strings = c("NA", "NaN", ""))
      if (!nrow(x) || anyNA(x$phenotype_id) || any(!nzchar(x$phenotype_id)) || anyDuplicated(x$phenotype_id))
        stop("Empty table or missing/duplicate phenotype IDs: ", path)
      if (!is.numeric(x$pval_nominal) || any(!is.na(x$pval_nominal) &
          (!is.finite(x$pval_nominal) | x$pval_nominal < 0 | x$pval_nominal > 1)))
        stop("Invalid nominal P values: ", path)
      record_input(path, cohort, cell, "permutation")
      x[is.finite(pval_nominal) & pval_nominal < p_cutoff,
        .(cell_type = cell, phenotype_id)]
    })
    rbindlist(genes)
  }
  donor_counts <- function(cohort, cells) {
    fam_path <- file.path(root, cohort, "sc-eQTL_R/genotype", paste0(cohort, "_EAS_indivpruned.fam"))
    fam <- fread(fam_path, header = FALSE, select = 2L, col.names = "individual_id")
    if (!nrow(fam) || anyNA(fam$individual_id) || anyDuplicated(fam$individual_id))
      stop("Missing/duplicate genotype donors: ", fam_path)
    record_input(fam_path, cohort, NA_character_, "genotype_fam")
    tables <- lapply(cells, function(cell) {
      path <- file.path(root, cohort, "h5ad_celltype_R",
                        paste0(cohort, "_rawQC_", cell, "_EAS_cov.csv"))
      x <- fread(path, select = c("project_sample", "cell_count"))
      if (!nrow(x) || anyNA(x$project_sample) || anyDuplicated(x$project_sample) ||
          !is.numeric(x$cell_count) || any(!is.finite(x$cell_count) | x$cell_count < 0 |
                                         x$cell_count != floor(x$cell_count)))
        stop("Invalid per-donor cell counts: ", path)
      record_input(path, cohort, cell, "expression_cell_counts")
      x[project_sample %in% fam$individual_id, .(individual_id = project_sample, cell_count)]
    })
    counts <- rbindlist(tables)[, .(total_cells = sum(cell_count)), by = individual_id]
    if (!setequal(fam$individual_id, counts$individual_id))
      stop("Some genotype donors have no expression cell counts: ", cohort)
    counts[, cohort := cohort]
    setcolorder(counts, c("cohort", "individual_id", "total_cells"))
    counts[order(individual_id)]
  }
  for (i in seq_len(nrow(pairs))) {
    reference <- pairs$reference[i]
    cohort <- pairs$cohort[i]
    ref_paths <- permutation_files(reference)
    query_paths <- permutation_files(cohort)
    if (!setequal(names(ref_paths), names(query_paths)))
      stop("Cell-type coverage differs between ", reference, " and ", cohort)
    cells <- memory_t_cells
    if (!all(cells %in% names(query_paths)))
      stop("Missing matched memory T-cell mappings for ", cohort, ": ",
           paste(setdiff(cells, names(query_paths)), collapse = ", "))
    ref_gene_cells <- significant_gene_cells(ref_paths[cells], reference)
    query_gene_cells <- significant_gene_cells(query_paths[cells], cohort)
    ref_genes <- unique(ref_gene_cells$phenotype_id)
    query_genes <- unique(query_gene_cells$phenotype_id)
    if (!length(ref_genes)) stop("No reference eGenes: ", reference)
    # Match significance within cell type before deduplicating genes. Genes
    # significant only in different cell types are not recovered.
    matched <- merge(ref_gene_cells, query_gene_cells, by = c("cell_type", "phenotype_id"))
    n_recovered <- uniqueN(matched$phenotype_id)
    matched[, `:=`(cohort = cohort, reference = reference)]
    setcolorder(matched, c("cohort", "reference", "cell_type", "phenotype_id"))
    matched_tables[[i]] <- matched
    # Preserve the total-cell x-axis definition; record this scope separately
    # from the matched memory T-cell scope used to calculate eGene recovery.
    count_cells <- sort(names(query_paths))
    counts <- donor_counts(cohort, count_cells)
    donor_tables[[i]] <- counts
    summary_tables[[i]] <- data.table(
      cohort = cohort, reference = reference, n_cell_types = length(cells),
      n_cell_types_for_cell_counts = length(count_cells),
      n_individuals = nrow(counts), total_cells = sum(counts$total_cells),
      median_cells_per_donor = median(counts$total_cells),
      mean_cells_per_donor = mean(counts$total_cells),
      n_wgs_unique_egenes = length(ref_genes), n_scrna_unique_egenes = length(query_genes),
      n_recovered_unique_same_cell = n_recovered,
      egene_recovery_percent = 100 * n_recovered / length(ref_genes),
      pval_column = "pval_nominal", p_threshold = p_cutoff,
      cell_types = paste(cells, collapse = ";"),
      cell_count_cell_types = paste(count_cells, collapse = ";"))
  }
  list(summary = rbindlist(summary_tables), donors = rbindlist(donor_tables),
       inputs = rbindlist(inputs), matched_gene_cells = rbindlist(matched_tables))
}
