#!/usr/bin/env Rscript
# Prepare allele-aligned AMR eQTLs for Peruvian tensorQTL analyses.
# Optional: Rscript Peruvian_benchmark_prepare.R PROJECT_ROOT OUTPUT.rds
# The plotting script also calls prepare_eqtls() if the prepared table is unavailable.

suppressPackageStartupMessages(library(data.table))

cell_labels <- c(CD4TNC = "CD4TNC", CD4TEM = "CD4TEM", CD8TRM = "CD8TEMTRM",
                 CD8TEMRA = "CD8TEMRA", Treg = "Treg")
platforms <- c("PeruvianWGS", "PeruvianscRNA")  # SNP-array and scRNA-seq, respectively
cutoff <- 5e-8

prepare_eqtls <- function(root) {
  results <- file.path(root, platforms, "sc-eQTL_R", "results")
  keys <- c("phenotype_id", "variant_id")
  tables <- lapply(names(cell_labels), function(ct) {
    message("Reading AMR eQTLs for ", ct)
    prefix <- paste0(ct, "_AMR_pc5")
    nominal <- lapply(results, function(path) {
      files <- file.path(path, "raw", paste0(prefix, ".cis_qtl_pairs.", 1:22, ".parquet"))
      arrow::open_dataset(files, format = "parquet")
    })
    significant <- lapply(seq_along(results), function(i) {
      sig <- fread(file.path(results[i], "sig", paste0(prefix, "_sig.tsv.gz")),
                   select = c(keys, "pval_nominal", "pval_nominal_threshold"))
      if (!nrow(sig) || any(sig$pval_nominal_threshold < cutoff, na.rm = TRUE)) {
        return(nominal[[i]] |>
          dplyr::filter(pval_nominal < cutoff) |>
          dplyr::select(tidyselect::all_of(keys)) |> dplyr::collect())
      }
      sig[pval_nominal < cutoff, ..keys]
    })
    # Fetch effects in both analyses for every pair significant in either one.
    pairs <- arrow::Table$create(unique(rbindlist(significant)))
    rbindlist(lapply(seq_along(platforms), function(i) {
      x <- nominal[[i]] |> dplyr::inner_join(pairs, by = keys) |> dplyr::collect() |> as.data.table()
      x[, c("Chr", "BP", "A2", "A1") := tstrsplit(variant_id, "_", keep = 1:4)]
      x[, `:=`(tissue = platforms[i], ancestry = "AMR", cell_type = ct)]
      x
    }))
  })

  all_pairs <- rbindlist(tables)
  rbindlist(lapply(platforms, function(platform) {
    path <- file.path(root, platform, "sc-eQTL_R", "genotype", paste0(platform, "_AMR_indivpruned.bim"))
    bim <- fread(path, header = FALSE, select = c(1, 4, 5, 6), colClasses = "character")
    x <- all_pairs[tissue == platform]
    by_x <- c("Chr", "BP", "A2", "A1")
    direct <- merge(x, bim, by.x = by_x, by.y = c("V1", "V4", "V6", "V5"))
    reverse <- merge(x, bim, by.x = by_x, by.y = c("V1", "V4", "V5", "V6"))
    reverse[, slope := -slope]
    rbind(direct, reverse)
  }))
}

# Match alleles across platforms. For panel G, use the union of their lead positions.
matched_effects <- function(eqtls, ct, lead_only = FALSE) {
  fields <- c("phenotype_id", "Chr", "BP", "A2", "A1", "pval_nominal", "slope")
  array <- eqtls[cell_type == ct & tissue == "PeruvianWGS", ..fields]
  scrna <- eqtls[cell_type == ct & tissue == "PeruvianscRNA", ..fields]
  position <- c("phenotype_id", "Chr", "BP")
  if (lead_only) {
    lead_positions <- function(x) x[order(pval_nominal), .SD[1L], by = phenotype_id][, ..position]
    leads <- unique(rbind(lead_positions(array), lead_positions(scrna)))
    array <- merge(array, leads, by = position)
    scrna <- merge(scrna, leads, by = position)
  }
  keys <- c(position, "A2", "A1")
  direct <- merge(array, scrna, by = keys, suffixes = c("_array", "_scrna"))
  reverse <- merge(array, scrna, by.x = keys, by.y = c(position, "A1", "A2"),
                   suffixes = c("_array", "_scrna"))
  reverse[, slope_scrna := -slope_scrna]
  rbind(direct, reverse)
}

# Compute r² between distinct lead variants using the SNP-array reference genotypes.
lead_ld <- function(leads, bfile) {
  leads <- copy(leads)
  leads[, r2 := fifelse(variant_array == variant_scrna, 1, NA_real_)]
  different <- unique(leads[variant_array != variant_scrna, .(variant_array, variant_scrna)])
  if (!nrow(different)) return(leads)
  work <- tempfile("s23-plink-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE))
  snps <- file.path(work, "variants.txt")
  prefix <- file.path(work, "dosages")
  writeLines(unique(c(different$variant_array, different$variant_scrna)), snps)
  plink <- Sys.which("plink")
  if (!nzchar(plink)) stop("Panel J requires PLINK 1.9 on PATH.")
  log <- system2(plink, c("--bfile", shQuote(bfile), "--extract", shQuote(snps),
                          "--recode", "A", "--out", shQuote(prefix)), stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(log, "status"))) stop(paste(tail(log, 5L), collapse = "\n"))
  dosage <- fread(paste0(prefix, ".raw"))
  # PLINK appends the counted allele to each dosage column, after the six sample fields.
  ids <- sub("_[^_]+$", "", names(dosage)[-(1:6)])
  different[, r2 := mapply(function(a, b) {
    if (is.na(a) || is.na(b)) return(NA_real_)
    cor(dosage[[a + 6L]], dosage[[b + 6L]], use = "pairwise.complete.obs")^2
  }, match(variant_array, ids), match(variant_scrna, ids))]
  leads[different, on = .(variant_array, variant_scrna), r2 := i.r2]
  leads
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 2L) stop("Usage: Rscript Peruvian_benchmark_prepare_public.R PROJECT_ROOT OUTPUT.rds")
  dir.create(dirname(args[2]), recursive = TRUE, showWarnings = FALSE)
  saveRDS(prepare_eqtls(args[1]), args[2])
}
