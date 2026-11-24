#!/usr/bin/env Rscript
# Prepare EUR cross-tissue mashr inputs. Run from the project root or supply its path.
# Usage: Rscript mashr_cross_tissue_prepare_public.R LINEAGE [B|T|NK] [project_root] [output_dir]
# Supply a subgroup only for Lymphoid. Use a fresh output directory for new input data.
# Inputs: crossTissue_analysis_R/tissue_ancestry_celltype_ss_cc_colors.txt and, per tissue,
# sc-eQTL_R/genotype/*_EUR_indivpruned.bim and sc-eQTL_R/results/{sig,raw} tensorQTL files.
# The sig directory contains eQTL_summary_5en8.tsv and <cell>_EUR_pc5_sig.tsv.gz.
# The raw directory contains <cell>_EUR_pc5_perm.tsv and nominal pairs for chromosomes 1–22.

suppressPackageStartupMessages(library(data.table))

cross_tissue_args <- function(args) {
  usage <- "Expected LINEAGE [B|T|NK] [project_root] [output_dir]; e.g. Lymphoid NK or Myeloid"
  if (!length(args)) stop(usage, call. = FALSE)
  group <- args[1]
  args <- args[-1]
  if (group == "Lymphoid") {
    if (!length(args) || !args[1] %in% c("B", "T", "NK")) stop(usage, call. = FALSE)
    group <- args[1]
    args <- args[-1]
  } else if (!group %in% c("Myeloid", "Epithelial", "Endothelial", "Stromal", "Other")) {
    stop(usage, call. = FALSE)
  }
  if (length(args) > 2L) stop(usage, call. = FALSE)
  root <- if (length(args)) args[1] else "."
  out_dir <- if (length(args) == 2L) args[2] else file.path(root, "crossTissue_analysis_R", "mashR", "cross_tissue")
  list(group = group, root = root, out_dir = out_dir)
}

# A completed checkpoint becomes visible after its write succeeds.
write_checkpoint <- function(object, path, writer = saveRDS) {
  tmp <- tempfile(tmpdir = dirname(path))
  on.exit(unlink(tmp))
  writer(object, tmp)
  if (!file.rename(tmp, path)) stop("Could not save ", path)
}

cross_tissue_contexts <- function(root, tissues, group) {
  meta <- fread(file.path(root, "crossTissue_analysis_R", "tissue_ancestry_celltype_ss_cc_colors.txt"))
  meta <- meta[ancestry == "EUR" & sample_size >= 50 & total_count >= 10000]
  meta[lineage %in% c("Erythroid", "Megakaryocytic", "Neural", "Progenitor"), lineage := "Other"]
  subtypes <- list(
    B = c("BIN", "BMemNS", "BMemS", "BMem", "Bf", "Plasma", "Plasmablasts"),
    T = c("CD4TEM", "CD4TEMRA", "CD4TNC", "CD4TNCInf", "CD8TEMRA", "CD8TNC", "CD8TRM",
          "gdT", "MAIT", "CD8TEMTRM", "Treg", "Th1", "Th17"),
    NK = c("NKp", "NK")
  )
  if (group %in% names(subtypes)) {
    meta <- meta[lineage == "Lymphoid" & cell_type_short %in% subtypes[[group]]]
  } else {
    meta <- meta[lineage == group]
  }
  # Preserve summary-file order: the first eligible Lung context supplies the random pairs.
  contexts <- rbindlist(lapply(tissues, function(ts) {
    results <- file.path(root, ts, "sc-eQTL_R", "results")
    summary <- fread(file.path(results, "sig", "eQTL_summary_5en8.tsv"))
    cells <- unique(summary[ancestry == "EUR" & cell_type %in% meta[tissue == ts, cell_type], cell_type])
    data.table(tissue = rep(ts, length(cells)), cell_type = cells,
               prefix = file.path(results, "raw", paste0(cells, "_EUR_pc5", recycle0 = TRUE)))
  }))
  contexts[file.exists(paste0(prefix, ".cis_qtl_pairs.22.parquet"))]
}

nominal_pairs <- function(prefix, chr) {
  path <- paste0(prefix, ".cis_qtl_pairs.", chr, ".parquet")
  arrow::open_dataset(path, format = "parquet")
}

prepare_cross_tissue <- function(group, root, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  checkpoint <- function(step) file.path(out_dir, paste0(step, "_", group, ".rds"))
  prepared <- checkpoint("step4_strong_random_for_mash")
  if (file.exists(prepared)) return(prepared)

  tissues <- c("Blood", "Lung", "Skin", "Liver", "Colon")
  contexts <- cross_tissue_contexts(root, tissues, group)
  if (!nrow(contexts)) stop("No eligible EUR contexts for ", group)
  fields <- c("phenotype_id", "variant_id", "slope", "slope_se", "pval_nominal")
  keys <- fields[1:2]
  cutoff <- 5e-8

  # Variants present in the EUR genotype reference for all five tissues.
  bim <- file.path(out_dir, "step0_overlap_EUR.bim")
  if (file.exists(bim)) {
    overlap <- fread(bim, header = FALSE)[[1]]
  } else {
    for (ts in tissues) {
      path <- file.path(root, ts, "sc-eQTL_R", "genotype", paste0(ts, "_EUR_indivpruned.bim"))
      ids <- fread(path, header = FALSE, select = 2L)[[1]]
      overlap <- if (ts == tissues[1]) unique(ids) else intersect(overlap, ids)
    }
    write_checkpoint(data.table(variant_id = sort(overlap, method = "radix")), bim,
                     function(x, path) fwrite(x, path, col.names = FALSE, quote = FALSE))
  }

  top_file <- checkpoint("step1_top_all")
  genes_file <- checkpoint("step1_gene_list")
  if (file.exists(top_file) && file.exists(genes_file)) {
    top <- readRDS(top_file)
    genes <- readRDS(genes_file)
  } else {
    top <- genes <- vector("list", nrow(contexts))
    for (i in seq_len(nrow(contexts))) {
      prefix <- contexts$prefix[i]
      perm <- fread(paste0(prefix, "_perm.tsv"))
      genes[[i]] <- unique(perm$phenotype_id)
      path <- file.path(dirname(dirname(prefix)), "sig", paste0(basename(prefix), "_sig.tsv.gz"))
      candidates <- fread(path, select = fields)
      # Recover pairs that can pass the nominal cutoff but fall outside the gene's FDR export.
      thresholds <- if ("pval_nominal_threshold" %in% names(perm)) perm else {
        fread(path, select = c("phenotype_id", "pval_nominal_threshold"))
      }
      tight_genes <- unique(thresholds[pval_nominal_threshold < cutoff, phenotype_id])
      if (length(tight_genes)) {
        extra <- lapply(1:22, function(chr) {
          nominal_pairs(prefix, chr) |>
            dplyr::filter(phenotype_id %in% tight_genes, pval_nominal < cutoff) |>
            dplyr::select(tidyselect::all_of(fields)) |>
            dplyr::collect()
        })
        candidates <- rbindlist(c(list(candidates), extra))
      }
      lead <- candidates[variant_id %in% overlap & !is.na(pval_nominal),
                         .SD[which.min(pval_nominal)], by = phenotype_id]
      lead[, `:=`(tissue = contexts$tissue[i], ancestry = "EUR")]
      top[[i]] <- lead
    }
    top <- rbindlist(top)
    write_checkpoint(top, top_file)
    write_checkpoint(genes, genes_file)
  }

  coverage <- table(unlist(genes, use.names = FALSE))
  keep_genes <- names(coverage)[coverage == max(coverage)]
  strong_file <- checkpoint("step2_top_all_cond")
  if (file.exists(strong_file)) {
    strong <- readRDS(strong_file)
  } else {
    strong <- top[phenotype_id %in% keep_genes, .SD[which.min(pval_nominal)], by = phenotype_id]
    strong <- strong[pval_nominal < cutoff]
    write_checkpoint(strong, strong_file)
  }

  # Sample up to 50 shared variants per gene in the first eligible Lung context.
  random_file <- checkpoint("step4_random_list")
  if (file.exists(random_file)) {
    random <- readRDS(random_file)
    # The original checkpoint stores one list of gene tables per chromosome.
    if (!is.data.frame(random)) random <- rbindlist(unlist(random, recursive = FALSE))
  } else {
    lung <- contexts[tissue == "Lung", prefix]
    if (!length(lung)) stop("No eligible Lung EUR context for random sampling in ", group)
    set.seed(1)
    random <- rbindlist(lapply(1:22, function(chr) {
      x <- nominal_pairs(lung[1], chr) |>
        dplyr::select(tidyselect::all_of(keys)) |>
        dplyr::collect() |> as.data.table()
      x <- x[phenotype_id %in% keep_genes & variant_id %in% overlap]
      x[, .(variant_id = variant_id[sample.int(.N, min(.N, 50L))]), by = phenotype_id]
    }))
    write_checkpoint(random, random_file)
  }

  # Rows are SNP-gene pairs; columns follow the original tissue/ancestry/cell-type order.
  pairs <- unique(rbindlist(list(random[, ..keys], strong[, ..keys])))
  setkeyv(pairs, keys)
  lookup <- arrow::Table$create(pairs[, .(phenotype_id, variant_id, pair_index = .I)])
  contexts[, condition := paste(tissue, "EUR", cell_type, sep = "_")]
  setorder(contexts, condition)
  b <- se <- z <- matrix(NA_real_, nrow(pairs), nrow(contexts), dimnames = list(NULL, contexts$condition))
  for (i in seq_len(nrow(contexts))) {
    message("Reading ", contexts$condition[i])
    for (chr in 1:22) {
      x <- nominal_pairs(contexts$prefix[i], chr) |>
        dplyr::inner_join(lookup, by = keys) |>
        dplyr::select(pair_index, slope, slope_se, pval_nominal) |>
        dplyr::collect()
      idx <- x$pair_index
      beta <- x$slope
      beta[is.na(beta)] <- 0
      b[idx, i] <- beta
      se[idx, i] <- x$slope_se
      signed_z <- sign(beta) * abs(qnorm(x$pval_nominal / 2))
      signed_z[is.na(signed_z)] <- 0
      z[idx, i] <- signed_z
    }
  }

  # Strict overlap removes any pair with a missing or unusable SE in any context.
  complete <- rowSums(!is.finite(se) | se <= 0 | se == 1000) == 0L
  extract <- function(selected) {
    idx <- sort(pairs[selected, on = keys, which = TRUE])
    idx <- idx[complete[idx]]
    if (!length(idx)) stop("No complete SNP-gene pairs remain for ", group)
    list(pair = as.data.frame(pairs[idx]), z = as.data.frame(z[idx, , drop = FALSE]),
         b = as.data.frame(b[idx, , drop = FALSE]), s = as.data.frame(se[idx, , drop = FALSE]))
  }
  out <- c(setNames(extract(random), paste0("random_", c("pair", "z", "b", "s"))),
           setNames(extract(strong), paste0("strong_", c("pair", "z", "b", "s"))))
  write_checkpoint(out, prepared)
  prepared
}

if (sys.nframe() == 0L) {
  message(do.call(prepare_cross_tissue, cross_tissue_args(commandArgs(trailingOnly = TRUE))))
}
