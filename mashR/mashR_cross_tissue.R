#!/usr/bin/env Rscript
# Fit cross-tissue mashr to EUR single-cell eQTLs for one lineage.
# Usage: Rscript mashr_cross_tissue_public.R LINEAGE [B|T|NK] [project_root] [output_dir]
# Examples from the project root:
#   Rscript script_mashr/mashr_cross_tissue_public.R Lymphoid NK
#   Rscript script_mashr/mashr_cross_tissue_public.R Myeloid
# The sibling preparation script can also be run on its own with the same arguments.
# Outputs default to crossTissue_analysis_R/mashR/cross_tissue.

script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
source(file.path(dirname(normalizePath(script)), "mashR_cross_tissue_prepare.R"))
args <- cross_tissue_args(commandArgs(trailingOnly = TRUE))
suppressPackageStartupMessages(library(mashr))
out <- readRDS(do.call(prepare_cross_tissue, args))
result_path <- function(step) file.path(args$out_dir, paste0(step, "_", args$group, ".rds"))

# Estimate the null correlation and fit mixture weights on the random pairs.
random <- mash_set_data(as.matrix(out$random_b), as.matrix(out$random_s))
V <- estimate_null_correlation_simple(random)
random <- mash_set_data(as.matrix(out$random_b), as.matrix(out$random_s), V = V)
strong <- mash_set_data(as.matrix(out$strong_b), as.matrix(out$strong_s), V = V)
set.seed(1)
n_pc <- if (args$group == "NK") 3L else 5L
U_pca <- cov_pca(strong, min(n_pc, ncol(out$strong_b), nrow(out$strong_b)))
U <- c(cov_ed(strong, U_pca), cov_canonical(random))
fit <- mash(random, Ulist = U, outputlevel = 1)
write_checkpoint(fit, result_path("step5_m"))

# Apply the fitted mixture to strong pairs and estimate sharing within a factor of two.
posterior <- mash(strong, g = get_fitted_g(fit), fixg = TRUE)
write_checkpoint(posterior, result_path("step5_m2"))
sharing <- get_pairwise_sharing(posterior, lfsr_thresh = 0.05, factor = 0.5)
write_checkpoint(sharing, result_path("step5_pairwise_PM"))
