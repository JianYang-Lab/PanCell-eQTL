library(data.table)
suppressMessages(library(tidyverse))
library(parallel)

setwd("/path/to")
tissue_list = c("Blood","Lung","Skin","Liver","Colon")
ancestry_list = c("EUR", "EAS", "AFR", "AMR")
ancestry_color_list = c(EUR="#66C2A5",EAS="#FC8D62",AFR="#8DA0CB",AMR="#E78AC3")
ancestries_all_list = list(Blood=c("EUR", "EAS", "AFR", "AMR"),
                           Lung=c("EUR","EAS"),
                           Skin=c("EUR"),
                           Colon=c("EUR"),
                           Liver=c("EUR","EAS"))

args <- commandArgs(TRUE)
if(length(args) < 2) stop("Please provide tissue (ts) and celltype (ct) arguments.")
ts <- args[1]
ct <- args[2]

# Genotype from the global population
POOLED_PLINK_PREFIX <- file.path(ts, "sc-eQTL_R", "genotype", paste0(ts,"_GLB_indivpruned"))

# 0. Helper function -----------------------
run_gxe_model <- function(expr_vec, geno_vec, pop_vec, covar_matrix) {
  df <- data.frame(
    expression = as.numeric(expr_vec),
    genotype = as.numeric(geno_vec),
    ancestry = as.factor(pop_vec) 
  )
  df <- cbind(df, covar_matrix)
  df <- na.omit(df)
  
  if(length(unique(df$genotype)) < 2 || length(unique(df$ancestry)) < 2) return(NA)
  
  covar_names <- colnames(covar_matrix)
  # Wrap covariate names in backticks to prevent formula evaluation errors (e.g., if names have hyphens)
  covar_str <- paste(sprintf("`%s`", covar_names), collapse = " + ")
  
  form_null <- as.formula(paste("expression ~ genotype + ancestry +", covar_str))
  form_full <- as.formula(paste("expression ~ genotype * ancestry +", covar_str))
  
  model_null <- lm(form_null, data = df)
  model_full <- lm(form_full, data = df)
  
  lrt <- anova(model_null, model_full)
  return(lrt$`Pr(>F)`[2])
}

# 1. Get the pairs to test ---------------------
print("Select gene-SNP pairs to test...")
eqtl_dir <- file.path(ts, "sc-eQTL_R", "results")
susie_path <- file.path(eqtl_dir, "raw", paste0(ct, "_GLB_pc5_susie.tsv"))
perm_path <- file.path(eqtl_dir, "raw", paste0(ct,"_GLB_pc5_perm.tsv"))

if(!file.exists(susie_path) || !file.exists(perm_path)) {
    stop("Input susie or perm files are missing.")
}

susie_df <- fread(susie_path)
egene <- unique(susie_df$phenotype_id)
perm_df <- fread(perm_path)

pairs_to_test <- perm_df[phenotype_id %in% egene, .(phenotype_id,variant_id)]

# 2. Get the genotype -------------------
variant_list_file <- file.path(eqtl_dir, "raw", paste0(ct, "_GLB_susie_variants_to_extract.txt"))
out_prefix <- file.path(eqtl_dir, "raw", paste0(ct, "_GLB_subset_dosages"))
  
fwrite(data.table(V1 = unique(pairs_to_test$variant_id)), variant_list_file, col.names = FALSE, quote = FALSE)
  
plink_cmd <- sprintf("plink2 --bfile %s --extract %s --export A --out %s", 
                       POOLED_PLINK_PREFIX, variant_list_file, out_prefix)
print("Extracting dosages with PLINK2...")
system(plink_cmd, ignore.stdout = TRUE) 
  
raw_file <- paste0(out_prefix, ".raw")
if(!file.exists(raw_file)) {
    print(paste("WARNING: PLINK extraction failed for", ct, "- skipping."))
    quit(status = 0) # 'next' is invalid here, exit gracefully instead
}
geno_matrix <- fread(raw_file)

# 3. Load expression and covariates for testing -------------------
bed_file <- file.path(ts, "sc-eQTL_R", "celltype_matrix", "bed", paste0(ct, "_GLB_expr.bed"))
if(file.exists(bed_file)) {
    expr_merged <- fread(bed_file)
} else {
    stop(paste("Expression BED file missing:", bed_file))
}
setnames(expr_merged, "TargetID", "phenotype_id")

cov_file <- file.path(ts, "sc-eQTL_R", "covariates", "final", paste0(ct, "_GLB_cov_final.txt"))
if(file.exists(cov_file)) {
    cov_merged <- fread(cov_file)
} else {
    stop(paste("Covariate file missing:", cov_file))
}

# Vectorized approach to replace spaces in column names (much faster than a loop)
setnames(cov_merged, colnames(cov_merged), gsub(" ", "_", colnames(cov_merged)))

print("GxE interaction tests start...")
res_list <- mclapply(1:nrow(pairs_to_test), function(j) {
    gene <- pairs_to_test$phenotype_id[j]
    snp <- pairs_to_test$variant_id[j]
    
    # PLINK appends the minor allele to the column name (e.g., "rs123_A")
    snp_col_idx <- grep(paste0("^", snp, "_"), colnames(geno_matrix))
    if (length(snp_col_idx) == 0) return(NULL) 
    
    # Align Samples
    valid_samples <- intersect(colnames(expr_merged)[-c(1:4)], geno_matrix$IID)
    valid_samples <- intersect(valid_samples, cov_merged$sample)
    if(length(valid_samples) < 50) return(NULL) # Skip if too few shared samples
    
    # Extract & safely unlist the expression data to avoid list-coercion errors
    y_raw <- expr_merged[phenotype_id == gene, ..valid_samples]
    if(nrow(y_raw) == 0) return(NULL)
    y <- as.numeric(unlist(y_raw))
    
    # Match order properly
    g_temp <- geno_matrix[match(valid_samples, geno_matrix$IID), ]
    g <- as.numeric(g_temp[[snp_col_idx]])
    
    meta_sub <- cov_merged[match(valid_samples, cov_merged$sample), ]
    pop <- meta_sub$Ancestry_pred
    
    # Extract only the top 5 PCs for the model to prevent overfitting
    pf_cols <- colnames(meta_sub)[str_detect(colnames(meta_sub), "PF")]
    pc_cols <- paste0("PC", 1:5)
    cond_cols <- colnames(meta_sub)[str_detect(colnames(meta_sub), "Condition_group")] 
    other_cols <- c("Age_pred","Gender_pred")
    col_select <- c(pf_cols, pc_cols, cond_cols, other_cols)
    
    # Crucial: Intersect to ensure we don't request columns that don't exist
    col_select <- intersect(col_select, colnames(meta_sub))
    covars <- meta_sub[, ..col_select]
    
    pval <- tryCatch({ run_gxe_model(y, g, pop, covars) }, error = function(e) NA)
    
    # Return a 1-row data.table so rbindlist works cleanly
    return(data.table(phenotype_id = gene, variant_id = snp, interaction_pval = pval))
  }, mc.cores = 6)

final_res <- rbindlist(res_list[!sapply(res_list, is.null)])
final_res <- final_res[!is.na(interaction_pval)]

if(nrow(final_res) > 0) {
    final_res[, pval_adj := p.adjust(interaction_pval, method="bonferroni")]
    out_dir <- file.path(ts, "sc-eQTL_R", "results", "interaction_tests")
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    out_file <- file.path(out_dir, paste0(ct, "_gxe_interaction.tsv"))
    
    fwrite(final_res, out_file, sep = "\t", quote = FALSE)
    print(paste("Saved:", out_file))
} else {
    print("No valid interaction results generated.")
}