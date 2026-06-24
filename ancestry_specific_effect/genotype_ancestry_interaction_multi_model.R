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

# 0. Helper function for multi-SNP model -----------------------
run_gxe_model_multi <- function(expr_vec, geno_df, pop_vec, covar_matrix, top_snp_col) {
  df <- data.frame(
    expression = as.numeric(expr_vec),
    ancestry = as.factor(pop_vec) 
  )
  # Combine expression, ancestry, multiple genotypes, and covariates
  df <- cbind(df, geno_df, covar_matrix)
  df <- na.omit(df)
  
  # Ensure sufficient variation
  if(nrow(df) < 50 || length(unique(df$ancestry)) < 2) return(NA)
  
  # Protect column names with backticks
  covar_str <- paste(sprintf("`%s`", colnames(covar_matrix)), collapse = " + ")
  snp_str <- paste(sprintf("`%s`", colnames(geno_df)), collapse = " + ")
  
  # Null Model: Expression ~ Ancestry + All Lead SNPs + Covariates
  form_null <- as.formula(paste("expression ~ ancestry +", snp_str, "+", covar_str))
  
  # Full Model: Null Model + Interaction(Ancestry x Top SNP)
  form_full <- as.formula(paste(
    "expression ~ ancestry +", 
    snp_str, "+", 
    sprintf("ancestry:`%s`", top_snp_col), "+", 
    covar_str
  ))
  
  # Use tryCatch to prevent crashing if a model is rank-deficient/aliased
  model_null <- tryCatch(lm(form_null, data = df), error = function(e) NULL)
  model_full <- tryCatch(lm(form_full, data = df), error = function(e) NULL)
  
  if (is.null(model_null) || is.null(model_full)) return(NA)
  
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
# Get one lead SNP per credible set
susie_df <- susie_df[order(-pip), .SD[1], by=.(phenotype_id, cs_id)] 

# 2. Get the genotype -------------------
variant_list_file <- file.path(eqtl_dir, "raw", paste0(ct, "_GLB_susie_multi_variants_to_extract.txt"))
out_prefix <- file.path(eqtl_dir, "raw", paste0(ct, "_GLB_multi_subset_dosages"))
  
fwrite(data.table(V1 = unique(susie_df$variant_id)), variant_list_file, col.names = FALSE, quote = FALSE)
  
plink_cmd <- sprintf("plink2 --bfile %s --extract %s --export A --out %s", 
                       POOLED_PLINK_PREFIX, variant_list_file, out_prefix)
print("Extracting dosages with PLINK2...")
system(plink_cmd, ignore.stdout = TRUE) 
  
raw_file <- paste0(out_prefix, ".raw")
if(!file.exists(raw_file)) {
    print(paste("WARNING: PLINK extraction failed for", ct, "- skipping."))
    quit(status = 0) 
}
geno_matrix <- fread(raw_file)

# 3. Load expression and covariates for testing -------------------
bed_file <- file.path(ts, "sc-eQTL_R", "celltype_matrix", "bed", paste0(ct, "_GLB_expr.bed"))
if(!file.exists(bed_file)) stop(paste("Expression BED file missing:", bed_file))
expr_merged <- fread(bed_file)
setnames(expr_merged, "TargetID", "phenotype_id")

cov_file <- file.path(ts, "sc-eQTL_R", "covariates", "final", paste0(ct, "_GLB_cov_final.txt"))
if(!file.exists(cov_file)) stop(paste("Covariate file missing:", cov_file))
cov_merged <- fread(cov_file)
setnames(cov_merged, colnames(cov_merged), gsub(" ", "_", colnames(cov_merged)))


# 4. Run Interactions -------------------
print("GxE interaction tests start...")

# Identify genes to test (Note: N>=2 restricts to multi-signal genes. Change to N>=1 to test all eGenes)
egene_list <- susie_df[, .N, by=.(phenotype_id)][N >= 2, phenotype_id]

res_list <- mclapply(1:length(egene_list), function(j) {
    gene <- egene_list[j]
    
    # Get all lead SNPs for this gene, sorted by highest PIP to determine the "top" SNP
    gene_snps <- susie_df[phenotype_id == gene][order(-pip)]
    top_snp <- gene_snps$variant_id[1]
    all_snps <- gene_snps$variant_id
    
    # Safely match all PLINK column names (e.g., rs123_A)
    regex_pattern <- paste0("^(", paste(all_snps, collapse="|"), ")_")
    valid_snp_cols <- grep(regex_pattern, colnames(geno_matrix), value = TRUE)
    
    if (length(valid_snp_cols) == 0) return(NULL) 
    
    # Identify the specific column name for the Top SNP
    top_snp_col <- grep(paste0("^", top_snp, "_"), valid_snp_cols, value = TRUE)
    if (length(top_snp_col) == 0) return(NULL) # Skip if the top SNP dosage is missing
    
    # Align Samples
    valid_samples <- intersect(colnames(expr_merged)[-c(1:4)], geno_matrix$IID)
    valid_samples <- intersect(valid_samples, cov_merged$sample)
    if(length(valid_samples) < 50) return(NULL) 
    
    # Extract Data
    y_raw <- expr_merged[phenotype_id == gene, ..valid_samples]
    if(nrow(y_raw) == 0) return(NULL)
    y <- as.numeric(unlist(y_raw))
    
    # Get multi-SNP matrix
    g_df <- geno_matrix[match(valid_samples, geno_matrix$IID), ..valid_snp_cols]
    
    meta_sub <- cov_merged[match(valid_samples, cov_merged$sample), ]
    pop <- meta_sub$Ancestry_pred
    
    # Covariates
    col_select <- c(colnames(meta_sub)[str_detect(colnames(meta_sub), "PF")],
                    paste0("PC", 1:5),
                    colnames(meta_sub)[str_detect(colnames(meta_sub), "Condition_group")],
                    "Age_pred", "Gender_pred")
    col_select <- intersect(col_select, colnames(meta_sub))
    covars <- meta_sub[, ..col_select]
    
    # Run Model
    pval <- tryCatch({ run_gxe_model_multi(y, g_df, pop, covars, top_snp_col) }, error = function(e) NA)
    
    return(data.table(phenotype_id = gene, 
                      top_variant_id = top_snp, 
                      n_snps_in_model = length(valid_snp_cols),
                      interaction_pval = pval))
  }, mc.cores = 6)

# 5. Save Results -------------------
final_res <- rbindlist(res_list[!sapply(res_list, is.null)])
final_res <- final_res[!is.na(interaction_pval)]

if(nrow(final_res) > 0) {
    final_res[, pval_adj := p.adjust(interaction_pval, method="bonferroni")]
    out_dir <- file.path(ts, "sc-eQTL_R", "results", "interaction_tests")
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    out_file <- file.path(out_dir, paste0(ct, "_multiSNP_gxe_interaction.tsv"))
    
    fwrite(final_res, out_file, sep = "\t", quote = FALSE)
    print(paste("Saved:", out_file))
} else {
    print("No valid interaction results generated.")
}