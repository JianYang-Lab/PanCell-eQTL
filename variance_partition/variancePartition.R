# Load necessary libraries
library(data.table)
library(parallel)
library(BiocParallel)
suppressMessages(library(tidyverse))
suppressMessages(library(variancePartition))

# 0. Create the gene expression matrix
args <- commandArgs(TRUE)
ts <- args[1]
anc <- args[2]
ct <- args[3]

# 1. Define File Paths
bed_file <- file.path(ts, "sc-eQTL_R/celltype_matrix/bed/", paste0(ct,"_",anc,"_expr.bed"))
cov_file <- file.path(ts, "sc-eQTL_R/covariates/final", paste0(ct,"_",anc,"_cov_final.txt"))
cellcount_file <- file.path(ts,"/h5ad_celltype_R/", paste0(ts,"_rawQC_",ct,"_",anc,"_cov.csv"))

# 2. Read the Normalized Expression BED File
bed_data <- fread(bed_file, data.table = FALSE)

# Extract just the gene IDs and the numeric sample matrix
rownames(bed_data) <- bed_data$TargetID 
expr_mat <- as.matrix(bed_data[, 5:ncol(bed_data)])

# 3. Read the Covariate Matrix
info <- fread(cov_file, data.table=FALSE)

# Ensure the first column is named consistently for merging later
colnames(info)[1] <- "project_sample"
rownames(info) <- info$project_sample

# 4. Regress out the covariates for each gene in expr_mat
# Define covariates to regress
col_to_resgress <- c("Age_pred", "Gender_pred", "PC1", "PC2", "PC3", "PC4", "PC5", 
                     colnames(info)[str_starts(colnames(info), "PF")], 
                     colnames(info)[str_starts(colnames(info), "Condition_group")])

cov_regress <- info[, col_to_resgress]
rownames(cov_regress) <- info$project_sample

# Align samples between expression matrix and covariate matrix before regression
common_samples_reg <- intersect(colnames(expr_mat), rownames(cov_regress))
expr_mat_aligned <- expr_mat[, common_samples_reg]
cov_regress_aligned <- cov_regress[common_samples_reg, ]

# Use lm.fit for fast matrix-based linear regression
# Create a design matrix (adding an intercept term implicitly via formula)
design_mat <- model.matrix(~ ., data = cov_regress_aligned)

# Fit the model (t() transposes expr_mat so genes are columns, which lm.fit expects)
fit <- lm.fit(design_mat, t(expr_mat_aligned))

# Extract residuals and transpose back (rows = genes, cols = samples)
# This is your new expression matrix with the covariates regressed out
expr_mat_regressed <- t(fit$residuals) 

# 5. Format Variables for the LMM
# Convert dummy variables for Condition_group back to Condition_group
info <- info %>% 
  pivot_longer(
    cols = starts_with("Condition_group"),
    names_to = "Condition_group", 
    names_prefix = "Condition_group", 
    values_to = "dummy_value"
  ) %>%
  filter(dummy_value == 1) %>%
  select(-dummy_value) %>%
  mutate(Condition_group = factor(Condition_group))

# Read cell counts and merge
cc <- fread(cellcount_file, data.table=FALSE)
# Left join drops rownames, so we must reassign them
info <- as.data.frame(info) %>% left_join(cc, by="project_sample")
rownames(info) <- info$project_sample

# Other variables 
info$Batch <- sapply(rownames(info), function(x) str_split(x, "_")[[1]][1])
info$Condition_group <- as.factor(info$Condition_group)
info$log10_cell_count <- log10(info$cell_count)
info$log10_n_genes <- log10(info$avg_n_genes)
info$log10_total_counts <- log10(info$avg_total_counts)
info$log10_Age_pred <- log10(info$Age_pred)

# 6. Run Variance Partitioning
## varPart on original expression matrix
form <- ~ (1|Batch) + (1|Condition_group) + log10_cell_count + log10_n_genes + log10_total_counts + log10_Age_pred + Gender_pred + PC1 + PC2 + PC3 + PC4 + PC5
varPart0 <- fitExtractVarPartModel(expr_mat_aligned, form, info, BPPARAM = MulticoreParam(workers = 4))
varPart0 <- sortCols(varPart0)
out_dir <- file.path(ts, "/sc-eQTL_R/results/varPart")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE) 
write.table(varPart0, file.path(out_dir, paste0("varPart_orig_expr_mat_", ts, "_", ct, "_", anc, ".txt")), sep='\t', row.names=FALSE, quote=FALSE, col.names=TRUE)

## varPart on regressed expression matrix
if(ts=="Lung"){
  info$Smoking_status_pred = ifelse(info$Smoking_status_pred=="Smoker",1,0)
  form <- ~ (1|Batch) + Smoking_status_pred + log10_cell_count + log10_n_genes + log10_total_counts 
}else{
  form <- ~ (1|Batch) + log10_cell_count + log10_n_genes + log10_total_counts 
}
varPart <- fitExtractVarPartModel(expr_mat_regressed, form, info, BPPARAM = MulticoreParam(workers = 4))
varPart <- sortCols(varPart)
write.table(varPart, file.path(out_dir, paste0("varPart_regressed_expr_mat_", ts, "_", ct, "_", anc, ".txt")), sep='\t', row.names=FALSE, quote=FALSE, col.names=TRUE)

# 7. Visualize the Results
# plotVarPart(varPart)
# ggsave("tmp.png", width=8, height=6)