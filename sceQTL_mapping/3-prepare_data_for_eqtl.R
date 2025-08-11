##############################################################################
# Script information
# Step 3: Prepare data for eQTL 
# Description: Prepare covariates and expression bed for eqtl
##############################################################################

library(tidyverse)


# Set the arguments to call from the bash script -----
args = commandArgs(TRUE) 
tissue = args[1]
celltype = args[2]
ancestry = args[3]
setwd(paste0("/path/to/", tissue, "/sc-eQTL/"))
# Get covariate files ---------------------------------
# Peer factors
factor_path = file.path("covariates", "peer_factors", 
                paste0(celltype, "_",ancestry, "_peer_factors.tsv"))
factors_df = read.csv(factor_path, sep='\t')

# Demographic covariates
cov_path = paste0("covariates/",tissue,"_all_sample_cov.txt")
covs = read.csv(cov_path, sep="\t")

# PCs
eig_path = file.path("genotype", paste0(tissue, "_", ancestry, "_indivpruned.eigenval"))
vec_path = file.path("genotype", paste0(tissue, "_", ancestry, "_indivpruned.eigenvec"))
eig = sqrt(read.table(eig_path)[,1])
vec = read.table(vec_path)
nPCS <- ncol(vec)-2
for(k in 1:nPCS){
	vec[,2+k] <- vec[,2+k] * eig[k]
}
vec = vec[,-1]
colnames(vec) = c("sample",paste("PC",1:20, sep="")) 

# Check all the individuals exist in both files
table(factors_df$sample %in% covs$sample)
table(factors_df$sample %in% vec$sample)
sample_isec = intersect(intersect(factors_df$sample,covs$sample),vec$sample)

# Merge all covariates
factors_cov = factors_df %>% dplyr::filter(sample %in% sample_isec) %>% 
	left_join(covs, by='sample') %>% 
        left_join(vec, by='sample')

factors_cov$Gender_pred = sapply(factors_cov$Gender_pred, function(x) ifelse(x=='female', 0, 1)) # F=0, M=1
factors_cov$Condition_group = as.factor(factors_cov$Condition_group)
dummy_matrix <- model.matrix(~ Condition_group - 1, data = factors_cov) # Create dummy variables for health conditions
factors_cov = cbind(factors_cov[,colnames(factors_cov)!="Condition_group"], as.data.frame(dummy_matrix))

out_path = file.path("covariates", "final", 
                paste0(celltype, "_", ancestry, "_cov_final.txt"))
write.table(factors_cov, file=out_path, sep='\t', col=T, row=F, quo=F)

# Create bed ----------------------------------------
expr_path = file.path("celltype_matrix", "norm_rint_tsv", 
                paste0(tissue,"_expr_norm_rint_", ancestry, "_",celltype, ".tsv"))
expr_df = read.csv(expr_path, sep='\t', check.names = FALSE)
colnames(expr_df)[1] = "gene_name"
expr_df = expr_df[,colnames(expr_df) %in% c("gene_name",sample_isec)]
n_samples = dim(expr_df)[2] - 1

# merge gene information
gene_bed = read.csv("/path/to/resource/genes.bed", sep='\t') %>%
            dplyr::mutate(chr=str_replace(chr, "chr", ""))
expr_df = expr_df %>% 
    left_join(gene_bed, by="gene_name") %>%
    dplyr::select(!gene_name) %>%
    dplyr::arrange(chr,start) 
#Chr    start   end     TargetID
bed_df_out = expr_df[,c((n_samples+1):(n_samples+4), 1:n_samples)]
colnames(bed_df_out)[1:4] = c("#Chr", "start", "end", "TargetID")
bed_df_out = bed_df_out[!is.na(bed_df_out$TargetID),]
bed_out_path = file.path("celltype_matrix", "bed",
                paste0(celltype, "_", ancestry, "_expr.bed"))
write.table(bed_df_out, file=bed_out_path, sep='\t', col=T, row=F, quo=F)

