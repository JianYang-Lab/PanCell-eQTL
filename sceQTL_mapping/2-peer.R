##############################################################################
# Script information
# Step 2: Identify PEER factors
# Description: PEER factors are identified on the normalized pseudobulk matrix
##############################################################################

# Note that you should conduct peer analysis in the peer environment 

# Import libraries
library(peer)

# Set the arguments to call from the bash script -----
args = commandArgs(TRUE) # nolint
tissue = args[1]
celltype = args[2]
ancestry = args[3]
setwd(paste0("/path/to/", tissue, "/sc-eQTL/"))


# Get expression file --------------------------------
expr_path = paste0("celltype_matrix/norm_rint_tsv/",tissue,"_expr_norm_rint_", ancestry, "_", celltype, ".tsv")
expr = read.csv(expr_path, sep="\t", check.names = FALSE)
rownames(expr) = expr[,1]
expr = t(expr[,-1])
sample_filter = read.csv(paste0("genotype/",tissue,"_",ancestry,"_indivpruned.fam"), sep='\t', header=FALSE)$V1
expr = expr[rownames(expr) %in% sample_filter,]
dim(expr)
samples = rownames(expr)
genes = colnames(expr)

# Perform PEER on RINT matrix --------------------------------
# Set PEER paramaters based on the instructions from PEER package website

model = PEER()

PEER_setPhenoMean(model, as.matrix(expr))

dim(PEER_getPhenoMean(model))


# PEER_setAdd_mean(model, TRUE)

# PEER_setNmax_iterations(model, 100)

# GTEx Nk recommendation: For eQTL analyses, the number of PEER factors was determined as 
# function of sample size (N): 15 factors for N<150, 30 factors for 150<N<250, 45 factors 
# for 250<N<350, and 60 factors for N>350
nsample <- dim(expr)[1]
if(nsample < 150 && nsample > 0){
     n_peer = 15
}else if(nsample >= 150 && nsample < 250){ 
     n_peer = 30
}else if(nsample >= 250 && nsample < 350){
     n_peer = 45
}else if (nsample >= 350) {
     n_peer = 60
}

PEER_setNk(model,n_peer) # Set to generate 60 PEER factors

PEER_getNk(model)

PEER_update(model)

# Calculate and save the PEER factors
factors = PEER_getX(model)
dim(factors)
factors_df = data.frame(factors)
colnames(factors_df) = paste0("PF", 1:n_peer)
factors_df$sample = samples
factors_df = factors_df[,c((n_peer+1),1:n_peer)]

write.table(factors_df, file=file.path("covariates", "peer_factors", paste0(celltype, "_",ancestry, "_peer_factors.tsv")),
     sep="\t", col=T, quo=F, row=F)

# Calculate and save the weights for each factor
weights = PEER_getW(model)
dim(weights)
weights_df = data.frame(weights)
colnames(weights_df) = paste0("PF", 1:n_peer)
weights_df$geneid = genes
weights_df = weights_df[,c((n_peer+1),1:n_peer)]
write.table(weights_df, file=file.path("covariates", "peer_factors", paste0(celltype, "_",ancestry, "_peer_factor_weights.tsv")),
     sep="\t", col=T, quo=F, row=F)

# Calculate and save the precision values
precision = PEER_getAlpha(model)
dim(precision)
precision_df = data.frame(precision)
precision_df$covariate = paste0("PF", 1:n_peer) 
write.table(precision_df, file=file.path("covariates", "peer_factors", paste0(celltype, "_",ancestry, "_peer_factor_precision.tsv")),
     sep="\t", col=T, quo=F, row=F)


# Calculate and save the residuals
# residuals = PEER_getResiduals(model)
# dim(residuals)
# residuals_df = data.frame(residuals)
# colnames(residuals_df) = genes 
# residuals_df = cbind(sample=samples, residuals_df)
# write.table(residuals_df, file=file.path("covariates", "peer_factors", paste0(celltype, "_",ancestry, "_peer_factor_precision.tsv")),
#      sep="\t", col=T, quo=F, row=F)

