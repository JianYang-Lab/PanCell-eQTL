# ==============================================================================
# R Script for creating feature matrix directory upon dropEst output
# ==============================================================================

library(Matrix)
library(tidyverse)

args = commandArgs(TRUE)
tissue=args[1]
project=args[2]
sample=args[3]

setwd(file.path("/path/to/",tissue,"raw_data",project))

data=readRDS(paste0(sample,"/",sample,".dropEst.rds"))

system(paste0("mkdir -p ", sample, "/outs/filtered_feature_bc_matrix"))

m = data$cm
writeMM(obj = m, file = paste0(sample,"/outs/filtered_feature_bc_matrix/matrix.mtx"))

# Write barcodes (cell identifiers)
write.table(colnames(m), file = paste0(sample,"/outs/filtered_feature_bc_matrix/barcodes.tsv"), sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

# Write features (gene identifiers), assuming the feature type is 'Gene Expression'
gene_gtf = read.table("/path/to/resource/genes.gtf", header=TRUE)
colnames(gene_gtf)[1:2] = c("feature_name","feature_id")
features <- data.frame(feature_name = rownames(m)) %>% left_join(gene_gtf[,1:2],by="feature_name")
features$feature_type = "Gene Expression"
write.table(features[,c(2,1,3)], file = paste0(sample,"/outs/filtered_feature_bc_matrix/features.tsv"), sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

system(paste0("gzip ", sample,"/outs/filtered_feature_bc_matrix/matrix.mtx"))
system(paste0("gzip ", sample,"/outs/filtered_feature_bc_matrix/barcodes.tsv"))
system(paste0("gzip ", sample,"/outs/filtered_feature_bc_matrix/features.tsv"))

print("Session completed")
