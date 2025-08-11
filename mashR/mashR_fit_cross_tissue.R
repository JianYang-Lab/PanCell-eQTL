##############################################################################
# Fit mashR model to all SNP-gene pairs for a given cell type across tissues 
# by chromosome.
##############################################################################

library(tidyverse)
library(data.table)
library(ashr)
library(mashr)
library(rmeta)
library(arrow)
library(parallel)
library(circlize)
library(corrplot)

# --------------------------------------------------------
# 0. Preparing ancestry, tissue, cell type list 
# --------------------------------------------------------
tissue_list = c("Blood","Lung","Skin","Colon","Liver")
ancestry_list = c("EUR", "EAS", "AFR", "AMR")
ancestry_color_list = c(EUR="#66C2A5",EAS="#FC8D62",AFR="#8DA0CB",AMR="#E78AC3")

ct_cross_tissue=read_delim("/path/to/crossTissue_analysis/cross_tissue_celltype.txt",delim="\t")

mashR_path = "/path/to/crossTissue_analysis/mashR/cross_tissue/" 

args = commandArgs(TRUE)
ct = args[1]
chrom = as.numeric(args[2])

## Fit all
ConvertP2Z <- function(pval, beta) {
        z <- abs(qnorm(pval / 2))
        z[which(beta < 0)] <- -1 * z[which(beta < 0)]
        return(z)
}

handle_nan_b = function(x) {
        x[which(is.nan(x) | is.na(x))] = 0
        return(x)
}

handle_nan_s = function(x) {
        x[which(is.nan(x) | is.infinite(x) | is.na(x) | x == 0)] = 1E3
        return(x)
}

handle_nan_z = function(z, b, s) {
        z[which(is.nan(s) | is.infinite(s) | is.na(s) | s == 0)] = 0
        z[which(is.nan(b) | is.na(b) | b == 0)] = 0
        z[which(is.nan(z) | is.na(z))] = 0
        return(z)
}

GetSS <- function(dat) {
        dat$"z" <- ConvertP2Z(dat$"pval_nominal", dat$"slope")
        dat[['z']] = handle_nan_z(dat[['z']], dat[['slope']], dat[['slope_se']])
        dat[['slope_se']] = handle_nan_s(dat[['slope_se']])
        dat[['slope']] = handle_nan_b(dat[['slope']])
        return(dat)
}

process_combinations <- function(ct_anc_combo) {
  anc <- ct_anc_combo$anc
  tissue <- ct_anc_combo$tissue
  ct <- ct_anc_combo$ct
  chrom <- ct_anc_combo$chrom
  tissue_path = paste0("/path/to/",tissue,"/sc-eQTL/results/")
  ct_anc_tmp <- paste0(ts,"_",ct)
    file_path <- file.path(tissue_path, "raw", paste0(ct, "_", anc, "_pc5.cis_qtl_pairs.", chrom, ".parquet"))
    if (file.exists(file_path)) {
      print(paste0(tissue, "_", ct))
      df <- read_parquet(file_path)
      df <- df[df$phenotype_id %in% gene_keep & df$variant_id %in% geno_overlap, ]
      df <- df[df$af >= 0.01 & df$af <= 0.99,]
      df$ancestry <- anc
      df$tissue <- tissue
      return(df)
    }
  
  return(NULL)  # Return NULL if the combination does not meet criteria
}


if (! file.exists(paste0(mashR_path, "step5_summary_statistics_",ct,"_chr",chrom,".rds"))){
  geno_overlap = read.table(paste0(mashR_path, "step0_overlap_EUR.bim"))$V1
  gene_list = readRDS(paste0(mashR_path, "step1_gene_list_", ct, ".rds"))
  count = table(unlist(gene_list))
  n_max = max(count)
  gene_keep = names(count)[count==n_max]

  print(paste0("Read in raw data for ",ct," chr",chrom, "."))

  # Create a list of ancestry and cell type combinations
  ct_anc_combinations <- list()
  anc="EUR"
    for (ts in tissue_list) {
      ct_raw = ct_cross_tissue %>% dplyr::filter(tissue==ts & celltype==ct) %>% pull(celltype_raw)
      ct_anc_combinations[[length(ct_anc_combinations) + 1]] <- list(
        anc = anc,
        ct = ct_raw,
        tissue = ts,
        chrom = chrom
      )
    }
    # Use mclapply to process combinations in parallel
  df_all_list <- mclapply(ct_anc_combinations, function(x) process_combinations(x), mc.cores=2)

  # Filter out NULL results
  df_all_list <- Filter(Negate(is.null), df_all_list)

  # Get summary statistics
  df_all_ss = mclapply(df_all_list, GetSS, mc.cores=6)
  df_all_ss = rbindlist(df_all_ss)

  df_all_z_wide <- dcast(df_all_ss[,.(phenotype_id,variant_id,z,ancestry,tissue)], phenotype_id+variant_id~tissue+ancestry, value.var=c("z"))
  df_all_slope_wide <- dcast(df_all_ss[,.(phenotype_id,variant_id,slope,ancestry,tissue)], phenotype_id+variant_id~tissue+ancestry, value.var=c("slope"))
  df_all_se_wide <- dcast(df_all_ss[,.(phenotype_id,variant_id,slope_se,ancestry,tissue)], phenotype_id+variant_id~tissue+ancestry, value.var=c("slope_se"))

  ncond=ncol(df_all_z_wide)
  df_all_pair = df_all_z_wide[,1:2]
  saveRDS(df_all_pair, paste0(mashR_path, "step5_all_pairs_",ct,"_chr",chrom,".rds"))
  saveRDS(list(df_all_slope_wide=df_all_slope_wide, df_all_se_wide=df_all_se_wide),  paste0(mashR_path, "step5_summary_statistics_",ct,"_chr",chrom,".rds"))
  }else{
  summary_stats <- readRDS(paste0(mashR_path, "step5_summary_statistics_",ct,"_chr",chrom,".rds"))
  df_all_slope_wide <- summary_stats$df_all_slope_wide
  df_all_se_wide <- summary_stats$df_all_se_wide
  ncond=ncol(df_all_slope_wide) 
}


out = readRDS(paste0(mashR_path, "step4_strong_random_for_mash_",ct,".rds"))
data.temp = mash_set_data(Bhat=as.matrix(out$random_b), Shat=as.matrix(out$random_s))
Vhat = estimate_null_correlation_simple(data.temp)
rm(data.temp)
data_STRONG = mash_set_data(as.matrix(df_all_slope_wide[,3:ncond]), as.matrix(df_all_se_wide[,3:ncond]), V = Vhat)

print("Read in fitted model")
m = readRDS(paste0(mashR_path, "step5_m_",ct,".rds"))

print("Computed posterior summaries")
m2 = mash(data_STRONG, g = get_fitted_g(m), fixg = TRUE)
saveRDS(m2, paste0(mashR_path, "step5_m2_all_posterior_",ct,"_chr",chrom,".rds"))

