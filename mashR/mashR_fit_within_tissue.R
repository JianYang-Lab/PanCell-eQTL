##############################################################################
# Fit mashR model to all SNP-gene pairs for in each tissue by chromosome.
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
ancestries_all_list = list(Blood=c("EUR", "EAS", "AFR", "AMR"),
                            Lung=c("EUR","EAS"),
                            Skin=c("EUR","EAS"),
                            Colon=c("EUR"),
                            Liver=c("EUR","EAS"))
# Color palettes for cell types -------------------------
celltypes_all_list = list(Blood=c('CD4TNC','CD4TEM','Treg','CD8TNC','CD8TEM','CD8TEMRA','MAIT','NKp','NKn','BIN','BMem','Plasma','Plasmablasts','MonoC','MonoNC','Nph','DC','pDC'),
Liver = c('CD4T',"CD8T","Treg","MAIT","Th","gdT",'Circulating_NK','Resident_NK','B','Plasma','Monocytes','Macrophages','DC','pDC','Neutrophils','Basophils','Hepatocytes','Cholangiocytes','Endothelium','Fibroblasts'),
Skin = c("Th","Tc", "Treg","NK","DC1","DC2","MigDC","Macro1","Macro2","MonoMac","Mast","KCdiff","KCundiff","Melanocyte","VE1","VE2","VE3","LE1","Pericyte1","Pericyte2","F2"),
Lung = c("CD4T","CD8T","NK","B","Monocytes","Mast","Macrophages","DC","AT1","AT2","Airway_Epi_Multiciliated","Airway_Epi_Basal","Airway_Epi_Secretory","EC_arterial","EC_capillary","EC_venous","LEC_mature","LEC_diff","Fibroblasts","Smooth_muscle","Mesothelium","Rare"),
Colon = c("CD4T","CD8T","Treg","Th","gdT","NKT","NK","ILC3","BIN","BMem","Bcyc","Plasma","Mono","Macro","cDC2","Mast","Colonocyte","GLoblet","Tuft","TA", "EEC" ,"ECcap","ECven","Stromal1","Stromal2","Myofibroblast","Glia"))

cell_type_cat_colors = read_delim("/path/to/crossTissue_analysis/cell_type_colors_new.txt", delim='\t') 
cell_type_short_dict = lapply(tissue_list,function(x) {
  tmp=cell_type_cat_colors$cell_type_short[cell_type_cat_colors$tissue==x]
  names(tmp) = cell_type_cat_colors$cell_type[cell_type_cat_colors$tissue==x]
  tmp 
  })
names(cell_type_short_dict) = tissue_list
cell_type_color_dict = cell_type_cat_colors$cell_type_colors
names(cell_type_color_dict) = cell_type_cat_colors$cell_type_short

smr_list = list()
for (t in tissue_list){
    smr = read_table(file.path("/path/to/",t,"sc-eQTL","results","sig","eQTL_summary_5en8_maf01.tsv"))
    smr$tissue = t
    smr_list[[t]] = smr
}

ts_ct_anc = read_delim("/path/to/crossTissue_analysis/tissue_celltype_ancestry_ss_new.txt",delim='\t')
smr_list_all = list_rbind(smr_list) %>% left_join(ts_ct_anc,by=c("tissue","cell_type","ancestry")) %>%
 dplyr::filter(!(tissue=="Blood" & cell_type %in% c("T","B","Eryth"))) %>%
  dplyr::filter(sample_size>=50) %>% inner_join(cell_type_cat_colors, by=c("tissue","cell_type")) %>%
    dplyr::mutate(ancestry=factor(ancestry,levels=ancestry_list), 
                ct_anc = paste0(cell_type_short,"_",ancestry), 
                total_count = round(sample_size*mean_count)) %>%
    dplyr::filter(total_count>=10000) 

# --------------------------------------------------------
# 0. Helper functions
# --------------------------------------------------------
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
  ct <- ct_anc_combo$ct
  anc <- ct_anc_combo$anc
  tissue <- ct_anc_combo$tissue
  chrom <- ct_anc_combo$chrom
  
  ct_anc_tmp <- paste0(cell_type_short_dict[[tissue]][ct], "_", anc)
  if (ct_anc_tmp %in% smr_list_all$ct_anc[smr_list_all$tissue == tissue]) {
    file_path <- file.path(tissue_path, "raw", paste0(ct, "_", anc, "_pc5.cis_qtl_pairs.", chrom, ".parquet"))
    
    if (file.exists(file_path)) {
      print(paste0(anc, "_", ct))
      df <- read_parquet(file_path)
      df <- df[df$phenotype_id %in% gene_keep & df$variant_id %in% geno_overlap, ]
      df <- df[df$af >= 0.01 & df$af <= 0.99,]
      df$ancestry <- anc
      df$celltype <- ct
      return(df)
    }
  }
  return(NULL)  # Return NULL if the combination does not meet criteria
}
# --------------------------------------------------------
# 1. Analysis starts
# --------------------------------------------------------
args = commandArgs(TRUE)
tissue = args[1]
chrom = args[2]

mashR_path = "/path/to/crossTissue_analysis/mashR/" 
tissue_path = paste0("/path/to/",tissue,"/sc-eQTL/results/")

if (! file.exists(paste0(mashR_path, "step5_summary_statistics_",tissue,"_chr",chrom,".rds"))){
  geno_overlap = read.table(paste0(mashR_path, "step3_",tissue,"_overlap.bim"))$V1
  top_all = readRDS(paste0(mashR_path, "step1_top_all_", tissue, ".rds"))
  gene_list = readRDS(paste0(mashR_path, "step1_gene_list_", tissue, ".rds"))
  count = table(unlist(gene_list))
  n_max = max(count)
  gene_keep = names(count)[count>=n_max-length(unique(top_all$ancestry))]
  #gene_keep = names(count)[count==n_max]
  rm(top_all)
  print(paste0("Read in raw data for ",tissue," chr",chrom, "."))

  # Create a list of ancestry and cell type combinations
  ct_anc_combinations <- list()
  for (anc in ancestries_all_list[[tissue]]) {
    for (ct in celltypes_all_list[[tissue]]) {
      ct_anc_combinations[[length(ct_anc_combinations) + 1]] <- list(
        anc = anc,
        ct = ct,
        tissue = tissue,
        chrom = chrom
      )
    }
  }

  # Use mclapply to process combinations in parallel
  df_all_list <- mclapply(
    ct_anc_combinations,
    process_combinations,
    mc.cores = 2  # Use all but one core
  )

  # Filter out NULL results
  df_all_list <- Filter(Negate(is.null), df_all_list)

  # Get summary statistics
  df_all_ss = mclapply(df_all_list, GetSS, mc.cores=2)
  df_all_ss = rbindlist(df_all_ss)

  df_all_z_wide <- dcast(df_all_ss[,.(phenotype_id,variant_id,z,ancestry,celltype)], phenotype_id+variant_id~celltype+ancestry, value.var=c("z"))
  df_all_slope_wide <- dcast(df_all_ss[,.(phenotype_id,variant_id,slope,ancestry,celltype)], phenotype_id+variant_id~celltype+ancestry, value.var=c("slope"))
  df_all_se_wide <- dcast(df_all_ss[,.(phenotype_id,variant_id,slope_se,ancestry,celltype)], phenotype_id+variant_id~celltype+ancestry, value.var=c("slope_se"))

  ncond=ncol(df_all_z_wide)
  df_all_pair = df_all_z_wide[,1:2]
  saveRDS(df_all_pair, paste0(mashR_path, "step5_all_pairs_",tissue,"_chr",chrom,".rds"))
  saveRDS(list(df_all_slope_wide=df_all_slope_wide, df_all_se_wide=df_all_se_wide),  paste0(mashR_path, "step5_summary_statistics_",tissue,"_chr",chrom,".rds"))
  
}else{
  summary_stats <- readRDS(paste0(mashR_path, "step5_summary_statistics_",tissue,"_chr",chrom,".rds"))
  df_all_slope_wide <- summary_stats$df_all_slope_wide
  df_all_se_wide <- summary_stats$df_all_se_wide
  ncond=ncol(df_all_slope_wide) 
}

out = readRDS(paste0(mashR_path, "step4_strong_random_for_mash_",tissue,".rds"))
data.temp = mash_set_data(Bhat=as.matrix(out$random_b), Shat=as.matrix(out$random_s))
Vhat = estimate_null_correlation_simple(data.temp)
rm(data.temp)
data_STRONG = mash_set_data(as.matrix(df_all_slope_wide[,3:ncond]), as.matrix(df_all_se_wide[,3:ncond]), V = Vhat)

print("Read in fitted model")
m = readRDS(paste0(mashR_path, "step5_m_",tissue,".rds"))

print("Computed posterior summaries")
m2 = mash(data_STRONG, g = get_fitted_g(m), fixg = TRUE)
saveRDS(m2, paste0(mashR_path, "step5_m2_all_posterior_",tissue,"_chr",chrom,".rds"))

