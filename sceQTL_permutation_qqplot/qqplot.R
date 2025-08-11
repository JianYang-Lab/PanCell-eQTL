library(tidyverse)
library(data.table)
library(cowplot)
library(RColorBrewer)
library(ggtext)
library(ggnewscale)
library(latex2exp)
library(arrow)
library(parallel)


tissue_list = c("Blood","Lung","Skin","Colon","Liver")
ancestry_list = c("EUR", "EAS", "AFR", "AMR")
ancestry_color_list = c(EUR="#66C2A5",EAS="#FC8D62",AFR="#8DA0CB",AMR="#E78AC3")

# Color palettes for cell types -------------------------
celltypes_all_list = list(Blood=c('CD4TNC','CD4TEM','Treg','CD8TNC','CD8TEM','CD8TEMRA','MAIT','NKp','NKn','BIN','BMem','Plasma','Plasmablasts','MonoC','MonoNC','Nph','DC','pDC'),
Liver = c('CD4T',"CD8T","Treg","MAIT","Th","gdT",'Circulating_NK','Resident_NK','B','Plasma','Monocytes','Macrophages','DC','pDC','Neutrophils','Basophils','Hepatocytes','Cholangiocytes','Endothelium','Fibroblasts'),
Skin = c("Th","Tc", "Treg","NK","DC1","DC2","MigDC","Macro1","Macro2","MonoMac","Mast","KCdiff","KCundiff","Melanocyte","VE1","VE2","VE3","LE1","Pericyte1","Pericyte2","F2"),
Lung = c("CD4T","CD8T","NK","B","Monocytes","Mast","Macrophages","DC","AT1","AT2","Airway_Epi_Multiciliated","Airway_Epi_Basal","Airway_Epi_Secretory","EC_arterial","EC_capillary","EC_venous","LEC_mature","LEC_diff","Fibroblasts","Smooth_muscle","Mesothelium","Rare"),
Colon = c("CD4T","CD8T","Treg","Th","gdT","NKT","NK","ILC3","BIN","BMem","Bcyc","Plasma","Mono","Macro","cDC2","Mast","Colonocyte","GLoblet","Tuft","TA", "EEC" ,"ECcap","ECven","Stromal1","Stromal2","Myofibroblast","Glia"))


cell_type_cat_colors = read_delim("crossTissue_analysis/cell_type_colors_new.txt", delim='\t') 

cell_type_short_dict = lapply(tissue_list,function(x) {
  tmp=cell_type_cat_colors$cell_type_short[cell_type_cat_colors$tissue==x]
  names(tmp) = cell_type_cat_colors$cell_type[cell_type_cat_colors$tissue==x]
  tmp 
  })
names(cell_type_short_dict) = tissue_list
cell_type_color_dict = cell_type_cat_colors$cell_type_colors
names(cell_type_color_dict) = cell_type_cat_colors$cell_type_short


# Load smr_list_all -----------------------
smr_list = list()
for (t in tissue_list){
    smr = read_table(file.path(t,"sc-eQTL","results","sig","eQTL_summary_5en8_maf01.tsv")) %>% distinct()
    smr$tissue = t
    smr_list[[t]] = smr
}

ts_ct_anc = read_delim("crossTissue_analysis/tissue_celltype_ancestry_ss_new.txt",delim='\t') %>% distinct()
smr_list_all = list_rbind(smr_list) %>% 
left_join(ts_ct_anc,by=c("tissue","cell_type","ancestry")) %>% 
  dplyr::filter(sample_size>=50) %>% inner_join(cell_type_cat_colors, by=c("tissue","cell_type")) %>%
    dplyr::mutate(ancestry=factor(ancestry,levels=ancestry_list), 
                ct_anc = paste0(cell_type_short,"-",ancestry),
                total_count = round(sample_size*mean_count)) %>% distinct()

smr_list_all = smr_list_all %>% dplyr::rename(neGenes=neGenes_5en8, neQTLs=neQTLs_5en8)
smr_list_all = smr_list_all %>% dplyr::filter(!(tissue=="Skin" & ancestry=="EAS"))

# Load pvalues from raw data -------------
args=commandArgs(TRUE)
Tissue=args[1]
Ancestry=args[2]

if (!file.exists(paste0("crossTissue_analysis/qqplot/quantile_data_",Tissue,"_",Ancestry,".rds"))){
  # Load pvalues
  if (!file.exists(paste0("crossTissue_analysis/qqplot/pval_list_",Tissue,"_",Ancestry,".rds"))){
  pval_list_all <- list()
  for (i in 1:nrow(smr_list_all)){
    anc=smr_list_all$ancestry[i]
    ct=smr_list_all$cell_type[i]
    tissue=smr_list_all$tissue[i]
    ss=smr_list_all$sample_size[i]
    if (anc==Ancestry & tissue==Tissue){
      if(file.exists(paste0(tissue,"/sc-eQTL/results/raw/",ct,"_",anc,"_pc5.cis_qtl_pairs.22.parquet"))){
    print(paste0(tissue,"_", ct, "_", anc))
    pval_list = list()
    pval_list <- mclapply(1:22, function(chr) {
      cat(sprintf("Chr%d processed\n", chr))
      parquet_file <- paste0(tissue,"/sc-eQTL/results/raw/",ct,"_",anc,"_pc5.cis_qtl_pairs.",chr,".parquet")
      pval=read_parquet(parquet_file, col_select=c("af","pval_nominal")) %>% dplyr::filter(af>=0.01 & af<=0.99) %>% pull(pval_nominal)
      return(pval)
    }, mc.cores=2) 
    pval_list_all[[paste0(tissue,"-",ct,"-",anc)]] = unlist(pval_list)
      } 
    }
  }
  saveRDS(pval_list_all, paste0("crossTissue_analysis/qqplot/pval_list_",Tissue,"_",Ancestry,".rds"))
  } else {
    pval_list_all = readRDS(paste0("crossTissue_analysis/qqplot/pval_list_",Tissue,"_",Ancestry,".rds"))
  }

  # Load permutated results
  anc=Ancestry
  ## permuted control cell type
  ct=ifelse(Tissue=="Blood","CD4TNC","CD4T")
  tissue=Tissue
  if(file.exists(paste0(tissue,"/sc-eQTL/results/perm/",ct,"_",anc,"_pc5_perm2.cis_qtl_pairs.22.parquet"))){
      print(paste0(tissue,"_", ct, "_", anc))
      pval_list <- mclapply(1:22, function(chr) {
        cat(sprintf("Chr%d processed\n", chr))
        parquet_file <- paste0(tissue,"/sc-eQTL/results/perm/",ct,"_",anc,"_pc5_perm2.cis_qtl_pairs.",chr,".parquet")
        pval=read_parquet(parquet_file, col_select=c("af","pval_nominal")) %>% dplyr::filter(af>=0.01 & af<=0.99) %>% pull(pval_nominal)
        return(pval)
    }, mc.cores=2) 
      pval_list_all[[paste0(tissue,"-",ct,"-",anc,"-perm")]] = unlist(pval_list)
  } 

  # Create qqplot statistics
  pval_list_names = names(pval_list_all)

  pval_all= lapply(pval_list_names, function(x) {
    print(x)
    pvals <- pval_list_all[[x]]
    # Remove NA and invalid p-values
    pvals <- pvals[!is.na(pvals) & pvals > 0 & pvals <= 1]
    expected = -log10(ppoints(length(pvals)))
    observed =  -log10(sort(pvals))
    group = ifelse(str_detect(x,"perm"),paste0(ct,"-perm"),cell_type_short_dict[[Tissue]][str_split(x,"-")[[1]][2]])
    data.frame(expected=expected,observed=observed, group=group)
    }
)

  rm(pval_list_all)
  gc()

  pval_all = rbindlist(pval_all)
  gc()

  # Define the number of quantiles 
  num_quantiles <- 1000000
  probabilities <- seq(0, 1, length.out = num_quantiles)

  # Calculate quantiles for each group
  quantile_data <- pval_all[, .(
    expected_quantile = quantile(expected, probs = probabilities, na.rm = TRUE),
    observed_quantile = quantile(observed, probs = probabilities, na.rm = TRUE)
  ), by = group]

  rm(pval_all); gc()
  saveRDS(quantile_data, paste0("crossTissue_analysis/qqplot/quantile_data_",Tissue,"_",Ancestry,".rds"))
  
}else{
  quantile_data = readRDS(paste0("crossTissue_analysis/qqplot/quantile_data_",Tissue,"_",Ancestry,".rds"))
}


# quantile QQ plot ------------------------------
quantile_data = quantile_data %>% 
dplyr::mutate(group = factor(group,levels=c(cell_type_short_dict[[Tissue]][celltypes_all_list[[Tissue]]],ifelse(Tissue=="Blood","CD4TNC-perm","CD4T-perm"))) )

if(Tissue=="Blood"){
  col_pal=c(cell_type_color_dict[levels(quantile_data$group)],`CD4TNC-perm`="#7b17b8")
}else{
  col_pal=c(cell_type_color_dict[levels(quantile_data$group)],`CD4T-perm`="#7b17b8")
}

quantile_data %>%
ggplot(aes(x = expected_quantile, y = observed_quantile, color=group)) +
  geom_point(size=1.5) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(title = paste0(Tissue,"-",Ancestry),
       x = expression(Expected~~-log[10](italic(p))),
       y = expression(Observed~~-log[10](italic(p))),
       color = "") + 
  scale_color_manual(values=col_pal) + theme_classic() +
  theme(legend.position="inside",
        legend.title=element_blank(),
        text=element_text(size=24,family="Helvetica"), 
        plot.title=element_text(hjust=0.5,size=28), 
        legend.position.inside=c(0.25,0.75),
        legend.text = element_text(size=12)) + 
  guides(color = guide_legend(ncol = ifelse(Tissue=="Blood" & Ancestry=="AMR",1,2)))

ggsave(paste0("crossTissue_analysis/qqplot/qqplot_",Tissue,"_",Ancestry,"_1mq_with_perm2.png"), height=7, width=7)
