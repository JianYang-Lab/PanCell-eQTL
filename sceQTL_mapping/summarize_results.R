library(tidyverse)
library(data.table)
library(parallel)

common_theme <- theme_classic(base_size=8) +
    theme(text = element_text(family = "Helvetica"), axis.text=element_text(color="black"))

setwd("/path/to/your/direcotry")
tissue_list = c("Blood","Lung","Skin","Colon","Liver")
ancestry_list = c("EUR", "EAS", "AFR", "AMR")
ancestry_color_list = c(EUR="#66C2A5",EAS="#FC8D62",AFR="#8DA0CB",AMR="#E78AC3")
ancestries_all_list = list(Blood=c("EUR", "EAS", "AFR", "AMR"),
                            Lung=c("EUR","EAS"),
                            Skin=c("EUR"),
                            Colon=c("EUR"),
                            Liver=c("EUR","EAS"))
celltypes_all_list = list(
  Blood = c('CD4TNC','CD4TNCInf','CD4TEM','Treg','CD4TEMRA','CD8TNC','CD8TRM','CD8TEMRA','MAIT','gdT','NKp','NKn','BIN','BMemNS','BMemS','Plasma','Plasmablasts','MonoC','MonoInf','MonoNC','DC1','DC2','Macro','Nph','pDC','Eryth','HSPC','Platelet'),
  
  Lung = c("CD4TNC","CD4TEM","Treg","CD8TEMRA","CD8TEMTRM","CD8TRM","Tfh","NKp","NKn","NKCyc","BIN","BMem","Plasma","MonoC","MonoNC","Macro","MacroCyc","MacroInt","MacroAlv","Mast","DC1","DC2","MigDC","pDC","Eryth","Basal","Suprabasal","Ciliated","Club","Goblet","Ion","AT1","AT2","AeroCap","ECArt","ECCap","ECVenPulm","ECVenSys","LECdiff","LECmat","FibroAdv","FibroAlv","Myofibro","Pericyte","Meso"),
  
  Skin = c("CD4TNC","CD4TEM","Treg","CD8TEMRA","CD8TEMTRM","ILC3","NKp","NKn","BMem","Plasma","MonoMac","Macro1","Macro2","MacroInf","Mast","DC1","DC2","moDC","pDC","LC","MigLC","MigDC","KCundiff","KCdiff","VEC1","VEC2","VEC3","LE1","LE2","Fibro1","Fibro2","Fibro3","Pericyte1","Pericyte2","Melano","Schwann1","Schwann2"),
  
  Liver = c("CD4TNC","CD4TEM","Treg","CD8TEMRA","CD8TRM","MAIT","NKp","NKn","BIN","BMem","Plasma","MonoC","MonoNC","Macro","KC","Mast","Nph","DC1","DC2","HepPP","HepIZ","HepPC","CGC","Endo_portal","LSEC_PP","LSEC_CV","Fibro","HSC"),
  
  Colon = c("CD4TNC","CD4TEM","Tfh","Th1","Th17","Treg","CD8TEMTRM","CD8TRM","gdT","ILC","NKp","NKn","Bf","BGC","BIN","BMem","Plasma","Mono","Macro","MacroIntest","Mast","DC2","Stem","TA","TA_CycSec","Entero","Entero_BEST4","Colono","GobletImm","Goblet","Paneth","Tuft","EEC","ECArt","ECCap","ECVen","LEC","Stromal1","Stromal2","Stromal3","Myofibro","Pericyte","Glia")
)

cell_type_cat_colors = read_delim("crossTissue_analysis_R/tissue_ancestry_celltype_ss_cc_colors.txt", delim='\t') %>% 
  dplyr::filter(sample_size>=50 & total_count >= 10000) %>% 
  dplyr::mutate(lineage = ifelse(lineage %in% c("Erythroid","Megakaryocytic","Neural","Progenitor"),"Other",lineage))

cell_type_short_dict = lapply(names(celltypes_all_list), function(x) {
  ref_order = celltypes_all_list[[x]] 
  tissue_data = cell_type_cat_colors %>% 
    dplyr::filter(tissue == x)
  ordered_cells = tissue_data %>%
    dplyr::filter(cell_type %in% ref_order) %>%
    dplyr::select(cell_type,cell_type_short) %>% distinct() %>%
    dplyr::mutate(cell_type = factor(cell_type, levels = ref_order)) %>%
    dplyr::arrange(cell_type) 
  tmp = ordered_cells$cell_type_short
  names(tmp) = ordered_cells$cell_type
  return(tmp)
})
names(cell_type_short_dict) = names(celltypes_all_list)
cell_type_color_dict = cell_type_cat_colors$color
names(cell_type_color_dict) = cell_type_cat_colors$cell_type_short

lineage_colors <- c(`Lymphoid` = "#8AB9D8", `Myeloid`="#e5a03d", Epithelial = "#8d4c28", Endothelial="#b84e3c", Stromal="#c76b85",Other="#c3e0e5")
celltype_cat_list = names(lineage_colors)

# Get the summary results --------------------------
#-----------------------------
# Nominal pvalue < 5e-8
#-----------------------------
sig_list <- list()
for (ts in tissue_list){
    print(ts)
    smr_list <- list()
    ctxt_files <- list.files(file.path(ts,"sc-eQTL_R","results","sig"),pattern="_sig.tsv.gz",full.names=TRUE)
    ct <- sapply(ctxt_files, function(x) str_split(basename(x),"_")[[1]][1])
    anc <- sapply(ctxt_files, function(x) str_split(basename(x),"_")[[1]][2])
    for (i in 1:length(ctxt_files)){
        sig <- fread(ctxt_files[i])
        ctxt <- str_replace(basename(ctxt_files[i]),"_pc5_sig.tsv.gz","")
        ctxt_split <- str_split(ctxt,"_")[[1]]
        ct <- paste(ctxt_split[-length(ctxt_split)],collapse="_")
        anc <- ctxt_split[length(ctxt_split)]
        print(paste(ts,anc,ct))
        if (min(sig$pval_nominal_threshold)<5e-8){
        parquet <- mclapply(1:22, function(chr) {
            f <- read_parquet(paste0(ts,"/sc-eQTL_R/results/raw/",ct,"_",anc,"_pc5.cis_qtl_pairs.",chr,".parquet"))
            f <- f %>% dplyr::filter(pval_nominal<5e-8)
            return(f)
        }, mc.cores=1) %>% rbindlist()
        }else{
            parquet <- sig[pval_nominal<5e-8]
        }
        n_egene <- length(unique(parquet$phenotype_id))
        n_eqtl <- nrow(parquet)
        smr_list[[length(smr_list)+1]] <- data.table(tissue=ts,cell_type=ct,ancestry=anc,nPC=5,neGenes=n_egene,neQTLs=n_eqtl)
        parquet$tissue <- ts; parquet$ancestry <- anc; parquet$cell_type=ct
        sig_list[[length(sig_list)+1]] <- parquet
    }
    smr_list_ts <- rbindlist(smr_list) 
    write.table(smr_list_ts, file.path(ts,"sc-eQTL_R","results","sig","eQTL_summary_5en8.tsv"), sep='\t', row=F,quo=F,col=T)
}
saveRDS(sig_list, "crossTissue_analysis_R/summary/sig_list.rds")

smr_list <- list()
for (t in tissue_list){
    smr <- read_table(file.path(t,"sc-eQTL_R","results","sig","eQTL_summary_5en8.tsv")) %>% distinct()
    smr$tissue <- t
    smr_list[[t]] <- smr
}

smr_list_all <- list_rbind(smr_list) %>% 
inner_join(cell_type_cat_colors,by=c("tissue","cell_type","ancestry")) %>% 
    dplyr::mutate(ancestry=factor(ancestry,levels=ancestry_list), 
                ct_anc = paste0(cell_type_short,"-",ancestry)) %>% distinct()

#-----------------------------
#Conditionally independent analysis
#-----------------------------
ind_res <- list()
for (i in 1:nrow(smr_list_all)){
  anc=smr_list_all$ancestry[i]
  ct=smr_list_all$cell_type[i]
  tissue=smr_list_all$tissue[i]
  
  if(file.exists(file.path(tissue,"sc-eQTL_R","results","raw",paste0(ct,"_",anc,"_pc5_indep.tsv")))){
    print(paste0(tissue,"_", ct, "_", anc))
    dt = fread(file.path(tissue,"sc-eQTL_R","results","raw",paste0(ct,"_",anc,"_pc5_indep.tsv")))
    df <- dt[, .(count = .N), by = phenotype_id]
    res <- data.frame(table(df$count))
    res$ancestry <- anc
    res$cell_type <- ct 
    res$tissue <- tissue
    ind_res[[length(ind_res)+1]] <- res 
  }
}

ind_res_all1 <- rbindlist(ind_res) %>% 
    left_join(cell_type_cat_colors, by=c("tissue","ancestry","cell_type")) %>%
    dplyr::mutate(ancestry=factor(ancestry,levels=ancestry_list), 
                ct_anc = paste0(cell_type_short,"-",ancestry))
saveRDS(ind_res_all1, "crossTissue_analysis_R/summary/ind_res_all1.rds")

eqtl_rank_count_plot1 <- function(data, t){
  df = data %>%dplyr::filter(tissue==t)
  df$cell_type_short = factor(df$cell_type_short, levels=cell_type_short_dict[[t]][celltypes_all_list[[t]]])
  df <- df %>% dplyr::mutate(
    ancestry=factor(ancestry,levels=ancestry_list),
    rank=factor(ifelse(as.numeric(Var1)>=5,">=5",Var1), levels=c(">=5",4,3,2,1))) %>% 
  group_by(rank,ancestry,cell_type,ct_anc,cell_type_short,tissue) %>% summarise(count=sum(Freq)) %>%
  arrange(cell_type_short,ancestry) 
  df <- df %>% dplyr::mutate(ct_anc = factor(ct_anc,levels=unique(df$ct_anc)))

  alpha_level = c(0.2,0.4,0.6,0.8,1) 
  names(alpha_level) = levels(df$rank) 
  p = df %>% ggplot(aes(fill=cell_type_short)) + 
    geom_bar(aes(x=ct_anc,y=count/1000, alpha=rank),stat="identity",position="stack")+
    facet_grid(cols=vars(cell_type_short),scales="free", space="free") + 
    theme_classic() +
    scale_fill_manual(values=cell_type_color_dict[levels(df$cell_type_short)])+
    scale_alpha_manual(values=alpha_level) + 
    scale_y_continuous(labels=scales::label_number(accuracy=0.1))+
    theme(text=element_text(size=9,family="Helvetica"),
      axis.title=element_text(size=9,color="black"), 
      axis.text=element_text(size=9,color="black"), 
      axis.text.x=element_text(size=8.5,angle=90,hjust=1,vjust=0.5,color="black"),
      strip.background.x = element_blank(), strip.text.x = element_blank(),
      axis.title.x=element_blank(),plot.title=element_text(size=12,hjust=0.5, face="bold"),
      plot.margin=margin(10,10,10,10),
      legend.key.size=grid::unit(2.5, "mm"),legend.position=c(0.35,0.78))+
    labs(x="", y="Number of eGenes (k)", alpha="# ind. QTLs", fill="Cell type", title=t) + 
    guides(fill = "none")
  return(p)
}

p_list = lapply(tissue_list, function(t) eqtl_rank_count_plot1(ind_res_all1,t))
plot_grid(
  p_list[[1]]+theme(legend.position=c(0.2,0.75)), 
  plot_grid(p_list[[2]]+theme(legend.position=c(0.3,0.75)), p_list[[3]]+theme(axis.title.y=element_blank()), nrow = 1, axis="l",align="v",rel_widths = c(1,0.5)),
  plot_grid(p_list[[5]]+theme(legend.position=c(0.2,0.75)), p_list[[4]]+theme(axis.title.y=element_blank(),legend.position=c(0.38,0.75)), nrow = 1, axis="l",align="v",rel_widths = c(1,0.55)),
  nrow = 3, 
  rel_heights= c(1, 1, 1)  
)
ggsave("crossTissue_analysis_R/summary/independent_eqtl_count_rank_5en8_horizontal.png", height=8,width=11) 