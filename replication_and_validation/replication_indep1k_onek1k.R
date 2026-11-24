#!/usr/bin/env Rscript
# Compare blood EUR cis-eQTLs from the reanalyzed OneK1K and independent 1K cohorts.
# Usage: Rscript 1k1k_benchmark_indep1k_public.R [project_root] [output_dir]
# Inputs: Blood/sc-eQTL_R/results/{sig,raw} and the cell-type metadata below.

suppressMessages({
    library(tidyverse)
    library(data.table)
    library(arrow)
    library(parallel)
    library(ggvenn)
    library(cowplot)
    library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1L) args[1] else "."
out_dir <- if (length(args) >= 2L) args[2] else file.path(project_root, "crossTissue_analysis_R", "1k1k_indep1k_benchmark")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
results_dir <- file.path(project_root, "Blood", "sc-eQTL_R", "results")
cores <- if (.Platform$OS.type == "windows") 1L else 4L

ct_list <- c('CD4TNC','CD4TEM','Treg','CD8TNC','CD8TRM','CD8TEMRA','MAIT','gdT','NKp','BIN','BMemNS','BMemS','Plasma','MonoC','MonoNC','DC2','pDC')
common_theme <- theme_classic(base_size = 12) +
    theme(axis.text=element_text(color="black"), text=element_text(family="Helvetica"),plot.title=element_text(hjust=0.5))

metadata <- read_delim(file.path(project_root, "crossTissue_analysis_R", "tissue_ancestry_celltype_ss_cc_colors.txt"),
                       delim = "\t", show_col_types = FALSE) %>%
    filter(tissue == "Blood", sample_size >= 50, total_count >= 5000, cell_type %in% ct_list) %>%
    distinct(cell_type, cell_type_short, color)
cell_type_short_dict <- setNames(metadata$cell_type_short, metadata$cell_type)
cell_type_color_dict <- setNames(metadata$color, metadata$cell_type_short)
ct_order <- unname(cell_type_short_dict[ct_list])

# Load results ----------------------
sig_list = list()
smr_list = list()
ts <- "Blood"
for (ct in ct_list){
    for (grp in c("1k1k","sub")){
        anc=paste0("EUR",grp)
        ctxt_file <- file.path(results_dir, "sig", paste0(ct,grp,"_",anc,"_pc5_sig.tsv.gz"))
        if (!file.exists(ctxt_file)) next
        sig <- fread(ctxt_file)
        print(paste(ts,anc,ct))
        if (min(sig$pval_nominal_threshold)<5e-8){
            parquet <- mclapply(1:22, function(chr) {
                f <- read_parquet(file.path(results_dir, "raw", paste0(ct,grp,"_",anc,"_pc5.cis_qtl_pairs.",chr,".parquet")))
                f <- f %>% dplyr::filter(pval_nominal<5e-8)
                return(f)
            }, mc.cores = cores) %>% rbindlist()
        }else{
            parquet <- sig[pval_nominal<5e-8]
        }
        n_egene <- length(unique(parquet$phenotype_id))
        n_eqtl <- nrow(parquet)
        smr_list[[length(smr_list)+1]] <- data.table(group=grp, tissue=ts,cell_type=ct,ancestry=substr(anc,1,3),nPC=5,neGenes=n_egene,neQTLs=n_eqtl)
        parquet$group <- grp
        parquet$tissue <- ts
        parquet$ancestry <- substr(anc,1,3)
        parquet$cell_type <- ct
        sig_list[[length(sig_list)+1]] <- parquet[,.(phenotype_id,variant_id,start_distance, af,ma_samples,ma_count,pval_nominal,slope,slope_se,group,tissue,ancestry,cell_type)]
    }
}

smr_list_ts <- rbindlist(smr_list)
sig_list_all <- rbindlist(sig_list)

smr_list_ts[,cell_type_short := cell_type_short_dict[cell_type]]
sig_list_all[,cell_type_short := cell_type_short_dict[cell_type]]


# Venn plot ---------------------------
p_list <- list()
for (ct in ct_list){
    sig <- sig_list_all[cell_type==ct]
    p_list[[ct]] <- ggvenn(
        list("OneK1K"=unique(sig[group=="1k1k"]$phenotype_id),
            "Indep1K"=unique(sig[group=="sub"]$phenotype_id)
           ),  auto_scale=TRUE,text_size=6,
        set_name_size=0,fill_color=c("#2646af","#b0ccff"),fill_alpha=0.8, show_percentage=FALSE
        ) + labs(title=sig$cell_type_short[1]) + theme(plot.title = element_text(hjust=0.5,size=18))
}
set_names <- c("OneK1K","Indep1K")
legend_plot <- ggplot(data.frame(x=1:2, y=1:2, Group=factor(set_names, levels=c("OneK1K","Indep1K"))), aes(x, y, fill=Group)) +
    geom_tile() +
    scale_fill_manual(values = c(OneK1K="#2646af",Indep1K="#b0ccff")) +
    theme_minimal() +
    theme(legend.position = "right",
          legend.title = element_blank(),
          legend.text = element_text(size = 12))

# Extract only the legend
shared_legend <- get_legend(legend_plot)

p_list[["legend"]] = shared_legend
p <- plot_grid(plotlist = p_list, nrow=3)
ggsave(file.path(out_dir, "eGene_compare.png"), plot = p, width=14,height=7)

# eGene count ----------------------------------------
p <- smr_list_ts %>% dplyr::mutate(cell_type=factor(cell_type_short, levels=ct_order), group=factor(group,levels=c("1k1k","sub"), labels=c("OneK1K","Indep1K"))) %>%
 ggplot(aes(x=cell_type,y=neGenes, group=group, fill=group)) +
    geom_bar(stat="identity", position="dodge", width=0.7) +
    scale_fill_manual(values = c(OneK1K="#2646af",Indep1K="#b0ccff")) +
    common_theme +
    theme(legend.position=c(0.7,0.75),
        axis.text=element_text(size=12,color="black"),
        axis.text.x=element_text(angle=45,hjust=1,vjust=1),
        axis.title=element_text(size=12,color="black")) +
    labs(x="Cell type", y="Number of eGenes", fill="")
ggsave(file.path(out_dir, "eGene_count.png"), plot = p, width=7,height=4)


p <- ggvenn(
        list("OneK1K"=unique(sig_list_all[group=="1k1k"]$phenotype_id),
            "Indep1K"=unique(sig_list_all[group=="sub"]$phenotype_id)
           ),  auto_scale=TRUE,
        set_name_size=0,fill_color=c("#2646af","#b0ccff"),fill_alpha=0.8, show_percentage=TRUE
        ) + theme(plot.title = element_text(hjust=0.5,size=18))

ggsave(file.path(out_dir, "eGene_total_venn.png"), plot = p, width=4,height=3)


# The overlap of sc-eQTLs identified in Indep1K with OneK1K ------------
calcu_std_b_se<-function(z,p,n){
    std_b_hat=z/sqrt(2*p*(1-p)*(n+z^2))
    std_se=1/sqrt(2*p*(1-p)*(n+z^2))
    res<-data.frame(std_b_hat,std_se);
    return(res)
}

sig_list = list()
ts <- "Blood"
for (ct in ct_list){
    print(ct)
    ## Load the sc-eQTLs in Indep1K
    anc_indep1k <- "EURsub"
    ctxt_file_indep1k <- file.path(results_dir, "sig", paste0(ct,"sub_",anc_indep1k ,"_pc5_sig.tsv.gz"))
    df_indep1k <- fread(ctxt_file_indep1k)
    if (min(df_indep1k$pval_nominal_threshold)<5e-8){
        parquet_indep1k <- mclapply(1:22, function(chr) {
            f <- read_parquet(file.path(results_dir, "raw", paste0(ct,"sub_",anc_indep1k,"_pc5.cis_qtl_pairs.",chr,".parquet")))
            f <- f %>% dplyr::filter(pval_nominal<5e-8)
            return(f)
        }, mc.cores = cores) %>% rbindlist()
    }else{
        parquet_indep1k <- df_indep1k[pval_nominal<5e-8]
    }
    rm(df_indep1k)

    anc_1k1k <- "EUR1k1k"
    indep_variants <- unique(parquet_indep1k$variant_id)
    parquet_1k1k <- mclapply(1:22, function(chr) {
        path <- file.path(results_dir, "raw", paste0(ct,"1k1k_",anc_1k1k,"_pc5.cis_qtl_pairs.",chr,".parquet"))
        f <- open_dataset(path, format = "parquet") %>%
            filter(variant_id %in% indep_variants) %>%
            collect() %>% as.data.table()
        merge(f, parquet_indep1k[, .(phenotype_id, variant_id)], by = c("phenotype_id", "variant_id"))
    }, mc.cores = cores) %>% rbindlist()

    parquet <- merge(parquet_indep1k, parquet_1k1k, by=c("phenotype_id","variant_id"),suffixes=c("_indep1k","_1k1k"))

    parquet[, maf_indep1k:=ifelse(af_indep1k>0.5,1-af_indep1k,af_indep1k)]
    parquet[, maf_1k1k:=ifelse(af_1k1k>0.5,1-af_1k1k,af_1k1k)]
    parquet[, c("slope_std", "slope_std_se") := {
        z <- slope_indep1k / slope_se_indep1k
        p <- maf_indep1k
        n <- round(ma_count_indep1k/maf_indep1k/2)
        res <- calcu_std_b_se(z, p, n)
        list(res$std_b_hat, res$std_se)
    }]
    parquet[, c("beta_std","beta_std_se") :={
        z <- slope_1k1k / slope_se_1k1k
        p <- maf_1k1k
        n <- round(ma_count_1k1k/maf_1k1k/2)
        res <- calcu_std_b_se(z, p, n)
        list(res$std_b_hat, res$std_se)
    }]

    parquet$cell_type <- ct
    sig_list[[length(sig_list)+1]] <- parquet
}
saveRDS(sig_list, file.path(out_dir, "eqtl_indep1k_snpgene_1k1k.rds"))

df_celltype_all = rbindlist(sig_list, use.names = TRUE, fill = TRUE) %>% dplyr::mutate(sig=ifelse(pval_nominal_1k1k<5e-8, "Significant","Non-significant"))
df_celltype_all$cell_type_short = factor(cell_type_short_dict[df_celltype_all$cell_type],levels=ct_order)
color_pal = c(cell_type_color_dict[levels(df_celltype_all$cell_type_short)],grey="#bebebe83")

p <- df_celltype_all %>% dplyr::mutate(cell_type_color=ifelse(sig=="Significant",as.character(cell_type_short),"grey")) %>%
 ggplot(aes(x=slope_std, y=beta_std,color=cell_type_color)) +
        geom_point(size=1) + geom_hline(yintercept=0) + geom_vline(xintercept = 0) +
        scale_color_manual(values=color_pal) +
        facet_wrap(~cell_type_short, nrow=3)+
        labs(x="Standardized effect size of sc-eQTLs in Indep1K",
            y="Standardized effect size in OneK1K",
            color="")  +
            theme_classic() + theme(text=element_text(size=12,family="Helvetica"), axis.text=element_text(size=12,color='black'),axis.title=element_text(size=12,color='black'), panel.grid = element_blank(),strip.background=element_rect(colour = "white"),legend.position='none',strip.text = element_text(size=14, face="bold"))+
            guides(color = guide_legend(override.aes = list(size = 5 ), keyheight = unit(1, "cm")  ))
ggsave(file.path(out_dir, "eqtl_indep1k_snpgene_1k1k.png"), plot = p, width=13,height=7)

sign_df = df_celltype_all %>% dplyr::mutate(sign_concordance = (sign(slope_std)==sign(beta_std))) %>%
    group_by(cell_type_short) %>%  summarise(Concordance=mean(sign_concordance, na.rm=TRUE), Disconcordance=1-mean(sign_concordance, na.rm=TRUE), n=length(sign_concordance))

p <- sign_df %>%
  pivot_longer(cols = c(Concordance, Disconcordance),
               names_to = "Concordance",
               values_to = "value") %>%
  mutate(Concordance = factor(Concordance, levels = c("Disconcordance", "Concordance"))) %>%
  ggplot(aes(y = value * 100,
             x = cell_type_short,
             fill = interaction(cell_type_short, Concordance))) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = c(
    # Disconcordance colors (grey for all)
    setNames(rep("#d1d1d1", length(unique(sign_df$cell_type_short))),
             paste0(unique(sign_df$cell_type_short), ".Disconcordance")),
    # Concordance colors (individual colors)
    setNames(cell_type_color_dict[levels(sign_df$cell_type_short)],
             paste0(levels(sign_df$cell_type_short), ".Concordance"))
  )) +
  labs(fill = "", x = "Cell type", y = "Sign concordance (%)") +
  theme_classic() +
  theme(text = element_text(size = 12,family="Helvetica"),legend.position="none",
        axis.text=element_text(size=12,color="black"),
        axis.text.x=element_text(angle=45,hjust=1,vjust=1),
        axis.title=element_text(size=12,color="black"))
ggsave(file.path(out_dir, "eqtl_indep1k_snpgene_1k1k_sign_conc.png"), plot = p, width=7,height=4)

sign_df_shared= df_celltype_all %>% dplyr::filter(pval_nominal_1k1k<5e-8) %>% dplyr::mutate(sign_concordance = (sign(slope_std)==sign(beta_std))) %>%
    group_by(cell_type_short) %>%  summarise(Concordance=mean(sign_concordance, na.rm=TRUE), Disconcordance=1-mean(sign_concordance, na.rm=TRUE), n=length(sign_concordance))

p <- sign_df_shared %>%
  pivot_longer(cols = c(Concordance, Disconcordance),
               names_to = "Concordance",
               values_to = "value") %>%
  mutate(Concordance = factor(Concordance, levels = c("Disconcordance", "Concordance"))) %>%
  ggplot(aes(y = value * 100,
             x = cell_type_short,
             fill = interaction(cell_type_short, Concordance))) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = c(
    # Disconcordance colors (grey for all)
    setNames(rep("#d1d1d1", length(unique(sign_df_shared$cell_type_short))),
             paste0(unique(sign_df_shared$cell_type_short), ".Disconcordance")),
    # Concordance colors (individual colors)
    setNames(cell_type_color_dict[levels(sign_df_shared$cell_type_short)],
             paste0(levels(sign_df_shared$cell_type_short), ".Concordance"))
  )) +
  labs(fill = "", x = "Cell type", y = "Sign concordance (%)") +
  theme_classic() +
  theme(text = element_text(size = 12,family="Helvetica"),legend.position="none",
        axis.text=element_text(size=12,color="black"),
        axis.text.x=element_text(angle=45,hjust=1,vjust=1),
        axis.title=element_text(size=12,color="black"))
ggsave(file.path(out_dir, "eqtl_indep1k_snpgene_1k1k_sign_conc_shared.png"), plot = p, width=7,height=4)
