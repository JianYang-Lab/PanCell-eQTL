# Create eSNP bed files for different tissues and ancestries
# This script processes significant eSNP files, extracts lead SNPs, and saves them in
# BED format for functional annotation.

# Load necessary libraries
library(data.table)
library(tidyverse)


# Create eSNP bed
for (tissue in c("Blood","Lung","Liver","Skin","Colon")){
    for (ancestry in c("EUR","EAS","AMR","AFR")){
    print(tissue)
    sig_path = paste0("/path/to/",tissue,"/sc-eQTL/results/sig")

    files=list.files(sig_path,pattern=paste0(ancestry,"_pc5_sig_5en8.tsv.gz"))
    

    if (length(files)>0){
        print(ancestry)
    
        df_all = lapply(files, function(f) {
            df=fread(file.path(sig_path,f))
            df[, c("chr","pos") := tstrsplit(variant_id, "_", keep=1:2)]
            df[, chr := paste0("chr",chr)]
            df[, start := as.numeric(pos)-1 ]
            df = df[,.(chr,start,pos,variant_id)]
            df
        }) %>% rbindlist()

        df_all = distinct(df_all)
        if (!dir.exists(paste0("/path/to/",tissue,"/sc-eQTL/results/functional_annotation/chromHMM/"))){
            dir.create( paste0("/path/to/",tissue,"/sc-eQTL/results/functional_annotation/chromHMM/"))
        }
        write.table(df_all, paste0("/path/to/",tissue,"/sc-eQTL/results/functional_annotation/chromHMM/",tissue,"_",ancestry,"_all_esnp_5en8.bed"),col=F,row=F,sep='\t',quo=F)
    }
    }
}