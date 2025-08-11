##############################################################################
# mashR for assessing the sharing and specificity of sc-eQTL
##############################################################################

library(tidyverse)
library(data.table)
library(ashr)
library(mashr)
library(rmeta)
library(arrow)
library(ComplexHeatmap)
library(UpSetR)
library(parallel)
library(circlize)
library(corrplot)

# 0. Preparing ancestry, tissue, cell type list ------------
# -------------------------------------------------------------------
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

#all_cell_type = read_delim("script_crossTissue/all_cell_type.txt", delim="\t")
ts_ct_anc = read_delim("/path/to/crossTissue_analysis/tissue_celltype_ancestry_ss_new.txt",delim='\t')
smr_list_all = list_rbind(smr_list) %>% left_join(ts_ct_anc,by=c("tissue","cell_type","ancestry")) %>%
  dplyr::filter(sample_size>=50) %>% inner_join(cell_type_cat_colors, by=c("tissue","cell_type")) %>%
    dplyr::mutate(ancestry=factor(ancestry,levels=ancestry_list), 
                ct_anc = paste0(cell_type_short,"_",ancestry), 
                total_count = round(sample_size*mean_count)) %>%
    dplyr::filter(total_count>=10000) 

mashR_path = "/path/to/crossTissue_analysis/mashR/" 

args = commandArgs(TRUE)
tissue = args[1]

tissue_path = paste0("/path/to/",tissue,"/sc-eQTL/results/")
print(paste0("Conducting mashR for ", tissue))

# 1.Find overlapping SNPs across ancestries ------------------------
# ------------------------------------------------------------------
print("Step 0: Find overlapping SNPs across ancestries")
if ( ! file.exists(paste0(mashR_path, "step3_",tissue,"_overlap.bim"))){
geno_list = list()
for (anc in ancestries_all_list[[tissue]]){
    print(anc)
    geno_list[[length(geno_list)+1]]=read.table(paste0(tissue_path, "../genotype/",tissue, "_",anc,"_indivpruned.bim"))$V2
}
geno_all = unlist(geno_list)
if (length(ancestries_all_list[[tissue]])==1){
    geno_overlap = geno_all
    rm(geno_all)
}else{
    geno_all_tab = table(geno_all)
    geno_overlap = names(geno_all_tab)[geno_all_tab==length(ancestries_all_list[[tissue]])]
}
write.table(geno_overlap, paste0(mashR_path, "step3_",tissue,"_overlap.bim"), sep='\n',col=F,row=F,quo=F)
}else{
geno_overlap = read.table(paste0(mashR_path, "step3_",tissue,"_overlap.bim"))$V1
}

# 1. Find the strong set (top sc-eQTL across conditions) ------------
# -------------------------------------------------------------------
# Load top sc-eQTL for each ancestry and cell type
print("Step 1: Find top eQTLs within each condtion")
if (!file.exists(paste0(mashR_path, "step1_top_all_", tissue, ".rds"))){
top_list <- list()
gene_list <- list()
for (anc in ancestries_all_list[[tissue]]){
    for (ct in celltypes_all_list[[tissue]]){
        ct_anc_tmp <- paste0(cell_type_short_dict[[tissue]][ct],"_",anc)
        if (ct_anc_tmp %in% smr_list_all$ct_anc[smr_list_all$tissue == tissue]) {
        if(file.exists(file.path(tissue_path, "sig",paste0(ct,"_",anc,"_pc5_sig_5en8.tsv.gz")))){
        print(paste0(anc,"_",ct))
        df_top <- fread(file.path(tissue_path, "sig", paste0(ct,"_",anc,"_pc5_sig_5en8.tsv.gz"))) # phenotype with all significant eQTLs 
        gene_list[[length(gene_list)+1]] <- fread(file.path(tissue_path, "raw", paste0(ct,"_",anc,"_pc5_perm.tsv")),select="phenotype_id")$phenotype_id
        df_top <- df_top[variant_id %in% geno_overlap,.(phenotype_id, variant_id, slope, slope_se, pval_nominal)] # only select variant in overlaping genotypes
        df_top <- df_top[, .SD[which.min(pval_nominal)], by = phenotype_id] # only keep the top eQTLs
        df_top$ancestry<-anc
        df_top$celltype<-ct
        top_list[[(length(top_list)+1)]] <- df_top
        }
    }
    }
}
top_all = rbindlist(top_list)
saveRDS(top_all, paste0(mashR_path, "step1_top_all_", tissue, ".rds"))
saveRDS(gene_list, paste0(mashR_path, "step1_gene_list_", tissue, ".rds"))
} else {
top_all = readRDS(paste0(mashR_path, "step1_top_all_", tissue, ".rds"))
gene_list = readRDS(paste0(mashR_path, "step1_gene_list_", tissue, ".rds"))
}


# 2.Find overlapping genes across conditions and find the top eQTLs -
# -------------------------------------------------------------------
# 2-1: Across all ancestries
print("Step 2: Find the strong set (top eQTLs across all condtions)")
count = table(unlist(gene_list))
n_max = max(count)
gene_keep = names(count)[count>=n_max-length(unique(top_all$ancestry))]
#gene_keep = names(count)[count==n_max]
if (! file.exists(paste0(mashR_path, "step2_top_all_cond_", tissue, ".rds"))){
top_all_cond = top_all[phenotype_id %in% gene_keep, .SD[which.min(pval_nominal)], by=phenotype_id]
saveRDS(top_all_cond, paste0(mashR_path, "step2_top_all_cond_", tissue, ".rds"))
}else{
	top_all_cond=readRDS(paste0(mashR_path, "step2_top_all_cond_", tissue, ".rds"))
}


# 4. Find the random set  ------------------------------------------
# ------------------------------------------------------------------
if (!file.exists(paste0(mashR_path, "step4_strong_random_for_mash_",tissue,".rds"))){
    print("Step 4: Find the random set")
    # Select random gene-SNP pairs from EUR pDC
    # keep only overlapping variants and only overlapping genes
    process_chromosome <- function(chromosome,ct) {
        print(chromosome)
        df <- read_parquet(file.path(tissue_path, "raw", paste0(ct,"_EUR_pc5.cis_qtl_pairs.", chromosome, ".parquet"))) %>%
            as.data.table()
        # Filter by phenotype_id and variant_id
        df <- df[phenotype_id %in% gene_keep & variant_id %in% geno_overlap,]
        # Function to process each gene for the given chromosome
        process_gene <- function(gene) {
            df_sub <- df[phenotype_id == gene,]
            if (nrow(df_sub)>0){
                df_sub <- df_sub[sample(1:nrow(df_sub), min(50,nrow(df_sub))), .(phenotype_id, variant_id)]
            return(df_sub)}
        }
        # Apply the gene-level processing in parallel
        df_random_list <- lapply(unique(df$phenotype_id), process_gene)
        return(df_random_list)
    }
    if(!file.exists(paste0(mashR_path, "step4_random_list_",tissue,".rds"))){
        # Apply the chromosome processing in parallel
        tmp=table(top_all_cond$celltype)
        ct=names(tmp)[which.max(tmp)]
        df_random_list_all <- mclapply(1:22, function(x) {process_chromosome(x,ct)}, mc.cores=2)
        # Combine all results into a single list
        df_random_list <- rbindlist(unlist(df_random_list_all, recursive = FALSE))
        saveRDS(df_random_list, paste0(mashR_path, "step4_random_list_",tissue,".rds"))
        rm(df_random_list_all)
    }else{
        df_random_list <- readRDS(paste0(mashR_path, "step4_random_list_",tissue,".rds"))
    }

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

    # Select the random gene-SNP pairs from all conditions
    # Select the strong gene-SNP pairs from all conditions
    df_random_all_list = list()
    df_strong_all_list = list()
    print("Step 4: Select the random and strong gene-SNP pairs from all conditions")
    for (anc in ancestries_all_list[[tissue]]){
        for (ct in celltypes_all_list[[tissue]]){
            ct_anc_tmp = paste0(cell_type_short_dict[[tissue]][ct],"_",anc)
            if (ct_anc_tmp %in% smr_list_all$ct_anc[smr_list_all$tissue == tissue]) {
            # Define the file path
            file_path = file.path(tissue_path, "raw", paste0(ct, "_", anc, "_pc5.cis_qtl_pairs.22.parquet"))
            # Check if the parquet file exists
            if (file.exists(file_path)) {
                print(paste0(anc, "_", ct))
                # Read the parquet file once
                df_list = mclapply(1:22, function(i) {
                    print(i)
                    # Read the parquet file
                    df = read_parquet(file.path(tissue_path, "raw", paste0(ct, "_", anc, "_pc5.cis_qtl_pairs.", i, ".parquet"))) %>% as.data.table()

                    # Random selection: merge with df_random_list
                    df_random_select = merge(df, df_random_list, by = c("phenotype_id", "variant_id"))
                    df_random_select = df_random_select[, .(phenotype_id, variant_id, slope, slope_se, pval_nominal)]
                    
                    # Strong selection: merge with top_all_cond
                    df_strong_select = merge(df, top_all_cond[, .(phenotype_id, variant_id)], by = c("phenotype_id", "variant_id"))
                    df_strong_select = df_strong_select[, .(phenotype_id, variant_id, slope, slope_se, pval_nominal)]
                    
                    return(list(random = df_random_select, strong = df_strong_select))
                }, mc.cores = 2)
                
                # Unlist the results from mclapply
                df_list_random = rbindlist(lapply(df_list, function(x) x$random))
                df_list_strong = rbindlist(lapply(df_list, function(x) x$strong))
                
                # Add ancestry and celltype information
                df_list_random$ancestry = anc
                df_list_random$celltype = ct
                df_list_strong$ancestry = anc
                df_list_strong$celltype = ct
                
                # Append the results to the main lists
                df_random_all_list[[length(df_random_all_list) + 1]] = df_list_random
                df_strong_all_list[[length(df_strong_all_list) + 1]] = df_list_strong
            }
            }
        }
    }

    df_random_all_ss = mclapply(df_random_all_list, GetSS, mc.cores=6)
    df_random_all_ss = rbindlist(df_random_all_ss)
    saveRDS(df_random_all_ss, paste0(mashR_path, "step4_random_all_ss_",tissue,".rds"))
    df_strong_all_ss = mclapply(df_strong_all_list, GetSS, mc.cores=6)
    df_strong_all_ss = rbindlist(df_strong_all_ss)
    saveRDS(df_strong_all_ss, paste0(mashR_path, "step4_strong_all_ss_",tissue,".rds"))

    df_random_all_z_wide <- dcast(df_random_all_ss[,.(phenotype_id,variant_id,z,ancestry,celltype)], phenotype_id+variant_id~celltype+ancestry, value.var=c("z"))
    df_random_all_slope_wide <- dcast(df_random_all_ss[,.(phenotype_id,variant_id,slope,ancestry,celltype)], phenotype_id+variant_id~celltype+ancestry, value.var=c("slope"))
    df_random_all_se_wide <- dcast(df_random_all_ss[,.(phenotype_id,variant_id,slope_se,ancestry,celltype)], phenotype_id+variant_id~celltype+ancestry, value.var=c("slope_se"))

    df_strong_all_z_wide <- dcast(df_strong_all_ss[,.(phenotype_id,variant_id,z,ancestry,celltype)], phenotype_id+variant_id~celltype+ancestry, value.var=c("z"))
    df_strong_all_slope_wide <- dcast(df_strong_all_ss[,.(phenotype_id,variant_id,slope,ancestry,celltype)], phenotype_id+variant_id~celltype+ancestry, value.var=c("slope"))
    df_strong_all_se_wide <- dcast(df_strong_all_ss[,.(phenotype_id,variant_id,slope_se,ancestry,celltype)], phenotype_id+variant_id~celltype+ancestry, value.var=c("slope_se"))

    ncond=ncol(df_random_all_z_wide)
    for (i in 3:ncond) {
        # Update NaN for random set
        df_random_all_slope_wide[[i]] = handle_nan_b(df_random_all_slope_wide[[i]])
        df_random_all_se_wide[[i]] = handle_nan_s(df_random_all_se_wide[[i]])
        df_random_all_z_wide[[i]] = handle_nan_z(df_random_all_z_wide[[i]],df_random_all_slope_wide[[i]], df_random_all_se_wide[[i]])
        # Update NaN for strong set
        df_strong_all_slope_wide[[i]] = handle_nan_b(df_strong_all_slope_wide[[i]])
        df_strong_all_se_wide[[i]] = handle_nan_s(df_strong_all_se_wide[[i]])
        df_strong_all_z_wide[[i]] = handle_nan_z(df_strong_all_z_wide[[i]],df_strong_all_slope_wide[[i]], df_strong_all_se_wide[[i]])
        }

    out = list(random_pair=as.data.frame(df_random_all_z_wide)[,1:2],
            random_z=as.data.frame(df_random_all_z_wide)[,3:ncond],
            random_b=as.data.frame(df_random_all_slope_wide[,3:ncond]),
            random_s=as.data.frame(df_random_all_se_wide[,3:ncond]),
            strong_pair=as.data.frame(df_strong_all_z_wide)[,1:2],
            strong_z=as.data.frame(df_strong_all_z_wide)[,3:ncond],
            strong_b=as.data.frame(df_strong_all_slope_wide[,3:ncond]),
            strong_s=as.data.frame(df_strong_all_se_wide[,3:ncond]))

    saveRDS(out,paste0(mashR_path, "step4_strong_random_for_mash_",tissue,".rds"))
}else{
    out=readRDS(paste0(mashR_path, "step4_strong_random_for_mash_",tissue,".rds"))
}


# 5. mashR  -------------------------------------------------------
# ------------------------------------------------------------------
print("Step 5: Start mashR")
## estimate null correlation structure
print("Estimate null correlation structure")
data.temp = mash_set_data(Bhat=as.matrix(out$random_b), Shat=as.matrix(out$random_s))
Vhat = estimate_null_correlation_simple(data.temp)

## transform to mashr format
print("Transform to mashR format")
#data_RANDOM1 = mash_set_data(as.matrix(out$random_b), as.matrix(out$random_s), V = Vhat)
data_RANDOM = mash_update_data(data.temp, V=Vhat)
data_STRONG = mash_set_data(as.matrix(out$strong_b), as.matrix(out$strong_s), V = Vhat)

## set up the data-driven covariance matrix
print("Set up the data-driven covariance matrix")
U.pca = cov_pca(data_STRONG, 5)

## extreme deconvolution
print("Extreme deconvolution")
U.ed = cov_ed(data_STRONG, U.pca)

## set up the canonical covariance matrix
print("Set up the canonical covariance matrix")
U.c = cov_canonical(data_RANDOM)

## covariance matrices
Ulist = c(U.ed, U.c) 


## run mash with both canonical and data-driven cov matrices -- fit the model 
## fit mash to the random tests using both data-driven and canonical covariances -- we have to fit using a random set of tests, and not a dataset that is enriched for strong tests 
## the outputlevel=1 option means that it will not compute posterior summaries for these tests (which saves time)
## - Likelihood calculations took 4504.71 seconds.
## - Model fitting took 13197.51 seconds.
print("Model fitting")
m = mash(data_RANDOM, Ulist = Ulist, outputlevel = 1)
saveRDS(m, paste0(mashR_path, "step5_m_",tissue,".rds"))

## use the fit from the previous run of mash by specifying g=get_fitted_g(m), fixg=TRUE to compute posterior summaries for any subset of tests
##  - Computing 4005 x 2381 likelihood matrix.
## - Likelihood calculations took 513.22 seconds.
## - Computing posterior matrices.
## - Computation allocated took 126.84 seconds.
print("Computed posterior summaries")
m2 = mash(data_STRONG, g = get_fitted_g(m), fixg = TRUE)
saveRDS(m2, paste0(mashR_path, "step5_m2_",tissue,".rds"))

## assess the proportion of significant signals shared by magnitude in each pair of conditions based on posterior means
print("Assess the proportion of significant signals shared by magnitude in each pair of conditions based on posterior means")
m.pairwise_PM <- get_pairwise_sharing(m2, lfsr_thresh = 0.1, factor = 0.5)
saveRDS(m.pairwise_PM, paste0(mashR_path,"step5_pairwise_PM_", tissue,".rds"))


