################################################################################
#
# Script: Generation of Main Figures 6-7 (GWAS Integration)
#
# Description:
# This script contains the R code to generate figures related to the
# integration of single-cell eQTLs with Genome-Wide Association Study (GWAS)
# summary statistics using SMR (Summary-data-based Mendelian Randomization)
# and COLOC (colocalization analysis).
#
# - Figure 6: Focuses on GWAS integration for Blood (EAS) sc-eQTLs.
# - Figure 7: Focuses on GWAS integration for Lung (EUR) sc-eQTLs.
#
# The script generates complex visualizations including scatterpie heatmaps,
# locus plots, and effect plots.

################################################################################


## 1. SETUP & GLOBAL PARAMETERS
# ------------------------------------------------------------------------------
print("--- Section 1: Setup & Global Parameters ---")

# --- 1.1: Load Libraries ---
library(tidyverse)
library(data.table)
library(arrow)
library(ggtext)
library(ggh4x)
library(cowplot)
library(ggvenn)
library(scatterpie)
library(EnsDb.Hsapiens.v86)
library(biomaRt)
library(locuszoomr)

# --- 1.2: User-Defined Paths and Parameters ---

# Set the base directory for the project
base_dir <- "/path/to/pancell_eqtl/" 

# Define output directory for figures and intermediate data
output_dir <- file.path(base_dir, "final_figures_and_data/")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Create a temporary directory for SMR/locus plot intermediate files
tmp_dir <- file.path(base_dir, "tmp_smr_files/")
dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)

# Other paths
smr_coloc_results_dir <- file.path(base_dir, "smr_coloc_results/")
metadata_dir <- file.path(base_dir, "metadata/")
gwas_lists_dir <- file.path(base_dir, "gwas_lists/") # Centralize GWAS list paths
ref_dir <- file.path(base_dir, "reference_data/") # Path to reference genotypes (e.g., 1000G)
resource_dir <- file.path(base_dir, "resources/") # Path to additional resources

# Path to the SMR executable
smr_executable <- "/path/to/smr-1.3.1-linux-x86_64/smr"

# --- 1.4: Load and Process Global Metadata ---
tissue_list <- c("Blood", "Lung", "Skin", "Colon", "Liver")
ancestry_list <- c("EUR", "EAS", "AFR", "AMR")
ancestry_color_list <- c(EUR = "#66C2A5", EAS = "#FC8D62", AFR = "#8DA0CB", AMR = "#E78AC3")

ancestries_all_list <- list(
  Blood = c("EUR", "EAS", "AFR", "AMR"), Lung = c("EUR", "EAS"), Skin = c("EUR"),
  Colon = c("EUR"), Liver = c("EUR", "EAS")
)

celltypes_all_list <- list(
  Blood = c('CD4TNC','CD4TEM','Treg','CD8TNC','CD8TEM','CD8TEMRA','MAIT','NKp','NKn','BIN','BMem','Plasma','Plasmablasts','MonoC','MonoNC','Nph','DC','pDC'),
  Liver = c('CD4T',"CD8T","Treg","MAIT","Th","gdT",'Circulating_NK','Resident_NK','B','Plasma','Monocytes','Macrophages','DC','pDC','Neutrophils','Basophils','Hepatocytes','Cholangiocytes','Endothelium','Fibroblasts'),
  Skin = c("Th","Tc", "Treg","NK","DC1","DC2","MigDC","Macro1","Macro2","MonoMac","Mast","KCdiff","KCundiff","Melanocyte","VE1","VE2","VE3","LE1","Pericyte1","Pericyte2","F2"),
  Lung = c("CD4T","CD8T","NK","B","Monocytes","Mast","Macrophages","DC","AT1","AT2","Airway_Epi_Multiciliated","Airway_Epi_Basal","Airway_Epi_Secretory","EC_arterial","EC_capillary","EC_venous","LEC_mature","LEC_diff","Fibroblasts","Smooth_muscle","Mesothelium","Rare"),
  Colon = c("CD4T","CD8T","Treg","Th","gdT","NKT","NK","ILC3","BIN","BMem","Bcyc","Plasma","Mono","Macro","cDC2","Mast","Colonocyte","GLoblet","Tuft","TA", "EEC" ,"ECcap","ECven","Stromal1","Stromal2","Myofibroblast","Glia")
)

# Load cell type color and category information
cell_type_cat_colors <- read_delim(file.path(metadata_dir,"cell_type_colors_new.txt"), delim = '\t') %>%
  mutate(lineage = case_when(
    cell_type_category == "T cells" ~ "Immune Lymphoid",
    cell_type_category == "B cells" ~ "Immune Lymphoid",
    cell_type_category == "ILC" ~ "Immune Lymphoid",
    cell_type_category == "Myeloid" ~ "Immune Myeloid",
    cell_type_category == "Endothelial" ~ "Endothelial",
    cell_type_category == "Epithelial" ~ "Epithelial",
    cell_type_category == "Fibroblasts/muscle cells" ~ "Mesenchymal",
    cell_type == "Mesothelium" ~ "Mesenchymal",
    .default = "Other"
  ))

# Create dictionaries for cell type names and colors
cell_type_short_dict <- lapply(tissue_list, function(x) {
  tmp <- cell_type_cat_colors$cell_type_short[cell_type_cat_colors$tissue == x]
  names(tmp) <- cell_type_cat_colors$cell_type[cell_type_cat_colors$tissue == x]
  tmp[celltypes_all_list[[x]]]
})
names(cell_type_short_dict) <- tissue_list
cell_type_color_dict <- cell_type_cat_colors$cell_type_colors
names(cell_type_color_dict) <- cell_type_cat_colors$cell_type_short
lineage_colors <- c(`Immune Lymphoid` = "#8AB9D8", `Immune Myeloid` = "#e5a03d", Epithelial = "#8d4c28", Endothelial = "#b84e3c", Mesenchymal = "#c76b85", Other = "#c3e0e5")

# Load and process summary statistics to define high-quality conditions
smr_list <- list()
for (t in tissue_list) {
  smr_file <- file.path(t, "sc-eQTL/results/sig/eQTL_summary_5en8_maf01.tsv")
  if(file.exists(smr_file)){
    smr <- read_table(smr_file) %>% distinct()
    smr$tissue <- t
    smr_list[[t]] <- smr
  }
}

ts_ct_anc <- read_delim(file.path(metadata_dir,"tissue_celltype_ancestry_ss_new.txt"), delim = '\t') %>% distinct()

smr_list_all <- list_rbind(smr_list) %>%
  left_join(ts_ct_anc, by = c("tissue", "cell_type", "ancestry")) %>%
  filter(sample_size >= 50) %>%
  inner_join(cell_type_cat_colors, by = c("tissue", "cell_type")) %>%
  mutate(
    ancestry = factor(ancestry, levels = ancestry_list),
    ct_anc = paste0(cell_type_short, "-", ancestry),
    total_count = round(sample_size * mean_count)
  ) %>%
  distinct() %>%
  filter(total_count > 10000) %>%
  rename(neGenes = neGenes_5en8, neQTLs = neQTLs_5en8)

gene_bed_path <- file.path(metadata_dir, "genes.bed")
gene_bed <- fread(gene_bed_path)

## 2. HELPER FUNCTIONS
# ------------------------------------------------------------------------------
print("--- Section 2: Defining Helper Functions ---")

# Helper function to query and process eQTL data
query_eqtl_data <- function(prefix, probe, besd_dir, tmp_dir, smr_executable) {
    output_file <- file.path(tmp_dir, paste0(prefix, "_", probe, ".txt"))
    if (!file.exists(output_file)) {
        command <- paste(
            shQuote(smr_executable),
            "--beqtl-summary", shQuote(file.path(besd_dir, prefix)),
            "--query 1 --probe", probe,
            "--out", shQuote(file.path(tmp_dir, paste0(prefix, "_", probe)))
        )
        system(command)
    }
    query_data <- fread(output_file)
    return(query_data)
}

# Helper function for loading GWAS SNPs and get positions
load_gwas_snp <- function(GWAS_list, trait) {
    gwas_path <- tail(GWAS_list[GWAS_list$X1 == trait, ]$X4, n = 1)
    snp_ids <- fread(gwas_path, select = "SNP", header = TRUE)$SNP
    return(snp_ids)
}

# Helper function for selecting GWAS query file
select_gwas_snp <- function(gene, snp_ids, GWAS_list, trait, gene_bed, tmp_dir) {
    output_file <- file.path(tmp_dir, paste0(trait, "_", gene, "_sub.txt"))
    if (!file.exists(output_file)) {
        gene_info <- gene_bed[gene_bed$gene_name == gene, ]
        chromosome <- gene_info$chr[1]
        tss_position <- gene_info$start[1]
        region_start <- max(tss_position - 1e6, 1)
        region_end <- tss_position + 1e6
        
        # NOTE: Assumes dbSNP files are in a specific path. Centralize this path.
        dbsnp_path <- file.path(ref_dir, "dbSNP_b151/split_into_chrAuto/")
        target_snps_pos_file <- file.path(tmp_dir, "target_snps_pos.txt")

        # Extract SNPs in the region of interest to save processing time and memory
        system(sprintf("awk '$2>=%s && $2<=%s {print $1,$2,$3}' %s > %s", region_start, region_end, shQuote(file.path(dbsnp_path, paste0("chr", chromosome, "_variants_hg38_b151.txt"))), shQuote(target_snps_pos_file)))
        
        target_snps <- fread(target_snps_pos_file, col.names = c("chr", "pos", "SNP"))
        temp_snp_file <- tempfile(tmpdir = tmp_dir)
        write.table(target_snps$SNP, temp_snp_file, row.names = FALSE, col.names = FALSE, quote = FALSE)
        
        gwas_path <- GWAS_list[GWAS_list$X1 == trait, ]$X4
        gwas_columns <- c("SNP", "A1", "A2", "freq", "b", "se", "p", "N")
        gwas_target <- fread(
            cmd = sprintf("grep -F -w -f %s %s", shQuote(temp_snp_file), shQuote(gwas_path)),
            col.names = gwas_columns
        )
        
        gwas_target <- merge(gwas_target, target_snps, by = "SNP", all.x = TRUE)
        write.table(gwas_target, output_file, quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
        
        unlink(c(target_snps_pos_file, temp_snp_file)) # Clean up intermediate files
    } else {
        gwas_target <- fread(output_file)
    }
    return(gwas_target)
}


# Helper function for reading SMR data from a plot file
ReadSMRData = function(plotfile)
{
    SMRData = list();
    key=c("$probe","$SNP","$GWAS","$eQTL");
    skiplines=0;
    keywords=scan(plotfile, what="", nlines=1, skip=skiplines);
    if(keywords[1]!=key[1])
    {
        print("ERROR: plot file is not correct!");
        quit();
    }
    nprobes=as.numeric(keywords[2]);
    SMRData$probeID=keywords[3];
  
    
    skiplines=skiplines+1;
    SMRData$SMR=read.table(plotfile, header=F, nrows=nprobes, skip=skiplines);
    skiplines=skiplines+nprobes;
    keywords=scan(plotfile, what="", nlines=1, skip=skiplines);
    if(keywords[1]!=key[2])
    {
        print("ERROR: plot file is not correct!");
        quit();
    }
    nrs=as.numeric(keywords[2]);
    skiplines=skiplines+1;
    SMRData$SNP=read.table(plotfile, header=F, nrows=nrs, skip=skiplines);
    skiplines=skiplines+nrs;
    keywords=scan(plotfile, what="", nlines=1, skip=skiplines);
    if(keywords[1]!=key[3])
    {
        print("ERROR: plot file is not correct!");
        quit();
    }
    ngwas=as.numeric(keywords[2]);
    skiplines=skiplines+1;
    SMRData$GWAS=read.table(plotfile, header=F, nrows=ngwas, skip=skiplines);
    skiplines=skiplines+ngwas;
    keywords=scan(plotfile, what="", nlines=1, skip=skiplines);
    if(keywords[1]!=key[4])
    {
        print("ERROR: plot file is not correct!");
        quit();
    }
    neqtl=as.numeric(keywords[2]);
    skiplines=skiplines+1;
    
    keywords=scan(plotfile, what="", nlines=1, skip=skiplines);
    prbname=keywords[1];
    neqtlsnp=as.numeric(keywords[2]);
    skiplines=skiplines+1;
    SMRData$eQTL=read.table(plotfile, header=F, nrows=neqtlsnp, skip=skiplines);
    SMRData$eQTL=cbind(prbname,SMRData$eQTL)
    skiplines=skiplines+neqtlsnp;
    if(neqtl>1)
    {
        for(i in 2:neqtl)
        {
            keywords=scan(plotfile, what="", nlines=1, skip=skiplines);
            prbname=keywords[1];
            neqtlsnp=as.numeric(keywords[2]);
            skiplines=skiplines+1;
            raweQTLtmp=read.table(plotfile, header=F, nrows=neqtlsnp, skip=skiplines);
            raweQTLtmp=cbind(prbname,raweQTLtmp);
            SMRData$eQTL=rbind(SMRData$eQTL,raweQTLtmp);
            skiplines=skiplines+neqtlsnp;
        }
    }
    
    keywords=scan(plotfile, what="", nlines=1, skip=skiplines);
    if(length(keywords)>0)
    {
        if(keywords[1]!="$Gene")
        {
            print("ERROR: plot file is not correct!");
            quit();
        }
        ngenes=as.numeric(keywords[2]);
        skiplines=skiplines+1;
        SMRData$Gene=read.table(plotfile, header=F, nrows=ngenes, skip=skiplines);
    }
    return(SMRData)
}

# Helper function for plotting SMR effect sizes
SMREffectggPlot <- function(data = SMRData, title = "", cisWindow = 2000, transWindow = 5000, pointsize = 10, r2_threshold=0.2) {
  # Extract and filter data
  snbuf <- which(as.character(data$eQTL[, 1]) == data$probeID)
  if (length(snbuf) == 0) {
    stop(paste("ERROR: no eQTL information found for probe", data$probeID, "!"))
  }
  
  plotData <- data$eQTL[snbuf, ]
  plotData <- plotData[!is.na(plotData[, 5]), ]
  
  snpbuf <- Reduce(intersect, list(as.character(plotData[, 2]), data$GWAS[, 1]))
  plotData <- plotData[match(snpbuf, as.character(plotData[, 2])), ]
  plotGWAS <- data$GWAS[match(snpbuf, as.character(data$GWAS[, 1])), ]

  # Compute effect sizes
  plotData <- transform(
    plotData,
    bZX = as.numeric(as.character(plotData[, 3])),
    seZX = as.numeric(as.character(plotData[, 4])),
    snpCorr = as.numeric(as.character(plotData[, 5]))
  )
  plotGWAS <- transform(
    plotGWAS,
    bZY = as.numeric(as.character(plotGWAS[, 2])),
    seZY = as.numeric(as.character(plotGWAS[, 3]))
  )

  plotData$bZY <- plotGWAS$bZY
  plotData$seZY <- plotGWAS$seZY

  # Determine plot limits
  xmin <- min(plotData$bZX - plotData$seZX, na.rm = TRUE)
  xmax <- max(plotData$bZX + plotData$seZX, na.rm = TRUE)
  ymin <- min(plotData$bZY - plotData$seZY, na.rm = TRUE)
  ymax <- max(plotData$bZY + plotData$seZY, na.rm = TRUE)

  # Define regions (cis/trans)
  idx <- which(data$SMR[, 1] == data$probeID)
  if (length(idx) != 1 || is.na(data$SMR[idx, 8])) {
    stop(paste("ERROR: no SMR result for probe", data$probeID, "!"))
  }

  probeChr <- as.numeric(as.character(data$SMR[idx, 2]))
  probeBP <- as.numeric(as.character(data$SMR[idx, 3]))

  snpChr <- as.numeric(as.character(data$SNP[match(as.character(plotData[, 2]), data$SNP[, 1]), 2]))
  snpBP <- as.numeric(as.character(data$SNP[match(as.character(plotData[, 2]), data$SNP[, 1]), 3]))

  plotData$region <- ifelse(
    probeChr == snpChr & abs(snpBP - probeBP) < cisWindow * 1000,
    "cis",
    "trans"
  )

  # Highlight top SNP
  plotData$highlight <- 0
  max_snp_idx <- which.max((plotData$bZX / plotData$seZX)^2)
  plotData$highlight[max_snp_idx] <- 1

  # Create ggplot
  p <- plotData %>% dplyr::mutate(highlight=factor(highlight, levels=c(1,0),labels = c("top cis-eQTL","cis-eQTL"))) %>% 
    dplyr::filter(snpCorr^2>r2_threshold) %>%
    ggplot(aes(x = bZX, y = bZY)) +
    geom_point(aes(color = highlight, size = snpCorr^2, shape = highlight), alpha = 0.8) +
    geom_errorbar(aes(ymin = bZY - seZY, ymax = bZY + seZY), width = 0, alpha = 0.5) +
    geom_errorbarh(aes(xmin = bZX - seZX, xmax = bZX + seZX), height = 0, alpha = 0.5) +
    geom_abline(slope = plotData$bZY[max_snp_idx] / plotData$bZX[max_snp_idx],
                intercept = 0, linetype = "dashed", color = "orange") +
    scale_color_manual(values = c("red","black")) +
    scale_shape_manual(values = c(17,16)) +
    scale_size_continuous(range = c(1, 5)) +
    theme_classic() +
    theme(text=element_text(size=14,family="Helvetica"),axis.text=element_text(size=14,color="black"),axis.title=element_text(size=14),legend.text=element_text(size=12),legend.title=element_text(size=12))+
    labs(
      color = "SNP",
      x = "eQTL effect sizes",
      y = "GWAS effect sizes",
      size = expression(italic(r)^2),
      title = title
    )+
    guides(
        shape = "none",  # Remove the shape legend
        color = guide_legend(title = "", override.aes = list(size = 4,shape=c(17,16)) )
    )

  return(p)
}

# Helper function for drawing locus plot
draw_locus_plot <- function(query, gene, chr, pos_range, chrom, pos, title, label_list,force_label=FALSE,genetrack=FALSE){
    loc <- locus(data = query, seqname=chr, xrange=pos_range, flank = 3e5, 
        chrom=chrom, pos=pos, LD = "r2", 
        ens_db = "EnsDb.Hsapiens.v86")
    #loc <- link_LD(loc, token = "222a9926836b")
    if (length(label_list)>1){
        label_list = c(loc$index,label_list)
    }else if (length(label_list)==1 & label_list!="index") {
        label_list = c(loc$index,label_list)
    }
    if (min(loc$data$p,na.rm=T)> 5e-8 | !force_label){
        label_list=NULL
    }
    
    pg <- gg_scatter(loc, labels = unique(label_list),showLD=TRUE, nudge_y=0, nudge_x = -0.15, legend_pos="none") 
    pg <- pg + labs(title=title) + 
        theme(text=element_text(size=18, family="Helvetica"),
            axis.title=element_text(size=18),
            axis.text=element_text(size=18,color="black"),
            title=element_text(size=18)
            )
    if (-log10(min(loc$data$p,na.rm=T))<12){pg <- pg + scale_y_continuous(limits = c(0,12)) + geom_hline(yintercept=-log10(5e-8),color="grey", linetype="dashed")}
    if (genetrack){
        g <- gg_genetracks(loc,cex.text=0.9) + theme(text=element_text(size=4,family="Helvetica",color="black"))
        return(list(loc=loc,p=pg,g=g))
    }else{
        return(list(loc=loc, p=pg))
    }
}

# Helper function for drawing effect plot
draw_effect_plot <- function(tissue,GWAS_list,trait,probe,gene,chr,ancestry,celltype,title,tmp_dir,ref_dir=ref_dir,SMR_path=smr_executable,besd_dir=besd_dir,cisWindow=2000, r2_threshold=0.1){
    plot_dir <- file.path(tmp_dir, "plot")
    dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
    out_prefix <- file.path(plot_dir, paste0(trait,"_",tissue,"_",ancestry,"_",celltype))

    if (! file.exists(file.path(out_prefix,".txt"))){
        gwas_path=GWAS_list[GWAS_list$X1==trait,]$X4[1]
        if (grepl("\\.gz$", gwas_path)) {
            gwas_local=file.path(tmp_dir,str_replace(basename(gwas_path),".gz",""))
            if (!file.exists(gwas_local)){system(sprintf("zcat %s > %s", shQuote(gwas_path), shQuote(gwas_local)))}
        } else {
            gwas_local = gwas_path
        }
        refdata=file.path(ref_dir,ancestry,paste0("1KGP_",ancestry,"_chr",chr,"_maf01_mind05_hwe1en6_msite05_hg38"))
        command=paste0(SMR_path," --bfile ", refdata, " --gwas-summary ", gwas_local, " --beqtl-summary ", ifelse(celltype=="GTEx_Lung",paste0(besd_dir,"/Lung_eQTL_all_chr",chr),paste0(besd_dir,"/",tissue,"_",celltype,"_",ancestry,"_sceQTL")), " --out ", out_prefix, " --plot --probe ", probe, " --probe-wind 500 --gene-list ", file.path(metadata_dir, "glist-hg38"), " --diff-freq-prop 0.5")
        system(command)
    }
    SMRdata=ReadSMRData(paste0(out_prefix, ".", probe,".txt"))
    p <- SMREffectggPlot(data=SMRdata, title=title,cisWindow=cisWindow, transWindow=5000, pointsize=8, r2_threshold=0.1)
    return(p)
}

## 3. FIGURE 6: GWAS INTEGRATION FOR BLOOD (EAS)
# ------------------------------------------------------------------------------
print("--- Section 3: Generating Figure 6 - Blood (EAS) ---")

# --- 3.1: Load and Pre-filter SMR/COLOC Results ---
gwas_list_eas <- read_delim(file.path(gwas_lists_dir, "GWAS_EAS_list.txt"), col_names = FALSE, delim = '\t')

smr_results_eas <- readRDS(file.path(smr_coloc_results_dir, "Blood_smr_results.rds"))$smr_results1 %>%
    filter(p_SMR_bhadj < 0.05 & p_HEIDI > 0.05)
coloc_results_eas <- readRDS(file.path(smr_coloc_results_dir, "Blood_coloc_results.rds"))$coloc_results1 %>%
    filter(pp4 > 0.8)

smr_results_eur <- readRDS(file.path(smr_coloc_results_dir, "Blood_smr_results.rds"))$smr_results2 %>%
    filter(p_SMR_bhadj < 0.05 & p_HEIDI > 0.05)
coloc_results_eur <- readRDS(file.path(smr_coloc_results_dir, "Blood_coloc_results.rds"))$coloc_results2 %>%
    filter(pp4 > 0.8)

trait_list <- c("CD", "IBD", "RA", "SLE", "Asthma", "AR", "COPD")

# --- 3.2: Gene-Trait-Celltype Heatmap (Scatterpie) ---
smr_coloc_eas <- smr_results_eas %>%
    inner_join(coloc_results_eas %>% filter(pp4 > 0.9), by = c("trait", "celltype", "Gene")) %>%
    mutate(trait = factor(trait, levels = trait_list))

trait_matrix_eas <- smr_coloc_eas %>%
  distinct(Gene, trait) %>% arrange(trait, Gene) %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = trait, values_from = present, values_fill = 0) %>%
  column_to_rownames("Gene")

gene_clusters <- hclust(dist(trait_matrix_eas, method = "binary"))
gene_order <- gene_clusters$labels[gene_clusters$order]

pie_data_eas <- smr_coloc_eas %>%
  count(celltype, Gene, trait) %>% arrange(trait) %>%
  pivot_wider(names_from = trait, values_from = n, values_fill = 0) %>%
  mutate(
    celltype = factor(celltype, levels = celltypes_all_list$Blood[celltypes_all_list$Blood %in% smr_coloc_eas$celltype]),
    Gene = factor(Gene, levels = gene_order),
    x = as.numeric(celltype), y = as.numeric(Gene)
  )

n_traits <- ncol(smr_coloc_eas) - 4  # subtract non-trait columns

grid_data <- expand.grid(
  x = unique(pie_data$x),
  y = unique(pie_data$y)
)

ggplot() +
  # Add a grid using geom_tile (covers all positions, including NA)
  geom_tile(aes(x = x, y = y), data = grid_data, fill = NA, color = "lightgrey") +
  geom_scatterpie(aes(x = x, y = y), 
                  data = smr_coloc_eas,pie_scale=1.5,
                  cols = colnames(smr_coloc_eas)[3:(3+n_traits-1)],
                  color = NA) +
  scale_x_continuous(breaks = unique(smr_coloc_eas$x),
                    labels = unique(smr_coloc_eas$celltype),
                    expand = c(0.03,0.03)) +
  scale_y_continuous(breaks = unique(smr_coloc_eas$y),
                    labels = unique(smr_coloc_eas$Gene),
                    expand = c(0.01,0.01)) +
  scale_fill_brewer(palette="Set2") +
  #scale_fill_manual(values=colors) +
  coord_fixed(ratio=1) +
  theme_classic() +
  theme(axis.text=element_text(size=16,color="black"), 
  axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5), 
  axis.title=element_text(size=16),
  text=element_text(size=18,family="Helvetica"),
  axis.text.y = element_text(face="italic"),
  plot.title=element_text(hjust=0.5,size=20)) + 
  labs(fill="Trait", x="", y="Gene")
ggsave(file.path(output_dir, "fig6_smr_coloc_blood_eas_heatmap.png"), width=20, height=20)

# --- 3.3: Compare Gene Overlap (EUR vs. EAS) ---

calculate_overlaps <- function(trait_data) {
  smr1_genes <- trait_data$SMR1
  coloc1_genes <- trait_data$COLOC1
  smr2_genes <- trait_data$SMR2
  coloc2_genes <- trait_data$COLOC2

  # Calculate overlaps
  smr_both <- intersect(smr1_genes, smr2_genes)
  smr1_only <- setdiff(smr2_genes, smr_both)
  coloc_both <- intersect(coloc1_genes, coloc2_genes)
  coloc1_only <- setdiff(coloc2_genes, coloc_both)
  
  data.frame(
    Category = c("smr_both", "smr1_only", "coloc_both", "coloc1_only"),
    Count = c(length(smr_both), length(smr1_only), length(coloc_both), length(coloc1_only))
  )
}

gene_list = list()
for (tr in trait_list){
    gene_list[[tr]] = list(
        SMR1=smr_results_eas %>% dplyr::filter(trait==tr) %>% dplyr::pull(Gene) %>% as.character() %>% unique(),
        SMR2=smr_results_eur %>% dplyr::filter(trait==tr) %>% dplyr::pull(Gene) %>% as.character() %>% unique(),
        COLOC1=coloc_results_eas %>% dplyr::filter(trait==tr) %>% dplyr::pull(Gene) %>% unique(),
        COLOC2=coloc_results_eur %>% dplyr::filter(trait==tr) %>% dplyr::pull(Gene) %>% unique()
    )
}

# Create a dataframe with all traits
plot_data <- lapply(names(gene_list), function(trait) {
  df <- calculate_overlaps(gene_list[[trait]])
  df$Trait <- trait
  df
}) %>% bind_rows()

plot_data %>% dplyr::filter(Category %in% c("smr_both","smr1_only") )  %>% group_by(Trait) %>% reframe(prop=Count/sum(Count),Category=Category)

plot_data %>% dplyr::filter(Category %in% c("coloc_both","coloc1_only") )  %>% group_by(Trait) %>% reframe(prop=round(Count/sum(Count),6),Category=Category) %>% as.data.frame()

p1 <- plot_data %>% dplyr::filter(Category %in% c("coloc_both","coloc1_only") ) %>% 
dplyr::mutate(Category=factor(Category,levels=c("coloc_both","coloc1_only"),labels=c("EUR & EAS prioritized","EAS-specific")),
            Trait=factor(Trait,levels=trait_list)) %>%
ggplot(aes(y = fct_rev(Trait), x = Count, fill = Category)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  scale_fill_brewer(palette="Set2") +
  labs(y = "Trait",
       x = "Number of GWAS loci",
       fill = "", title="COLOC") +
  theme_classic() +
  theme(text = element_text(size=14,family="Helvetica"),
  axis.text=element_text(color="black",size=14),
  plot.title = element_text(hjust=0.5),
        legend.position = "none")

p2 <- plot_data %>% dplyr::filter(Category %in% c("smr_both","smr1_only") ) %>% 
dplyr::mutate(Category=factor(Category,levels=c("smr_both","smr1_only"),labels=c("EUR & EAS prioritized","EAS-specific")),
            Trait=factor(Trait,levels=trait_list)) %>%
ggplot(aes(y = fct_rev(Trait), x = Count, fill = Category)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  scale_fill_brewer(palette="Set2") +
  labs(y = "Trait",
       x = "Number of GWAS loci",
       fill = "",title="SMR") +
  theme_classic() +
  theme(text = element_text(size=14,family="Helvetica"),
        axis.text=element_text(color="black",size=14),
        plot.title = element_text(hjust=0.5),
        legend.position="bottom")


legend <- get_plot_component(p2+theme(legend.key.spacing.x = unit(1,"cm"),legend.text=element_text(family="Helvetica",size=14),legend.title=element_text(family="Helvetica",size=14)), pattern="guide-box-bottom", return_all=TRUE)

plot_grid(plot_grid(p1,p2+theme(legend.position="none",axis.text.y=element_blank())+labs(y=""),nrow=1,axis="b",align="h",rel_widths=c(1,0.8)),
    legend, nrow=2,align="v",rel_heights=c(1,0.1))
ggsave(file.path(output_dir, "fig6_smr_coloc_eur_vs_eas_overlap.png"), width=5.6, height=4)

# --- 3.4: Example Locus Analysis (SLC15A4 - SLE) ---
# --- A. Effect plots ---
tissue <- "Blood"
gene <- "SLC15A4"
anc <- "EAS"
tr <- "SLE"
besd_dir <- file.path(base_dir, tissue, "sc-eQTL", "results", "besd")
# Find intersecting cell types efficiently
intersecting_celltypes <- smr_results_eas %>%
    filter(trait == tr, Gene == gene) %>%
    inner_join(coloc_results_eas, by = c("trait", "celltype", "Gene")) %>%
    pull(celltype) %>%
    unique()
    
celltypes <- c(intersecting_celltypes, "JCTF_Japanese_blood", "eQTLGen")

# Locus plot generation with proper path handling
gwas_target <- select_gwas_snp(gene, snp_ids, GWAS_list_eas, tr, gene_bed, tmp_dir)
leadsnp <- gwas_target[which.min(p)]$SNP
ld_file <- file.path(tmp_dir, paste0("Blood_", anc, "_", leadsnp, ".ld"))
# Check if LD file exists, if not, generate it
if (!file.exists(ld_file)) {
    plink_bfile <- file.path(base_dir, "Blood/sc-eQTL/genotype",paste0("Blood_", anc, "_indivpruned_updated"))
    system(paste0("plink --bfile ", plink_bfile, " --ld-snp ", leadsnp, " --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0 --r2 --out ", file.path(tmp_dir, paste0("Blood_", anc, "_", leadsnp))))
}
ld <- fread(ld_file)
gwas_target <- merge(gwas_target, ld[,.(SNP_B,R2)], by.x = "SNP", by.y = "SNP_B", all.x = TRUE)
setnames(gwas_target,"R2","ld")
gwas_target <- gwas_target[!is.na(ld)]

# Define chromosome and position for the locus plot
chr=as.numeric(smr_results1[smr_results1$Gene==gene,]$ProbeChr[1])
probe = smr_results1[smr_results1$Gene==gene,]$probeID[1]
pos = gene_bed[gene_bed$Gene==probe,]$end
query_gene <- list()

# Loop over cell types and query eQTL data with error handling
for (c in celltypes[1:length(celltypes)]) {
    if(! c %in% c("JCTF_Japanese_blood","eQTLGen")){
    for (a in c("EAS","EUR")){
    result <- tryCatch({
        prefix_tmp = paste0(tissue,"_", c, "_",a,"_sceQTL")
        besd_dir_tmp = besd_dir 
        probe_tmp = probe
        df <- query_eqtl_data(
            prefix = prefix_tmp,
            probe = probe_tmp,
            besd_dir = besd_dir_tmp
        )
        if (min(df$p)<5e-8){
            leadsnp_tmp <- df[which.min(p)]$SNP
            if (!file.exists(paste0("Blood_",a,"_", leadsnp_tmp,".ld"))){
                tryCatch({
                system(paste0("plink --bfile ",paste0(base_dir,"Blood/sc-eQTL/genotype/Blood_",a,"_indivpruned_updated")," --ld-snp ", leadsnp_tmp, " --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0 --r2 --out ",paste0(tmp_dir,"Blood_",a,"_", leadsnp_tmp)))
                
                }, error=function(e){message("No valid variant")})
            }
            if (file.exists(file.path(tmp_dir, "Blood_",a,"_", leadsnp_tmp,".ld"))){
                ld <- fread(file.path(tmp_dir, "Blood_",a,"_", leadsnp_tmp,".ld"))
                ld_tmp <- fread(file.path(tmp_dir, "Blood_",a, "_", leadsnp_tmp,".ld"))
                df <- merge(df, ld_tmp[,.(SNP_B,R2)], by.x = "SNP", by.y = "SNP_B", all.x = TRUE)
                setnames(df,"R2","ld")
                df[!is.na(ld)]
            }else{
                df
            }
            
        }else{
            df
        }
    }, error = function(e) {
        message(paste("Error in query_eqtl_data for cell type", c, ":", e$message))
        NULL  # Return NULL on error
    })
    query_gene[[paste0(c,"-",a)]] <- result  # Store the result (or NULL) indexed by cell type
    }
    }else{
        if (c == "JCTF_Japanese_blood"){
            prefix_tmp = "JCTF_blood_eQTL.N1405"
            besd_dir_tmp = file.path(resource_dir, "eQTL/JCTF_Japanese_blood_eQTL/besd")
            epi = fread(paste0(besd_dir_tmp,"/",prefix_tmp,".epi"),select=2)$V2
            probe_tmp = epi[str_detect(epi,as.character(probe))]  
            a = "EAS"
        }else if (c=="eQTLGen"){
            prefix_tmp = paste0("eQTLGen_chr",chr)
            besd_dir_tmp = file.path(resource_dir, "eQTL/eQTLGen/cis-eQTL/GRCh38/")
            epi = fread(paste0(besd_dir_tmp,prefix_tmp,".epi"))
            probe_tmp = epi$V2[str_detect(epi$V2,as.character(probe))]
            a = "EUR"
        }
        result <- tryCatch({
        df <- query_eqtl_data(
            prefix = prefix_tmp,
            probe = probe_tmp,
            besd_dir = besd_dir_tmp
        )
        if (min(df$p)<5e-8){
            leadsnp_tmp <- df[which.min(p)]$SNP
            if (!file.exists(file.path(tmp_dir, "ref_",a,"_", leadsnp_tmp,".ld"))){
                tryCatch({
                system(paste0("plink --bfile ",file.path(ref_dir,"1000G","1KGP3_3202_GRCh38", a, paste0("1kGP_high_coverage_Illumina.chr",chr,".filtered.QCed.maf01.rmDup.",a))," --ld-snp ", leadsnp_tmp, " --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0 --r2 --out ",paste0(tmp_dir,"ref_",a,"_", leadsnp_tmp)))
                }, error=function(e){message("No valid variant")})
            }
            if (file.exists(file.path(tmp_dir, "ref_",a,"_", leadsnp_tmp,".ld"))){
                ld <- fread(file.path(tmp_dir, "ref_",a,"_", leadsnp_tmp,".ld"))
                ld_tmp <- fread(file.path(tmp_dir, "ref_",a, "_", leadsnp_tmp,".ld"))
                df <- merge(df, ld_tmp[,.(SNP_B,R2)], by.x = "SNP", by.y = "SNP_B", all.x = TRUE)
                setnames(df,"R2","ld")
                df[!is.na(ld)]
            }else{
                df
            }
        }else{
            df
        }
    }, error = function(e) {
        message(paste("Error in query_eqtl_data for cell type", c, ":", e$message))
        NULL  # Return NULL on error
    })
    query_gene[[paste0(c,"-",a)]] <- result  # Store the result (or NULL) indexed by cell type
    }
}


if(!file.exists(file.path(output_dir, paste0("fig6_locusplot_locusplot_",tissue,"_",anc,"_", tr, "_", gene, "_no_label.png")))){

        # Create the GWAS locus plot with error handling
        pg <- tryCatch({
            pg_result <- draw_locus_plot(gwas_target, chr = chr, pos_range = c(max(0, pos - 5e5), pos + 5e5),
            chrom = "chr", pos = "pos", title = paste0(tr," GWAS-EAS"), label_list = "index", force_label = FALSE, genetrack = TRUE)
            pg_result  # Return the plot object
        }, error = function(e) {
            message(paste("Error in draw_locus_plot for GWAS:", e$message))
            NULL  # Return NULL on error
        })

        pqs <- list()
        # Loop over query_gene and create locus plots with error handling
        for (x in seq_along(query_gene)) {
            ctype <- names(query_gene)[x]
            df <- query_gene[[x]]
            
            if (!is.null(df)) {
            plot_result <- tryCatch({
                p <- draw_locus_plot(df, chr = chr, pos_range = c(max(0, pos - 5e5), pos + 5e5),
                chrom = "Chr", pos = "BP", title = paste0(tissue,"-",sub("_blood","",ctype)), label_list = pg$loc$index_snp, force_label = FALSE,genetrack = FALSE)
                p  # Return the plot object
            }, error = function(e) {
                message(paste("Error in draw_locus_plot for cell type", ctype, ":", e$message))
                NULL  # Return NULL on error
            })
            pqs[[ctype]] <- plot_result
            } else {
            message(paste("No data available for cell type", ctype))
            pqs[[ctype]] <- NULL
            }
        }

        plot_list <- list()
        if (!is.null(pg)) {
            plot_list <- c(plot_list, list(pg$p),list(pg$g))
        }
        # Extract plots from pqs and exclude NULLs
        pqs_plots <- lapply(pqs, function(p) if (!is.null(p)) p$p else NULL)
        pqs_plots <- Filter(Negate(is.null), pqs_plots)
        
        if (length(pqs_plots) > 0) {
            plot_list <- c(plot_list, pqs_plots)
        }

        for (i in 1:floor(length(plot_list)/2) ){
            plot_list[[i*2]] = plot_list[[i*2]] + theme(axis.title.y=element_blank(),) 
        }
        for (i in 1:(length(plot_list)-2)){
            plot_list[[i]] = plot_list[[i]] + theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) 
        }
        plot_list = lapply(plot_list, function(x) x + theme(plot.margin=margin(b=0.5,t=0.8,unit="cm")))
        # Check if there are plots to display
        if (length(plot_list) > 0) {
            # Plot and save the combined figure
            plot_grid(plotlist = plot_list, ncol = ifelse(length(plot_list)>5,2,1), align = "v", rel_heights=c(rep(1,ceiling(length(plot_list)/2)-1),1.15))
            ggsave(
            filename = file.path(output_dir, paste0("fig6_locusplot_locusplot_",tissue,"_",anc,"_", tr, "_", gene, "_no_label.png")),
            width = ifelse(length(plot_list)>5,10,5), height = ceiling(length(plot_list)/2)*3
            )
        } else {
            message("No plots were generated due to errors.")
        }
}

## 4. FIGURE 7: GWAS INTEGRATION FOR LUNG (EUR)
# ------------------------------------------------------------------------------
print("--- Section 4: Generating Figure 7 - Lung (EUR) ---")
# --- 4.1: Load SMR/COLOC Results ---
trait_list <- c("IPF","COPD","CPD","LuC","COVID19-HGI-A2", "COVID19-HGI-B1", "COVID19-HGI-B2", "COVID19-HGI-C2", "CAD","Pleurisy","Bronchiectasis","Chronic Bronchitis","Bronchitis","Chronic Sinusitis","Ped_asthma","Asthma","Pneumoconiosis","Pneumonia","Pneumothorax","PTB","SAS","Tonsillitis")

colors <-  c("dodgerblue", "firebrick", "forestgreen", "gold", "darkorchid",
              "orange", "darkturquoise", "deeppink", "limegreen", 
              "mediumpurple", "royalblue", "darkorange", "seagreen", 
              "violet", "chocolate", "steelblue", "darkmagenta", "olivedrab",
              "hotpink", "sienna", "deepskyblue", "darkcyan", "darkgoldenrod")

names(colors) <- trait_list
ancestry="EUR"
smr_results <- readRDS(file.path(smr_coloc_results_dir,"Lung_smr_results.rds"))
smr_results <- smr_results %>% dplyr::filter(celltype != "T") %>% left_join(smr_list_all %>% rename(celltype=cell_type) %>% dplyr::filter(tissue%in% c("Lung","Blood")) %>% dplyr::select(tissue,celltype,cell_type_short) %>% distinct(),by=c("tissue","celltype"))
smr_results$cell_type_short[smr_results$celltype=="GTEx_Lung"] = "GTEx_Lung"
smr_results <- smr_results %>% dplyr::filter(!is.na(cell_type_short))
smr_results$trait[smr_results$trait=="Chronic_bronchitis"] = "Chronic Bronchitis" 
smr_results$trait[smr_results$trait=="Chronic_Sinusitis"] = "Chronic Sinusitis" 

coloc_results <- readRDS(file.path(smr_coloc_results_dir,"Lung_coloc_results.rds"))
coloc_results <- coloc_results %>% dplyr::filter(celltype != "T")%>% left_join(smr_list_all %>% rename(celltype=cell_type) %>% dplyr::filter(tissue%in% c("Lung","Blood")) %>% dplyr::select(tissue,celltype,cell_type_short) %>% distinct(),by=c("tissue","celltype"))
coloc_results$cell_type_short[coloc_results$celltype=="GTEx_Lung"] = "GTEx_Lung"
coloc_results <- coloc_results %>% dplyr::filter(!is.na(cell_type_short))
colnames(coloc_results)[2] <- "Gene"
coloc_results$trait[coloc_results$trait=="Chronic_bronchitis"] = "Chronic Bronchitis" 
coloc_results$trait[coloc_results$trait=="Chronic_Sinusitis"] = "Chronic Sinusitis" 

smr_results <- smr_results %>% dplyr::filter(p_SMR_bhadj<0.05 & p_HEIDI>0.05)
coloc_results <- coloc_results %>% dplyr::filter(pp4>0.8)

# --- 4.2: Gene-Trait-Celltype Heatmap (Scatterpie) ---
smr_coloc <- smr_results %>% inner_join(coloc_results %>% dplyr::filter(pp4>0.9), by=c("tissue","trait","celltype","Gene","cell_type_short")) %>% dplyr::filter(celltype != "GTEx_Lung" ) %>% dplyr::mutate(cell_type_short=factor(cell_type_short,levels=c(unique(c(cell_type_short_dict[['Blood']],cell_type_short_dict[['Lung']])),"GTEx_Lung"))) %>% 
arrange(tissue,cell_type_short) 
trait_matrix <- smr_coloc %>%
  filter(celltype != "GTEx_Lung") %>%
  dplyr::mutate(trait=factor(trait,levels=trait_list))%>%
  distinct(Gene, trait) %>%          # Remove duplicate gene-trait pairs
  arrange(trait,Gene) %>%
  mutate(present = 1) %>%            # Binary indicator
  pivot_wider(
    names_from = trait,
    values_from = present,
    values_fill = 0
  ) %>%
  column_to_rownames("Gene")         # Convert to matrix format

gene_clusters <- hclust(dist(trait_matrix, method = "binary"))  # Jaccard distance
gene_order <- gene_clusters$labels[gene_clusters$order]         # Extract order

# Prepare data in wide format for scatterpie
p_list=list()
for (ts in c("Lung","Blood")){
    pie_data <- smr_coloc%>% 
    dplyr::filter(celltype != "GTEx_Lung" & tissue==ts) %>%
    count(tissue,cell_type_short, Gene, trait) %>% 
    dplyr::mutate(trait=factor(trait,levels=trait_list),
            cell_type_short=factor(cell_type_short, levels=cell_type_short_dict[[ts]][cell_type_short_dict[[ts]] %in% smr_results$cell_type_short])) %>% 
            arrange(trait) %>%
    tidyr::pivot_wider(names_from = trait, values_from = n, values_fill = 0) %>%
    arrange(tissue,cell_type_short) %>%
    mutate(ts_anc_ct = paste0(tissue,"-",ancestry,"-",cell_type_short)) 

    pie_data <- pie_data %>% dplyr::mutate(
            ts_anc_ct = factor(ts_anc_ct,levels=unique(pie_data$ts_anc_ct)),
            Gene=factor(Gene,levels=gene_order),
            x = as.numeric(ts_anc_ct),
            y = as.numeric(Gene))

    n_traits <- ncol(pie_data) - 6  # subtract non-trait columns

    grid_data <- expand.grid(
    x = unique(pie_data$x),
    y = 1:length(levels(pie_data$Gene))
    )

    p_list[[ts]]=ggplot() +
    # Add a grid using geom_tile (covers all positions, including NA)
    geom_tile(aes(x = x, y = y), data = grid_data, fill = NA, color = "lightgrey") +
    geom_scatterpie(aes(x = x, y = y), 
                    data = pie_data, pie_scale=ifelse(ts=="Blood",1,1.4),
                    cols = colnames(pie_data)[4:(4+n_traits-1)],
                    color = NA) +
    scale_x_continuous(breaks = unique(pie_data$x),
                        labels = unique(pie_data$cell_type_short),
                        expand = c(0.03,0.03)) +
    scale_y_continuous(breaks = unique(grid_data$y),
                        labels = levels(pie_data$Gene),
                        expand = c(0.01,0.01)) +
    #scale_fill_brewer(palette="Set3") +
    facet_wrap(~tissue, labeller = as_labeller(c(`Blood`="Blood-EUR",`Lung`="Lung-EUR")) ) +
    scale_fill_manual(values=colors) +
    coord_fixed(ratio=1) +
    theme_classic() +
    theme(axis.text.x = element_text(size=18,color="black",angle = 90, hjust =1,vjust=0.5), 
            axis.text.y = element_text(size=18,color="black",face="italic"),
            text=element_text(family="Helvetica"),
            axis.title = element_text(size=20),
            legend.position="none",
            strip.background = element_blank(),strip.text = element_blank()) + 
    labs(fill="Trait", x="", y="Gene")
}
plot_grid(p_list[["Lung"]],
    p_list[["Blood"]] + theme(axis.text.y=element_blank(),
        axis.title.y=element_blank(),axis.ticks.y=element_blank(),
        plot.margin=margin(l=0,unit="cm")),nrow=1, align="h", axis="lrbt")

ggsave(file.path(output_dir, "fig7_smr_coloc_lung_blood_heatmap.png"), width=12, height=15)

# --- 4.3: eGene Count Comparison & Overlap ---
smr_coloc_negene <- rbind(
smr_results %>% dplyr::filter(tissue=="Lung" & celltype!="GTEx_Lung" ) %>%
  group_by(cell_type_short) %>% 
  reframe(nGenes = n_distinct(Gene), method="SMR") %>% 
  dplyr::mutate(cell_type_short=factor(cell_type_short,levels=c(cell_type_short_dict[["Lung"]][celltypes_all_list[["Lung"]]],"GTEx_Lung"))) ,
coloc_results %>% dplyr::filter(tissue=="Lung" & celltype!="GTEx_Lung") %>%
  group_by(cell_type_short) %>% 
  reframe(nGenes = n_distinct(Gene),method="COLOC") %>% 
  dplyr::mutate(cell_type_short=factor(cell_type_short,levels=c(cell_type_short_dict[["Lung"]][celltypes_all_list[["Lung"]]],"GTEx_Lung")))
)

smr_coloc_negene %>%
  ggplot(aes(x=nGenes, y=fct_rev(cell_type_short), color=cell_type_short)) + geom_point(size=2.5) + 
  scale_color_manual(values=c(cell_type_color_dict[cell_type_short_dict[["Lung"]]],GTEx_Lung="grey")) +
  theme_classic() + labs(x="Number of eGenes",y="Cell type") + 
  facet_wrap(~method,axes="all", axis.labels="margins",scales="free_x") +
  theme(text=element_text(size=14,family="Helvetica"),
        axis.text=element_text(size=14,color="black"),
        axis.title=element_text(size=14),
        axis.ticks.y=element_blank(),
        strip.text=element_text(size=14),
        strip.background=element_rect(color="white",fill="white"),
        legend.position="none", 
        panel.grid.major.y=element_line(linetype="dashed"))
ggsave(file.path(output_dir, "fig7_smr_coloc_lung_negene_counts.png"),width=5,height=4)

smr_coloc <- smr_results %>% inner_join(coloc_results, by=c("tissue","trait","celltype","Gene","cell_type_short")) %>% dplyr::filter(celltype != "GTEx_Lung" ) %>% dplyr::mutate(cell_type_short=factor(cell_type_short,levels=c(unique(c(cell_type_short_dict[['Blood']],cell_type_short_dict[['Lung']])),"GTEx_Lung"))) %>% 
arrange(tissue,cell_type_short) 

gene_list = list()
for (tr in unique(smr_coloc$trait)){
    gene_list[[tr]] = list(
        blood=smr_coloc %>% dplyr::filter(trait==tr & celltype != "GTEx_Lung" & tissue=="Blood") %>% dplyr::mutate(tr_gene=paste0(trait,"-",Gene)) %>% dplyr::pull(tr_gene) %>% unique(),
        lung=smr_coloc %>% dplyr::filter(trait==tr & celltype != "GTEx_Lung" & tissue=="Lung") %>% dplyr::mutate(tr_gene=paste0(trait,"-",Gene)) %>% dplyr::pull(tr_gene) %>% unique()
    )
}

tr_gene_list = list(Lung= lapply(gene_list,function(x) x$lung) %>% unlist() %>% unique(),Blood=lapply(gene_list,function(x) x$blood) %>% unlist() %>% unique())
ggvenn(tr_gene_list, fill_color=c("#3561db","#efab45"),set_name_size=0,text_size=5.2,fill_alpha=0.8,stroke_size=1.2) + theme(text=element_text(family="Helvetica"))
ggsave(file.path(output_dir, "fig7_lung_vs_blood_gene_venn.png"), width=4.8,height=4.8)

# --- 4.4: Example Locus Analysis (FAM13A - IPF) ---
anc="EUR"
GWAS_list=read_delim(file.path(gwas_lists_dir, "Blood/script_post_eqtl/GWAS_",anc,"_list.txt"), col_names = F, delim='\t') 
ts="Lung"
tr="IPF"
gene="FAM13A"
print(tr)
snp_ids <- load_gwas_snp(GWAS_list,tr)

print(gene)
ct1 <- smr_results %>% dplyr::filter(trait==tr & Gene==gene) %>% distinct()
ct2 <- coloc_results %>% dplyr::filter(trait==tr & Gene==gene)  %>% distinct()
celltypes = ct1 %>% inner_join(ct2 %>% dplyr::rename(probeID=probe_name)%>%dplyr::select(probeID,Gene,pp4,trait,celltype,tissue,cell_type_short))

#### SMR Effect plot ------------------
chr=celltypes$ProbeChr[1]
probe = celltypes$probeID[1] 
pos = gene_bed[gene_bed$Gene==probe,]$end

for (x in 1:nrow(celltypes)){
    ts = celltypes$tissue[x]
    ct = celltypes$celltype[x]
    anc = celltypes$ancestry[x]
    besd_dir=file.path(base_dir,ts,"/sc-eQTL/results/besd")
        draw_effect_plot(ts,GWAS_list,trait=tr,probe=probe,gene=gene,ancestry=anc,celltype=ct,chr=chr,title=bquote(paste(.(ts), "-", .(cell_type_short_dict[[ts]][ct]),"-",.(anc))),tmp_dir=tmp_dir,r2_threshold=0.1, cisWindow=2000,ref_dir=ref_dir,SMR_path=SMR,besd_dir=ifelse(celltypes[x] == "GTEx_Lung", "eQTL/gtex_resources_besd/eQTL_hg38/eQTL_besd_hg38", besd_dir))
        ggsave(file.path(output_dir, paste0("fig7_SMR_effectplot_",tr,"_",ts,"_",anc,"_",gene,"_",ct,".png")), width=4.5, height=3)
}


#### Locus plot ------------------
gwas_target <- select_gwas_snp(gene,snp_ids,GWAS_list,tr,gene_bed)
leadsnp <- gwas_target[which.min(p)]$SNP
ld_out_prefix <- file.path(tmp_dir, paste0("Lung_EUR_", leadsnp))
if (!file.exists(paste0(ld_out_prefix, ".ld"))){
    system(paste0("plink --bfile Lung/sc-eQTL/genotype/Lung_EUR_indivpruned_updated --ld-snp ", leadsnp, " --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0 --r2 --out", ld_out_prefix))
    #ref_bfile = paste0("1000G/1KGP3_3202_GRCh38/", anc, "/1kGP_high_coverage_Illumina.chr",chr,".filtered.QCed.maf01.rmDup.",anc)
    #system(paste0("plink --bfile ",ref_bfile, " --ld-snp ", leadsnp, " --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0 --r2 --out ref_",anc,"_", leadsnp))
}
ld <- fread(paste0(ld_out_prefix,".ld"))
#ld <- fread(paste0("ref_",anc,"_", leadsnp,".ld"))
gwas_target <- merge(gwas_target, ld[,.(SNP_B,R2)], by.x = "SNP", by.y = "SNP_B", all.x = TRUE)
setnames(gwas_target,"R2","ld")
gwas_target <- gwas_target[!is.na(ld)]
query_gene <- list()
if(!"GTEx_Lung" %in% celltypes$celltype){
    celltypes = bind_rows(celltypes,tibble(`probeID`=probe,`ProbeChr`=chr,`Gene`=gene,`tissue`="Lung",`celltype`="CD4T",`cell_type_short`="CD4T",ancestry="EUR"))
    celltypes = bind_rows(celltypes,tibble(`probeID`=probe,`ProbeChr`=chr,`Gene`=gene,`tissue`="Blood",`celltype`="Treg",`cell_type_short`="Treg",ancestry="EUR"))
    celltypes = bind_rows(celltypes,tibble(`probeID`=probe,`ProbeChr`=chr,`Gene`=gene,`tissue`="Lung",`celltype`="GTEx_Lung",`cell_type_short`="GTEx_Lung",ancestry="EUR"))
    celltypes = bind_rows(celltypes,tibble(`probeID`=probe,`ProbeChr`=chr,`Gene`=gene,`tissue`="Blood",`celltype`="eQTLGen",`cell_type_short`="eQTLGen",ancestry="EUR"))
}
# Loop over cell types and query eQTL data with error handling
for (x in 1:nrow(celltypes)) {
    ts = celltypes$tissue[x]
    c = celltypes$cell_type[x]
    a = celltypes$ancestry[x]
    result <- tryCatch({
        if(c=="GTEx_Lung"){
            prefix_tmp = paste0("Lung_eQTL_all_chr",chr)
            besd_dir_tmp = file.path(resource_dir,"eQTL/gtex_resources_besd/eQTL_hg38/eQTL_besd_hg38/")
            epi = fread(paste0(besd_dir_tmp,prefix_tmp,".epi"))
            probe_tmp = epi$V2[str_detect(epi$V2,as.character(probe))]
        }else if(c=="eQTLGen") {
            prefix_tmp = paste0("eQTLGen_chr",chr)
            besd_dir_tmp =file.path(resource_dir,"eQTL/eQTLGen/cis-eQTL/GRCh38/" )
            epi = fread(paste0(besd_dir_tmp,prefix_tmp,".epi"))
            probe_tmp = epi$V2[str_detect(epi$V2,as.character(probe))]
        }else{
            prefix_tmp = paste0(ts,"_", c, "_",a,"_sceQTL")
            besd_dir_tmp = file.path(base_dir,ts,"/sc-eQTL/results/besd/") 
            probe_tmp = gene_bed[gene_bed$gene_name==gene,]$Gene
        }
        df <- query_eqtl_data(
            prefix = prefix_tmp,
            probe = probe_tmp,
            besd_dir = besd_dir_tmp
        )
        #df <- df[SNP %in% query_gene[[1]]$SNP]
        if (min(df$p)<5e-8){
            #a <- "EUR"
            leadsnp_tmp <- df[which.min(p)]$SNP
            if (!file.exists(paste0(ts,"_",a,"_", leadsnp_tmp,".ld"))){
                system(paste0("plink --bfile ",ts,"/sc-eQTL/genotype/",ts,"_",a,"_indivpruned_updated --ld-snp ", leadsnp_tmp, " --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0 --r2 --out ",ts,"_",a,"_", leadsnp_tmp))
            }
            ld <- fread(file.path(tmp_dir,paste0(ts,"_",a,"_", leadsnp_tmp,".ld")))
            ld_tmp <- fread(file.path(tmp_dir,paste0(ts,"_",a, "_", leadsnp_tmp,".ld")))
            df <- merge(df, ld_tmp[,.(SNP_B,R2)], by.x = "SNP", by.y = "SNP_B", all.x = TRUE)
            setnames(df,"R2","ld")
            df[!is.na(ld)]
        }else{
            df
        }
    }, error = function(e) {
        message(paste("Error in query_eqtl_data for cell type", c, ":", e$message))
        NULL  # Return NULL on error
    })
    cout = ifelse(c %in% c("GTEx_Lung","eQTLGen"),c,cell_type_short_dict[[ts]][c])
    aout = ifelse(c %in% c("GTEx_Lung","eQTLGen"),"EUR",a) 
    query_gene[[paste0(ts,"-",a,"-",cout)]] <- result  # Store the result (or NULL) indexed by cell type
}



if(!file.exists(file.path(output_dir, paste0("fig7_locusplot_Lung_EUR_", tr, "_", gene, ".png")))){
    # Create the GWAS locus plot with error handling
    #gwas_target = gwas_target[SNP %in% query_gene[[1]]$SNP]
    idx_list = lapply(query_gene, function(q) {
        if (min(q$p)<5e-8){q$SNP[which.min(q$p)]}
    }) %>% unlist()
    anc="EUR"
    pg <- tryCatch({
        pg_result <- draw_locus_plot(gwas_target, chr = chr, pos_range = c(max(0, pos - 5e5), pos + 5e5),
        chrom = "chr", pos = "pos", title = paste0(tr," GWAS-",anc), label_list = "index", force_label = FALSE, genetrack = TRUE)
        pg_result  # Return the plot object
    }, error = function(e) {
        message(paste("Error in draw_locus_plot for GWAS:", e$message))
        NULL  # Return NULL on error
    })
     
    pqs <- list()
    # Loop over query_gene and create locus plots with error handling
    for (x in seq_along(query_gene)) {
        ts <- str_split(names(query_gene)[x],"-")[[1]][1]
        anc <-  str_split(names(query_gene)[x],"-")[[1]][2] 
        ctype <- sub(paste0(ts,"-",anc,"-"), "",names(query_gene)[x])
        
        #ctype <- ifelse(ctype!="GTEx_Lung",paste0("Lung-",cell_type_short_dict[["Lung"]][ctype]),ctype)
        df <- query_gene[[x]]
        
        if (!is.null(df)) {
        plot_result <- tryCatch({
            p <- draw_locus_plot(df, chr = chr, pos_range = c(max(0, pos - 5e5), pos + 5e5),
            chrom = "Chr", pos = "BP", title = names(query_gene)[x], label_list = pg$loc$index_snp,force_label = FALSE, genetrack = FALSE)
            p  # Return the plot object
        }, error = function(e) {
            message(paste("Error in draw_locus_plot for cell type", ctype, ":", e$message))
            NULL  # Return NULL on error
        })
        pqs[[paste0(ts,"-",anc,"-",ctype)]] <- plot_result
        } else {
        message(paste("No data available for cell type", ctype))
        pqs[[paste0(ts,"-",anc,"-",ctype)]] <- NULL
        }
    }

    plot_list <- list()
    if (!is.null(pg)) {
        plot_list <- c(plot_list, list(pg$p), list(pg$g))
    }
    # Extract plots from pqs and exclude NULLs
    pqs_plots <- lapply(pqs, function(p) if (!is.null(p)) p$p else NULL)
    pqs_plots <- Filter(Negate(is.null), pqs_plots)
    
    if (length(pqs_plots) > 0) {
        plot_list <- c(plot_list, pqs_plots)
    }
    #if (!is.null(pg)) {
    #    plot_list <- c(plot_list)
    #}
    for (i in 1:floor(length(plot_list)/2) ){
        plot_list[[i*2]] = plot_list[[i*2]] + theme(axis.title.y=element_blank(),) 
    }
    for (i in 1:(length(plot_list)-2)){
        plot_list[[i]] = plot_list[[i]] + theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) 
    }
    # Check if there are plots to display
    if (length(plot_list) > 0) {
        # Plot and save the combined figure
        plot_list = lapply(plot_list, function(p) p+theme(plot.margin=margin(b=0.5,t=0.8,unit="cm")))
        plot_grid(plotlist = plot_list, ncol = ifelse(length(plot_list)>5,2,1), align = "v", rel_heights=c(rep(1,ceiling(length(plot_list)/2)-1),1.15))
        ggsave(
            filename = file.path(output_dir, paste0("fig7_locusplot_Lung_EUR_", tr, "_", gene, ".png")),
            width = ifelse(length(plot_list)>5,10,5), height = ceiling(length(plot_list)/2)*3
        )

    } else {
        message("No plots were generated due to errors.")
    }
}
```
## Expression heatmap
#### expression + SMR + COLOC heatmap ---------------------
gene="FAM13A"
tr="IPF"
for (ts in c("Blood","Lung")){
    gene_expr_df = fread(file.path(base_dir,"expression_decile",paste0(ts,"_EUR_cross_celltype.txt")))
    expr = gene_expr_df[gene_name==gene,-1]
    ct <- sapply(colnames(expr), function(x) str_replace(str_replace(x,"expr_mean_",""),"_EUR",""))
    stats <- data.frame(cell_type_short=cell_type_short_dict[[ts]][ct],expr=t(expr[1]))
    stats <- stats %>% left_join(smr_results %>% dplyr::filter(trait==tr & Gene==gene & tissue==ts)%>%dplyr::select(cell_type_short,p_SMR),by="cell_type_short")
    stats <- stats %>% left_join(coloc_results %>% dplyr::filter(trait==tr & Gene==gene & tissue==ts)%>%dplyr::select(cell_type_short,pp4),by="cell_type_short")

    long_data <- stats %>%
    pivot_longer(cols = -cell_type_short, names_to = "variable", values_to = "value") %>%
    mutate(
        # Convert p-values to -log10 scale for better visualization
        value = ifelse(variable == "p_SMR" & !is.na(value), -log10(value), value),
        # Format cell types as factors to preserve order
        cell_type_short = factor(cell_type_short, levels = rev(cell_type_short_dict[[ts]]))
    )


    long_data <- long_data %>% dplyr::mutate(variable=factor(variable,levels=c("expr","pp4","p_SMR"),labels=c("Mean expression","coloc PP4","-log_10(P_SMR)")))  

    p <- ggplot(long_data) +
    # First layer: expr (orange scale)
    geom_tile(
        data = . %>% filter(variable == "Mean expression"),
        aes(x = 1, y = cell_type_short, fill = value),
        color = "black", na.rm = TRUE
    ) +
    scale_fill_gradientn(
        colours = c("white", "#4965b0"),
        na.value = "white",
        name = "Mean expression",
        limits = c(0, 0.6),
        guide = guide_colorbar(barwidth = 5, barheight = 0.8, title.position="top",order=1) # Control bar size
    ) +

    # Second layer: pp4 (purple scale)
    ggnewscale::new_scale_fill() +
    geom_tile(
        data = . %>% filter(variable == "coloc PP4"),
        aes(x = 1, y = cell_type_short, fill = value),
        color = "black", na.rm = TRUE
    ) +
    scale_fill_gradientn(
        colours = c("white", "#a30543"),
        na.value = "white",
        limits = c(0, 1),
        name = "coloc PP4",
        guide = guide_colorbar(barwidth = 5, barheight = 0.8, title.position="top",order=2) # Control bar size
    ) +

    # Third layer: p_SMR (pink scale)
    ggnewscale::new_scale_fill() +
    geom_tile(
        data = . %>% filter(variable == "-log_10(P_SMR)"),
        aes(x = 1, y = cell_type_short, fill = value),
        color = "black", na.rm = TRUE
    ) +
    scale_fill_gradientn(
        colours = c("white", "#f36f43"),
        na.value = "white",
        limits = c(1, 7),
        name = "-log<sub>10</sub>(P<sub>SMR</sub>)",
        guide = guide_colorbar(barwidth = 5, barheight = 0.8, title.position="top",order=3) # Control bar size
    ) +

    # Faceting and theming
    facet_wrap(~variable, nrow = 1,
                labeller = as_labeller(c(
                `Mean expression` = "Mean expression",
                `coloc PP4` = "coloc PP4",
                `-log_10(P_SMR)` = "-log<sub>10</sub>(P<sub>SMR</sub>)"
                )),
                strip.position = "top") +
    theme_minimal() +
    theme(
        text = element_text(size = 12, family = "Helvetica"),
        axis.text=element_text(size=12,color="black"),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        panel.grid = element_blank(),
        plot.title=element_text(size=16,hjust=0.5),
        legend.position = "bottom",
        legend.text = element_text(size=8,angle=90),
        legend.title = element_markdown(size = 10, margin = margin(r = 10),family="Helvetica"), # Add space between title and bar
        strip.text = element_markdown(size=12,angle = 90, hjust = 0, vjust = 0.5, family="Helvetica"),
        #-- This is the key change for stacking legends --#
        legend.box = "vertical",
        legend.box.just = "left", # Align the block of legends to the left
        legend.spacing = unit(0.15, "lines") # Adjust vertical space between legends
    ) + 
    coord_fixed(ratio = 1)

    if (ts=="Lung"){ p <- p+theme(legend.position="none")}

    ggsave(file.path(output_dir, paste0("fig7_expr_smr_coloc_heatmap_",ts,"_EUR_",tr,"_",gene,"_cross_celltype.png")),width=4.5, height=7)
}
