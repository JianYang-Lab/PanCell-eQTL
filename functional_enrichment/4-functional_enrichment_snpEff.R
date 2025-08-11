##############################################################################
# Description: Conduct enrichment test on eSNPs VS null
##############################################################################
library(tidyverse)
library(data.table)
library(dplyr)
library(stringr)
library(parallel)

library(GenomicRanges)

preprocess_controls <- function(index_snps, gene_annotations, background_snps,
                                maf_breaks = c(0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5),
                                tss_breaks = c(0, 1e5, 2e5, 3e5, 4e5, 5e5, 6e5, 7e5, 8e5, 9e5, 10e5),
                                n_cores = 8) {
  
  # Create genomic ranges once
  genes_gr <- makeGRangesFromDataFrame(
    gene_annotations,
    seqnames.field = "chr",
    start.field = "TSS",
    end.field = "TSS"
  )
  
  index_gr <- makeGRangesFromDataFrame(
    index_snps,
    seqnames.field = "chr",
    start.field = "pos",
    end.field = "pos"
  )
  
  # Precompute nearest genes and bins for index SNPs
  nearest_tss_index <- distanceToNearest(index_gr, genes_gr)
  index_snps[, phenotype_id := gene_annotations[subjectHits(nearest_tss_index), phenotype_id]]
  index_snps[, tss_distance := mcols(nearest_tss_index)$distance]
  index_snps[, maf_bin := cut(afreq, maf_breaks, include.lowest = TRUE)]
  index_snps[, tss_bin := cut(tss_distance, tss_breaks, include.lowest = TRUE)]
  
  # Precompute control pools for all genes
  unique_genes <- unique(index_snps$phenotype_id)
  gene_control_pools <- mclapply(unique_genes, function(gene_id) {
    gene_info <- gene_annotations[phenotype_id == gene_id]
    chr_controls <- background_snps[chr == gene_info$chr]
    
    if (nrow(chr_controls) == 0) return(NULL)
    
    chr_controls[, tss_distance := abs(as.numeric(pos) - gene_info$TSS)] %>%
      .[, `:=`(
        maf_bin = cut(afreq, maf_breaks, include.lowest = TRUE),
        tss_bin = cut(tss_distance, tss_breaks, include.lowest = TRUE)
      )] %>%
      .[, .(variant_id, maf_bin, tss_bin)] %>% .[!is.na(tss_bin)]
  }, mc.cores = n_cores) %>%
    set_names(unique_genes)
  
  return(list(
    index_attributes = index_snps[, .(variant_id, phenotype_id, maf_bin, tss_bin)],
    control_pools = gene_control_pools

  ))
}

# Optimized matching using precomputed data
fast_match_controls <- function(precomputed, n_cores = 8) {
  index_attrs <- precomputed$index_attributes
  control_pools <- precomputed$control_pools

  # Split index SNPs by gene for parallel processing
  matched <- mclapply(split(index_attrs, by = "phenotype_id"), function(gene_group) {
    gene_id <- unique(gene_group$phenotype_id)
    pool <- control_pools[[gene_id]]
    
    if (is.null(pool)) return(gene_group[, .(variant_id)])
    
    gene_group[, control_id := {
      matched_pool <- pool[maf_bin == .BY$maf_bin & tss_bin == .BY$tss_bin, variant_id]
      if (length(matched_pool) > 0) sample(matched_pool, .N, replace = TRUE) else NA_character_
    }, by = .(maf_bin, tss_bin)]
    
    gene_group

  }, mc.cores = n_cores) %>%
    rbindlist()
  
  return(matched)
}

# Modified enrichment calculation using precomputed data
calculate_chromatin_enrichment_optimized <- function(precomputed, chromatin_annotation,
                                                     n_iterations = 1000, n_cores = 8) {
  # Precompute observed proportion
  index_annotated <- merge(precomputed$index_attributes,
                          chromatin_annotation,
                          by = "variant_id",
                          all.x = TRUE)
  observed_prop <- sum(!is.na(index_annotated$Annotation)) / nrow(index_annotated)
  
  # Parallel bootstrap using precomputed control pools
  expected_prop <- mclapply(1:n_iterations, function(i) {
    if (i %% 10 == 0) message(paste("Completed", i, "iterations"))
    sampled_controls <- fast_match_controls(precomputed, n_cores = 1)  # Single-core for bootstrap

    sum(sampled_controls$control_id %in% chromatin_annotation$variant_id, na.rm = TRUE) / nrow(sampled_controls)
  }, mc.cores = n_cores) %>% unlist()
  
  # Calculate statistics
  expected_mean <- mean(expected_prop)
  expected_var <- var(expected_prop) 
  print(paste0("Successful iterations: ",length(expected_prop)))
  data.table(
    chromatin_state = unique(chromatin_annotation$Annotation),
    observed = observed_prop,
    expected = expected_mean,
    fold_enrichment = observed_prop / expected_mean,
    enrichment_se = sqrt((observed_prop/expected_mean)^2 * expected_var * (1 / observed_prop^2 + 1 / (length(expected_prop) * expected_mean^2)))
  )
}


# Example usage -----------------------------------------------------------
# Load required data:

# enrichment_results <- calculate_chromatin_enrichment(
#   index_snps = index_variants,
#   control_snps = matched_controls,
#   chromatin_annotations = chromhmm_annotations,
#   n_cores = 8,
#   n_iterations = 1000

# )

# 1. Load required data
args = commandArgs(TRUE)
tissue = args[1]
ancestry = args[2]
n_cores = as.numeric(args[3])
print(paste0("Start loading data for ", tissue))

# 2. Precompute
snpeff_path = paste0("/path/to/scRNA/",tissue,"/sc-eQTL/results/functional_annotation/snpEff/")
chromhmm_path = paste0("/path/to/scRNA/",tissue,"/sc-eQTL/results/functional_annotation/chromHMM/") 
if (!file.exists(paste0(chromhmm_path,"precomputed_matching_controls_",ancestry,"_5en8.rds"))){
  print(paste0("Start loading data for ", tissue))
  # Load gene annotations
  gene_annotations = fread(paste0("/path/to/scRNA/Blood/resource/genes.bed"), select=c(1,3,4))
  setnames(gene_annotations, c("end","gene_id"), c("TSS","phenotype_id"))
  gene_annotations = gene_annotations[chr %in% paste("chr",1:22, sep="")]

  # Load background SNPs (all tested vcf)
  background_path = paste0("/path/to/scRNA/",tissue,"/sc-eQTL/genotype/",tissue,"_",ancestry,"_indivpruned.frq")
  background_snps = fread(background_path, select=c(2,5))
  colnames(background_snps) = c("variant_id","afreq")
  background_snps[,c("chr","pos") := tstrsplit(variant_id,"_",keep=1:2)]
  background_snps[, afreq := ifelse(afreq>0.5,1-afreq,afreq)]
  background_snps[,chr := paste0("chr",chr)]

  # Load eSNPs
  index_snps = fread(paste0(chromhmm_path,tissue,"_",ancestry,"_all_esnp_5en8.bed"), select=c(1,3,4))
  colnames(index_snps) = c("chr","pos","variant_id")
  index_snps$pos = as.character(index_snps$pos)
  index_snps = background_snps[index_snps, .(variant_id,chr,pos,afreq), on=c("variant_id","chr","pos")]

  print("Precomputing matching control SNPs")
  start.time <- Sys.time()
  precomputed <- preprocess_controls(
    index_snps = index_snps,
    gene_annotations = gene_annotations,
    background_snps = background_snps,
    n_cores = n_cores
  )
  end.time <- Sys.time()
  end.time-start.time
  saveRDS(precomputed, paste0(chromhmm_path,"precomputed_matching_controls_",ancestry,"_5en8.rds"))
  rm(background_snps)
} else{
  precomputed <- readRDS(paste0(chromhmm_path,"precomputed_matching_controls_",ancestry,"_5en8.rds"))
}

# 3. Load SNP annotations
annotation_list = read.table(file.path(snpeff_path,"annotation.list"))$V1

overlap_files = list.files(snpeff_path, pattern='overlaps')

a <- annotation_list[as.numeric(args[4])]
print(a)

fs = overlap_files[str_detect(overlap_files,ancestry) & str_detect(overlap_files,a) & str_detect(overlap_files, "null")]
chromatin_annotation <- lapply(fs, function(f) {
    df=fread(file.path(snpeff_path, f),select=1, header=FALSE)
    colnames(df) = "variant_id"
    df[,Annotation := a]
  }) %>% rbindlist()

# 4. Perform functional enrichment for each annotation
func_enrich <- calculate_chromatin_enrichment_optimized(
    precomputed = precomputed,
    chromatin_annotation = chromatin_annotation,
    n_iterations = 1000,
    n_cores =  n_cores # Already parallelized at higher level
)
print("Writing out data")
write.table(func_enrich, paste0("/path/to/scRNA/crossTissue_analysis/functional_enrichment_5en8/functional_enrichment_snpEff_5en8_",tissue,"_",ancestry,"_", a, ".txt"), row=F,sep="\t",col=T,quo=F)

