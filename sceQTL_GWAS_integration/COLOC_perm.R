#!/usr/bin/env Rscript
##############################################################################
# Script information
# Colocalization analysis of eQTL and GWAS summary statistics
##############################################################################

suppressMessages(library("optparse"))
option_list = list(
  make_option("--sumstats", action="store", default=NA, type='character',
              help="Path to summary statistics [required]"),
  make_option("--GWASN", action="store", default=NA, type='integer',
              help="Sample size of GWAS [optional]"),
  make_option("--besd", action="store", default=NA, type='character',
              help="Path to QTL summary data in BESD format [required]"),
  make_option("--QTLName", action="store", default=NA, type='character',
              help="Name of QTL summary data [required]"),
  make_option("--QTLN", action="store", default=NA, type='integer',
              help="Sample size of QTL summary data [required]"),
  make_option("--expr_bed", action="store", default=NA, type='character',
              help="Path to gene expression bed file"),
  make_option("--probe_path", action="store", default=NA, type='character',
              help="Path to probes to be analyzed, specify when --perm=TRUE"),
  make_option("--perm", action="store", default=FALSE, type='logical',
              help="Whether to use permutation probes"),
  make_option("--chr", action="store", default=NA, type='character',
              help="chromosome"),
  make_option("--out", action="store", default=NA, type='character',
              help="Path to result file"),
  make_option("--SMR", action="store", default=NA, type='character',
              help="Path to SMR software [required]"),
  make_option("--R_functions", action="store", default="CommonFunc.r", type='character',
              help="Path to R_functions [required]"),
  make_option("--tmp_dirt", action="store", default=NA, type='character',
              help="Temporary directory of cis and query results [required]")
)
opt = parse_args(OptionParser(option_list=option_list))

suppressMessages({
  library("coloc");
  library(dplyr);
  library(data.table);
  library(parallel)
})
source(opt$R_functions) 
if (!file.exists(opt$R_functions)) {
  stop("R_functions file does not exist. Please provide a valid path.")
}
# Check required options
if (is.na(opt$sumstats) || is.na(opt$besd) || is.na(opt$QTLName) || is.na(opt$QTLN) || is.na(opt$expr_bed) || is.na(opt$SMR) || is.na(opt$tmp_dirt)) {
  stop("Please provide all required options: --sumstats, --besd, --QTLName, --QTLN, --expr_bed, --SMR, --tmp_dirt")
}   
if (!file.exists(opt$sumstats)) {
  stop("Summary statistics file does not exist. Please provide a valid path.")
}
if (!file.exists(opt$besd)) {
  stop("QTL summary data file does not exist. Please provide a valid path.")
}
if (!file.exists(opt$expr_bed)) {
  stop("Expression bed file does not exist. Please provide a valid path.")
}
if (!file.exists(opt$SMR)) {
  stop("SMR software does not exist. Please provide a valid path.")
}       
    
SMR = opt$SMR
tmp_dirt = opt$tmp_dirt
dir.create(paste0(tmp_dirt,"/query"), recursive=TRUE, showWarnings=FALSE)

# ------------------ GWAS Input Treatment ---------------------
gwas = fread(opt$sumstats, head=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
colnames(gwas) = c("SNP","A1","A2","Freq","b","se","P","N")
rm_index = which(gwas$se == 0 | is.na(gwas$se))
if(length(rm_index) > 0){
  gwas = gwas[-rm_index,]
}

if (sum(is.na(gwas$N)) > 0){
  if (is.na(opt$GWASN)){
    stop("Please specify GWAS sample size since there is missing data in the summary statistics.")
  } else {
    gwas[is.na(gwas$N),"N"] = opt$GWASN # If N is NA, replace it by the input
  }
}

# ------------------ QTL Input Treatment ---------------------
if (!opt$perm){
  command1 = paste0(SMR, " --beqtl-summary ", opt$besd, " --chr ", opt$chr, " --descriptive-cis --out ", tmp_dirt, "/query/", opt$QTLName, "_chr", opt$chr)
  system(command1)
  qtl = fread(paste0(tmp_dirt, "/query/", opt$QTLName, "_chr", opt$chr, ".cis.summary.txt"), head=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
  qtl = qtl[qtl$Probe_Chr == opt$chr, ]
  if (!is.null(opt$probe_path) && file.exists(opt$probe_path)) {
    pbs = fread(opt$probe_path, header=FALSE, stringsAsFactors=FALSE)$V1
  } else {
    pbs = unique(qtl$Probe)
  }
} else {
  snp_pbs = fread(opt$probe_path, header=FALSE, stringsAsFactors=FALSE, data.table=FALSE)
  pbs = unique(snp_pbs[,2])
}

bed = fread(opt$expr_bed, header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
bed = bed[bed$TargetID %in% pbs, ]

# ------------------ Coloc Analysis in Parallel ------------------------

process_probe <- function(probe) {
  tryCatch({
    # Query SNPs based on the probe
    query_output = paste0(tmp_dirt, "/query/", opt$QTLName, "_", probe, ".query")
    command = paste0(
      SMR, " --beqtl-summary ", opt$besd,
      " --query 1 --probe ", probe,
      " --out ", query_output
    )
    system(command)
    
    query_file = paste0(query_output, ".txt")
    if (!file.exists(query_file)) {
      message("Query file not found for probe ", probe, ". Skipping.")
      return(NULL)
    }
    
    # Read and filter query results
    query = fread(query_file, head=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
    query = query[query$SE != 0, ]
    if (nrow(query) == 0) return(NULL)
    
    # Extract gene name
    gene_name = query$Gene[1]
    
    # Match query SNPs with GWAS data
    matched_indices = match(query$SNP, gwas$SNP)
    valid_indices = which(!is.na(matched_indices))
    if (length(valid_indices) < 10) return(NULL)
    
    query_com = query[valid_indices, ]
    gwas_com = gwas[matched_indices[valid_indices], ]
    
    # Harmonize allele coding between GWAS and query data
    gwas_com$A1 = toupper(gwas_com$A1)
    gwas_com$A2 = toupper(gwas_com$A2)
    query_com$A1 = toupper(query_com$A1)
    query_com$A2 = toupper(query_com$A2)
    
    # Flip alleles where necessary
    # Flip matching where alleles are reversed
    reverse_indices = which((query_com$A1 == gwas_com$A2) & (query_com$A2 == gwas_com$A1))
    if (length(reverse_indices) > 0){
      gwas_com[reverse_indices, c("A1", "A2")] = query_com[reverse_indices, c("A1", "A2")]
      gwas_com[reverse_indices, "b"] = -gwas_com[reverse_indices, "b"]  # Flip effect size
      if (!all(is.na(gwas_com$Freq))) {
        gwas_com[reverse_indices, "Freq"] = 1 - gwas_com[reverse_indices, "Freq"]  # Flip frequency
      }
    }

    # Remove SNPs where alleles do not match
    allele_match = which((query_com$A1 == gwas_com$A1) & (query_com$A2 == gwas_com$A2))
    if (length(allele_match) < 10) return(NULL)
    
    query_com = query_com[allele_match, ]
    gwas_com = gwas_com[allele_match, ]
    
    snpbuf = query_com$SNP
    
    # Prepare GWAS summary statistics
    freq1 = as.numeric(gwas_com$Freq)
    pval1 = as.numeric(gwas_com$P)
    n1 = gwas_com$N
    bzy_hat = gwas_com$b
    bzy_se = gwas_com$se
    
    # Check and fix invalid frequencies
    invalid_freq1 = which(is.na(freq1) | freq1 <= 0 | freq1 >= 1)
    if (length(invalid_freq1) > 0) {
      freq1[invalid_freq1] = query_com$Freq[invalid_freq1]
      std_eff1 = calcu_std_b_se(bzy_hat[invalid_freq1] / bzy_se[invalid_freq1], freq1[invalid_freq1], n1[invalid_freq1])
      bzy_hat[invalid_freq1] = std_eff1[, 1]
      bzy_se[invalid_freq1] = std_eff1[, 2]
    }
    
    data1 = data.frame(
      pvalues = pval1,
      N = n1,
      MAF = freq1,
      beta = bzy_hat,
      varbeta = bzy_se^2,
      snp = snpbuf,
      stringsAsFactors = FALSE
    )
    
    # Prepare QTL summary statistics
    pval2 = query_com$p
    n2 = opt$QTLN
    freq2 = query_com$Freq
    bzx_hat = query_com$b
    bzx_se = query_com$SE

    # Check and fix invalid frequencies
    invalid_freq2 = which(is.na(freq2) | freq2 <= 0 | freq2 >=1)
    if (length(invalid_freq2) > 0) {
      freq2[invalid_freq2] = freq1[invalid_freq2]
      z2 = bzx_hat[invalid_freq2] / bzx_se[invalid_freq2]
      std_eff2 = calcu_std_b_se(z2, freq2[invalid_freq2], n2)
      bzx_hat[invalid_freq2] = std_eff2[, 1]
      bzx_se[invalid_freq2] = std_eff2[, 2]
    }
    
    # sdY calculation
    bed_row = bed[bed$TargetID == probe, -c(1:4)]
    sdY = sd(as.numeric(bed_row))
    
    data2 = data.frame(
      pvalues = pval2,
      N = n2,
      MAF = freq2,
      beta = bzx_hat,
      varbeta = bzx_se^2,
      snp = snpbuf,
      stringsAsFactors = FALSE
    )
    
    # Filter valid MAF and SNP data
    valid_indices = which(
      !is.na(data1$MAF) & !is.na(data2$MAF) &
      data1$MAF > 0 & data2$MAF > 0 &
      data1$MAF < 1 & data2$MAF < 1
    )
    data1 = data1[valid_indices, ]
    data2 = data2[valid_indices, ]

    if (nrow(data1) < 10) return(NULL)
    data1 = as.list(data1)
    data2 = as.list(data2)
    data1$sdY = sdY
    data1$type = "cc"
    data2$type = "quant"
   
    # Perform coloc analysis
    resbuf = coloc.abf(data1,data2)
    
    # Extract coloc results
    pp_values = resbuf$summary
    pp_snps = resbuf$results
    #topsnp = pp_snps[which.max(pp_snps$SNP.PP.H4),]$snp
    #topsnp_pp = pp_snps[which.max(pp_snps$SNP.PP.H4),]$SNP.PP.H4
    
    pp0 = as.numeric(pp_values[2])
    pp1 = as.numeric(pp_values[3])
    pp2 = as.numeric(pp_values[4])
    pp3 = as.numeric(pp_values[5])
    pp4 = as.numeric(pp_values[6])
    #PP3_PP4 = pp3 + pp4
    
    res = data.frame(probe_name = probe, gene_name = gene_name, pp0, pp1, pp2, pp3, pp4, stringsAsFactors = FALSE)
    
    # Clean up temporary files
    file.remove(query_file)
    
    return(res)
  }, error = function(e) {
    message("Error processing probe ", probe, ": ", conditionMessage(e))
    return(NULL)
  })
}

# Set the number of cores for parallel processing
numCores =  6  # Leave one core for the OS
if (numCores < 1) numCores = 1  # Ensure at least one core is used

# Run the analysis in parallel
results_list = mclapply(pbs, process_probe, mc.cores = numCores)

# Combine results while removing NULLs
results_list = results_list[!sapply(results_list, is.null)]
summary = do.call(rbind, results_list)

# Append or write results
output_file = paste0(opt$out, ".coloc")
write.table(summary, output_file, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
