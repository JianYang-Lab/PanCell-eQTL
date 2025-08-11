##############################################################################
# Create eQTL query file from parquet for SMR and HEIDI analysis
##############################################################################


import pandas as pd
import gzip
import os
import sys
import gc

# Inputs
tissue = sys.argv[1]
ancestry = sys.argv[2]
celltype = sys.argv[3]
pc_num = 5

dir_local=f"/data/cc_qr_{tissue}_{ancestry}_{celltype}"
os.system(f"mkdir -p {dir_local}")
os.chdir(dir_local)

def load_genes_gtf(annotation_gtf="/path/to/Blood/resource/genes.gtf"): 
    os.system(f"cp {annotation_gtf} {dir_local}")
    gene_gtf = pd.read_csv(annotation_gtf, sep='\t')
    gene_gtf.columns = ['gene_name','phenotype_id','gene_chr','gene_start','gene_end','strand'] 
    gene_gtf['gene_tss'] = gene_gtf.apply(lambda row: row['gene_start'] if row['strand']=="+" else row['gene_end'], axis=1)
    gene_gtf = gene_gtf[['phenotype_id', 'gene_name', 'gene_tss', 'strand']]
    return(gene_gtf)

def load_parquet_chr(celltype, ancestry, pc_num, chrom, gtf, parquet_dir=f"/path/to/{tissue}/sc-eQTL/results/raw"):
    parquet_file = f'{celltype}_{ancestry}_pc{pc_num}.cis_qtl_pairs.{chrom}.parquet'
    os.system(f"cp {parquet_dir}/{parquet_file} {dir_local}")
    df = pd.read_parquet(parquet_file)
    df = df[(df['af']>=0.01) & (df['af']<=0.99)]
    df = df.merge(gtf, on='phenotype_id',how='left')
    df[['CHR', 'POS', 'REF', 'ALT']] = df['variant_id'].str.split("_", expand=True).iloc[:, :4]
    print(f'Chr{chrom} processed') 
    return(df)

def add_rsid_with_flip(df,chrom, variant_files_dir="/storage/yangjianLab/sharedata/GATK_resource/dbSNP_b151/split_into_chrAuto"):
    # Process each chromosome present in the GWAS data
    variant_file = f"chr{chrom}_variants_hg38_b151.txt"
    os.system(f"cp {variant_files_dir}/{variant_file} {dir_local}")
    print(f"Adding RSID for chromosome {chrom} using file {variant_file}...")
    # Load the variant file for this chromosome
    try:
        variants = pd.read_csv(variant_file, sep="\t", header=None, names=['CHR', 'POS', 'RSID', 'REF', 'ALT'], 
                                usecols=[0, 1, 2, 3, 4], dtype={'CHR': str, 'POS': str, 'RSID': str, 'REF': str, 'ALT': str})
        variants = variants.assign(ALT=variants['ALT'].str.split(',')).explode('ALT').reset_index(drop=True)
    except FileNotFoundError:
        print(f"File {variant_file} not found. Skipping chromosome {chrom}.") 
        return(df)
    # Memory-efficient flipped REF/ALT
    variants_flipped = variants.copy(deep=False)
    variants_flipped.rename(columns={'REF': 'ALT', 'ALT': 'REF'}, inplace=True)
    # Merge original REF/ALT
    merged = pd.merge(df, variants, on=['CHR', 'POS', 'REF', 'ALT'], how='left')
    del variants  # Free memory
    gc.collect()
    # Merge flipped REF/ALT
    merged_flipped = pd.merge(df, variants_flipped, on=['CHR', 'POS', 'REF', 'ALT'], how='left')
    del variants_flipped  # Free memory
    gc.collect()
    # Combine both merges
    merged['RSID'] = merged['RSID'].combine_first(merged_flipped['RSID'])
    del merged_flipped  # Free memory
    gc.collect()
    print(f"Chromosome {chrom} RSID added")
    return merged

def process_chromosome(chrom, celltype, ancestry, pc_num, gtf):
    df = load_parquet_chr(celltype, ancestry, pc_num, chrom, gtf)
    df = add_rsid_with_flip(df, chrom)
    df['RSID'] = df['RSID'].fillna(df['variant_id'])
    df = df[['RSID', 'CHR', "POS", "ALT", "REF", "af", "phenotype_id", "CHR", "gene_tss", "gene_name", "strand", "slope", "slope_se", "pval_nominal"]]
    df.columns = ["SNP", "Chr", "BP", "A1", "A2", "Freq", "Probe", "Probe_Chr", "Probe_bp", "Gene", "Orientation", "b", "se", "p"]
    # Filter by Freq
    df = df[(df.Freq>=0.01) & (df.Freq<=0.99)]
    # Handle duplciate SNPs for each probe (SNPs with different ALT alleles)
    new_snp = df['Chr'].astype(str) + "_" + df['BP'].astype(str) + "_" + df['A2'] + "_" + df['A1'] + "_b38"
    duplicated_mask = df.duplicated(subset=['Probe', 'SNP'], keep=False)
    df.loc[duplicated_mask, 'SNP'] = new_snp[duplicated_mask]
    # Handle SNPs with multiple positions (Normally DELINs that are merged into other RSID)
    snp_with_multiple_pos = df.groupby('SNP')['BP'].nunique()
    duplicates = snp_with_multiple_pos[snp_with_multiple_pos > 1].index
    duplicated_mask = df.SNP.isin(duplicates)
    df.loc[duplicated_mask, 'SNP'] = new_snp[duplicated_mask]  
    return df

# Load gene gtf file
gtf = load_genes_gtf()

# Sequentially process chromosomes to save memory
output_file = f"{tissue}_{celltype}_{ancestry}_sceQTL_query.txt"
with open(output_file, 'w') as f:
    for i in range(1, 23):
        print(f"Processing chromosome {i}...")
        df = process_chromosome(i, celltype, ancestry, pc_num, gtf)
        # Append to the output file incrementally
        df.to_csv(f, sep='\t', index=False, header=f.tell() == 0, na_rep="NA")  # Add header only for the first chunk
        # Free memory
        del df
        gc.collect()

os.makedirs(f"/path/to/{tissue}/sc-eQTL/results/query", exist_ok=True)
os.system(f"cp {output_file} /path/to/{tissue}/sc-eQTL/results/query/")
os.system("cd /data")
os.system(f"rm -rf {dir_local}")

print(f"Results saved to {output_file}")
