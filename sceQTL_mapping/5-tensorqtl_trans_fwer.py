##############################################################################
# Step 5: Empirical trans-eQTL FWER threshold from joint donor permutations
# Usage: python 5-tensorqtl_trans_fwer.py EUR CD4TNCInf 5 1000 \
#          --base Blood/sc-eQTL_R --gene-mappability resources/hg38_gene_mappability.txt.gz
# Dependencies: numpy, pandas, scipy, torch, tensorqtl.
# Covariates: PF*, PC1..PC{pc_num}, Condition_group*, Age_pred, Gender_pred.
# Uses all SuSiE credible-set variants.
# Each permutation scans all eligible trans SNP-gene pairs and retains one
# minimum nominal P. Their 5th percentile is the global nominal P cutoff.
# Exclusion criteria: >=5-cell RINT, MAF>=0.01, MAC>=10 and +/-5Mb.
##############################################################################

import argparse
import hashlib
import json
import os
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
import torch
import tensorqtl
from tensorqtl import core

# Input loading, genotype QC, and scanning helpers are included here so this
# file does not depend on other project scripts. No beta fits or q-values.
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")


def log(message):
    print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] '+message, flush=True)


def load_phenotype_and_covariates(ancestry, celltype, pc_num, mappability_file, cell_counts_file):
    expression_bed = f"celltype_matrix/bed/{celltype}_{ancestry}_expr.bed"
    covariates_file = f"covariates/final/{celltype}_{ancestry}_cov_final.txt"
    phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)
    assert phenotype_df.index.is_unique and phenotype_df.columns.is_unique
    # Keep the existing RINT matrix intact; fail if its donors violate the QC.
    cell_cov = pd.read_csv(cell_counts_file, index_col='project_sample')
    if not cell_cov.loc[phenotype_df.columns, 'cell_count'].ge(5).all():
        raise ValueError('The existing RINT matrix includes donors with fewer than 5 cells.')

    covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0)
    pf_cols = list(covariates_df.columns[covariates_df.columns.str.startswith('PF')])
    pc_cols = [f'PC{i}' for i in range(1, pc_num+1)]
    cd_cols = list(covariates_df.columns[covariates_df.columns.str.startswith('Condition_group')])
    other_cols = ['Age_pred', 'Gender_pred']
    covariates_df = covariates_df.loc[phenotype_df.columns, pf_cols+pc_cols+cd_cols+other_cols]
    # Centered one-hot condition columns sum to zero; omit one reference level.
    drop_cols = []
    if len(cd_cols) > 1 and np.allclose(covariates_df[cd_cols].sum(axis=1), 1):
        drop_cols.append(cd_cols[0])
    drop_cols += [c for c in covariates_df if covariates_df[c].nunique() == 1 and c not in drop_cols]
    covariates_df = covariates_df.drop(columns=drop_cols)
    cov = covariates_df.to_numpy()
    if not np.isfinite(cov).all() or np.linalg.matrix_rank(cov-cov.mean(axis=0)) != cov.shape[1]:
        raise ValueError('Covariates must be finite and have full centered rank.')
    if len(cov)-cov.shape[1]-2 <= 0:
        raise ValueError('Insufficient residual degrees of freedom.')
    if drop_cols:
        log('Dropped redundant/constant covariates: '+', '.join(drop_cols))

    gene_mappability = pd.read_csv(mappability_file, sep='\t', header=None,
                                  names=['phenotype_id', 'mappability'])
    gene_mappability['phenotype_id'] = gene_mappability['phenotype_id'].str.split('.').str[0]
    gene_keep = gene_mappability.loc[gene_mappability['mappability'] >= 0.8, 'phenotype_id']
    phenotype_df = phenotype_df.loc[phenotype_df.index.isin(gene_keep)]
    if not np.isfinite(phenotype_df.to_numpy()).all():
        raise ValueError('Expression contains missing or nonfinite values.')
    phenotype_df = phenotype_df.loc[np.ptp(phenotype_df.to_numpy(), axis=1) > 0]
    phenotype_pos_df = phenotype_pos_df.loc[phenotype_df.index].copy()
    phenotype_pos_df['chr'] = phenotype_pos_df['chr'].astype(str).str.replace(r'^chr', '', regex=True)
    if phenotype_df.empty:
        raise ValueError('No nonconstant, mappable genes remain.')
    return(phenotype_df, phenotype_pos_df, covariates_df)


def filter_genotypes(genotype_df, variant_df, covariates_df, batch_size=10000):
    """Apply MAF/MAC to called alleles, before mean imputation can inflate MAC."""
    genotype_df = genotype_df.loc[:, covariates_df.index]
    keep = np.zeros(len(genotype_df), dtype=bool)
    af, mac = np.zeros(len(genotype_df)), np.zeros(len(genotype_df))
    residualizer = core.Residualizer(torch.tensor(covariates_df.to_numpy(), dtype=torch.float32, device=device))
    for start in range(0, len(genotype_df), batch_size):
        stop = min(start+batch_size, len(genotype_df))
        g = torch.tensor(genotype_df.iloc[start:stop].to_numpy(), dtype=torch.float32, device=device)
        if not torch.isfinite(g).all() or not ((g == -9) | ((g >= 0) & (g <= 2))).all():
            raise ValueError('Expected diploid genotypes in [0,2], with -9 for missing calls.')
        called = g != -9
        allele_count = torch.where(called, g, 0).sum(dim=1)
        allele_number = 2*called.sum(dim=1)
        batch_af = allele_count/allele_number.clamp_min(1)
        batch_mac = torch.minimum(allele_count, allele_number-allele_count)
        passed = (batch_mac/allele_number.clamp_min(1) >= 0.01) & (batch_mac >= 10)
        ix = torch.where(passed)[0]
        if len(ix):
            subset = g[ix].clone()
            core.impute_mean(subset)
            # Exclude constant genotypes, including those explained by covariates.
            valid = residualizer.transform(subset).var(dim=1) > 1e-8
            keep[start+ix[valid].cpu().numpy()] = True
        af[start:stop], mac[start:stop] = batch_af.cpu().numpy(), batch_mac.cpu().numpy()
    variant_df = variant_df.loc[genotype_df.index].copy()
    variant_df['af'], variant_df['mac'] = af, mac
    if not keep.any():
        raise ValueError('No variants pass MAF>=0.01, MAC>=10 and residual-variance QC.')
    log(f'{keep.sum():,}/{len(keep):,} variants passed MAF>=0.01, MAC>=10 and residual-variance QC')
    return(genotype_df.loc[keep], variant_df.loc[keep])


def load_genotype(ancestry, covariates_df, batch_size=10000, prefix=None):
    """Decode selected donors from SNP-major BED in small read-only batches."""
    prefix = str(prefix) if prefix is not None else f'genotype/Blood_{ancestry}_indivpruned_allcissig_mappability_filtered'
    variant = pd.read_csv(prefix+'.bim',sep=r'\s+',header=None,
        names=['chrom','variant_id','cm','pos','a1','a2'],dtype={'chrom':str}).set_index('variant_id')
    family = pd.read_csv(prefix+'.fam',sep=r'\s+',header=None,dtype=str)[1]
    index = pd.Index(family).get_indexer(covariates_df.index)
    if not family.is_unique or (index < 0).any() or not variant.index.is_unique:
        raise ValueError('Genotype variants/donors must be unique and include every expression donor.')
    width = (len(family)+3)//4
    with open(prefix+'.bed','rb') as handle:
        if handle.read(3) != b'\x6c\x1b\x01' or os.path.getsize(prefix+'.bed') != 3+len(variant)*width:
            raise ValueError('Expected a complete SNP-major PLINK BED file.')
    packed = np.memmap(prefix+'.bed',mode='r',dtype=np.uint8,offset=3,shape=(len(variant),width))
    values = np.empty((len(variant),len(index)),dtype=np.int8)
    # Codes 00/01/10/11 encode A1 dosage 2/missing/1/0, as in TensorQTL.
    lookup = np.array([2,-9,1,0],dtype=np.int8)
    for start in range(0,len(variant),batch_size):
        code = (packed[start:start+batch_size][:,index//4] >> (2*(index%4))) & 3
        values[start:start+len(code)] = lookup[code]
    variant['chrom'] = variant.chrom.str.replace(r'^chr','',regex=True)
    genotype = pd.DataFrame(values,index=variant.index,columns=covariates_df.index)
    return filter_genotypes(genotype,variant[['chrom','pos']],covariates_df,batch_size)


class TransScanner:
    """Cache normalized residual genotypes once; scan batched phenotype rows."""
    def __init__(self, genotype_df, variant_df, positions, covariates_df,
                 variant_batch=10000, phenotype_batch=2048, window=5000000):
        self.device = device
        self.variant_batch, self.phenotype_batch, self.window = variant_batch, phenotype_batch, window
        self.residualizer = core.Residualizer(torch.tensor(covariates_df.to_numpy(), dtype=torch.float32, device=self.device))
        self.dof = len(covariates_df)-len(covariates_df.columns)-2
        chroms = pd.Index(pd.concat([variant_df.chrom, positions['chr'].astype(str)]).unique())
        self.vchr = torch.tensor(chroms.get_indexer(variant_df.chrom), device=self.device)
        self.pchr = torch.tensor(chroms.get_indexer(positions['chr'].astype(str)), device=self.device)
        self.vpos = torch.tensor(variant_df.pos.to_numpy(), device=self.device)
        self.left = torch.tensor(positions['pos' if 'pos' in positions else 'start'].to_numpy(), device=self.device)
        self.right = torch.tensor(positions['pos' if 'pos' in positions else 'end'].to_numpy(), device=self.device)
        size = genotype_df.shape[0]*genotype_df.shape[1]*4
        cache_device = self.device
        if self.device.type == 'cuda' and size > .6*torch.cuda.mem_get_info()[0]:
            cache_device = torch.device('cpu')
        self.g = torch.empty(genotype_df.shape, dtype=torch.float32, device=cache_device)
        for start in range(0, len(genotype_df), variant_batch):
            g = torch.tensor(genotype_df.iloc[start:start+variant_batch].to_numpy(), dtype=torch.float32, device=self.device)
            core.impute_mean(g)
            normalized = core.center_normalize(self.residualizer.transform(g), dim=1)
            if not torch.isfinite(normalized).all():
                raise ValueError('Invalid normalized genotype after production QC.')
            self.g[start:start+len(g)] = normalized.to(cache_device)
        log(f'Cached {len(self.g):,} residual genotype rows on {cache_device}')

    def residuals(self, values):
        x = torch.tensor(values, dtype=torch.float32, device=self.device)
        return self.residualizer.transform(x).cpu().numpy()

    def scan(self, values, genes):
        """Return full-panel and +/-5Mb-excluded maximum R2 for every row."""
        genes = torch.as_tensor(genes, dtype=torch.long, device=self.device)
        if len(values) != len(genes):
            raise ValueError('Every phenotype row needs its target gene index.')
        full = torch.zeros(len(values), device=self.device)
        eligible = torch.full((len(values),), -1., device=self.device)
        for left in range(0, len(values), self.phenotype_batch):
            right = min(left+self.phenotype_batch, len(values))
            p = torch.tensor(values[left:right], dtype=torch.float32, device=self.device)
            p = self.residualizer.transform(p)
            norm = p.square().sum(dim=1, keepdim=True).sqrt()
            p = torch.where(norm > 1e-8, p/norm.clamp_min(1e-8), torch.zeros_like(p))
            gi = genes[left:right]
            for start in range(0, len(self.g), self.variant_batch):
                stop = min(start+self.variant_batch, len(self.g))
                r2 = (self.g[start:stop].to(self.device) @ p.T).square_().clamp_(0, 1)
                full[left:right] = torch.maximum(full[left:right], r2.amax(dim=0))
                cis = ((self.vchr[start:stop, None] == self.pchr[gi]) &
                       (self.vpos[start:stop, None] >= self.left[gi]-self.window) &
                       (self.vpos[start:stop, None] <= self.right[gi]+self.window))
                r2.masked_fill_(cis, -1)
                eligible[left:right] = torch.maximum(eligible[left:right], r2.amax(dim=0))
        if (eligible < 0).any():
            raise ValueError('A target gene has no eligible trans SNPs.')
        return full.cpu().numpy(), eligible.cpu().numpy()


def nominal_pvalue(r2, dof):
    """Two-sided regression P, using the same residual df for every pair."""
    r2 = np.clip(np.asarray(r2, dtype=float), 0, 1)
    with np.errstate(divide='ignore'):
        return 2 * stats.t.sf(np.sqrt(dof * r2 / (1 - r2)), dof)


def permutation_minima(scanner, values, nperms=1000, seed=20261112, batch_size=4,
                       method='rint'):
    """Jointly shuffle raw RINT or nuisance residuals, then reproject covariates."""
    if method == 'freedman-lane':
        # Freedman-Lane: Y* = fitted(Y | covariates) + permuted residuals.
        # scan() projects out the same covariates, removing the fitted term
        # exactly. Only the residuals need to be shuffled and passed to scan().
        values = scanner.residuals(values)
    elif method != 'rint':
        raise ValueError(f'Unknown permutation method: {method}')
    rng = np.random.RandomState(seed)
    gene_index = np.arange(len(values))
    maxima = np.empty(nperms, dtype=float)
    for start in range(0, nperms, batch_size):
        stop = min(start + batch_size, nperms)
        # One common donor order per permutation preserves gene dependence.
        shuffled = [values[:, rng.permutation(values.shape[1])] for _ in range(stop - start)]
        _, trans_r2 = scanner.scan(np.concatenate(shuffled), np.tile(gene_index, len(shuffled)))
        # With common df, max R2 across all genes/SNPs gives minimum nominal P.
        maxima[start:stop] = trans_r2.reshape(len(shuffled), len(values)).max(axis=1)
        log(f'Joint trans permutations: {stop}/{nperms}')
    if not np.isfinite(maxima).all() or ((maxima < 0) | (maxima > 1)).any():
        raise ValueError('Invalid permutation maximum R2.')
    return pd.DataFrame(dict(permutation=np.arange(1, nperms + 1),
                             max_r2=maxima, min_pval=nominal_pvalue(maxima, scanner.dof)))


def global_cutoff(minima, alpha=.05):
    """Quantile is taken on nominal P values, after reducing each permutation."""
    p = np.asarray(minima, dtype=float)
    if p.ndim != 1 or len(p) < 2 or not np.isfinite(p).all() or ((p < 0) | (p > 1)).any():
        raise ValueError('Expected one finite minimum P per permutation.')
    if not 0 < alpha < 1:
        raise ValueError('alpha must be within (0,1).')
    return float(np.quantile(p, alpha, method='linear'))


def main():
    parser = argparse.ArgumentParser(description='Global trans-eQTL cutoff from permutation minimum nominal P values.')
    parser.add_argument('ancestry')
    parser.add_argument('celltype')
    parser.add_argument('pc_num', type=int)
    parser.add_argument('nperms', nargs='?', type=int, default=1000)
    parser.add_argument('--alpha', type=float, default=.05)
    parser.add_argument('--seed', type=int, default=20261112)
    parser.add_argument('--permutation-method', choices=['rint', 'freedman-lane'], default='rint',
                        help='Shuffle raw RINT (default) or covariate-model residuals (Freedman-Lane).')
    parser.add_argument('--permutation-batch', type=int, default=4)
    parser.add_argument('--variant-batch', type=int, default=10000)
    parser.add_argument('--phenotype-batch', type=int, default=2048)
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--base', type=Path, default=Path(os.environ.get('TRANS_EQTL_BASE', '.')),
                        help='Analysis directory containing genotype, expression, and covariate inputs.')
    parser.add_argument('--gene-mappability', type=Path, required=True,
                        help='Headerless gene-ID/mappability TSV (optionally gzipped).')
    parser.add_argument('--cell-counts', type=Path,
                        help='CSV with project_sample and cell_count; defaults to the original relative layout.')
    parser.add_argument('--bfile', type=Path,
                        help='PLINK prefix; default: the Blood ancestry-specific SuSiE credible-set panel.')
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    if args.pc_num < 0 or min(args.permutation_batch, args.variant_batch, args.phenotype_batch, args.threads) < 1:
        parser.error('Batch/thread counts must be positive; pc_num must be nonnegative.')
    if args.nperms < 20 or not 0 < args.alpha < 1:
        parser.error('Use at least 20 permutations and alpha within (0,1).')
    args.base = args.base.resolve()
    args.gene_mappability = args.gene_mappability.resolve()
    args.cell_counts = (args.cell_counts or args.base /
        f'Blood_rawQC_{args.celltype}_{args.ancestry}_cov.csv').resolve()
    for path in [args.gene_mappability, args.cell_counts]:
        if not path.is_file():
            parser.error(f'Missing input file: {path}')
    genotype_prefix = (args.bfile or args.base /
        f'genotype/Blood_{args.ancestry}_indivpruned_allcissig_mappability_filtered_susie_cs').resolve()
    for extension in ['bed', 'bim', 'fam']:
        if not Path(str(genotype_prefix) + '.' + extension).is_file():
            parser.error(f'Missing genotype file: {genotype_prefix}.{extension}')
    suffix = '_freedman_lane' if args.permutation_method == 'freedman-lane' else ''
    output = (args.out or args.base / 'results/trans_eqtl/fwer_validation/susie_cs' /
              f'{args.celltype}_{args.ancestry}_pc{args.pc_num}_B{args.nperms}{suffix}').resolve()
    if any((output / name).exists() for name in ['summary.json', 'permutation_minima.tsv']):
        parser.error('Output already contains results; choose another --out directory.')
    output.mkdir(parents=True, exist_ok=True)
    torch.set_num_threads(args.threads)
    os.chdir(args.base)

    phenotype, positions, covariates = load_phenotype_and_covariates(args.ancestry, args.celltype, args.pc_num,
        args.gene_mappability, args.cell_counts)
    log(f'Genotype panel: {genotype_prefix}')
    genotype, variants = load_genotype(args.ancestry, covariates, args.variant_batch,
                                       prefix=genotype_prefix)
    scanner = TransScanner(genotype, variants, positions, covariates,
                           args.variant_batch, args.phenotype_batch)
    del genotype
    log(f'Permutation method: {args.permutation_method}')
    minima = permutation_minima(scanner, phenotype.to_numpy(dtype=np.float32),
                                args.nperms, args.seed, args.permutation_batch, args.permutation_method)
    cutoff = global_cutoff(minima.min_pval, args.alpha)
    minima.to_csv(output / 'permutation_minima.tsv', sep='\t', index=False)
    log(f'{args.celltype}: global nominal P cutoff = {cutoff:.8g} (alpha={args.alpha}, B={args.nperms})')


if __name__ == '__main__':
    main()
