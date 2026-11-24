"""PanCell-eQTL conditional mapping extension for tensorQTL."""

import time
from collections import OrderedDict

import numpy as np
import pandas as pd
import torch
from tensorqtl import genotypeio
from tensorqtl.cis import calculate_cis_nominal, map_cis, map_independent, map_nominal
from tensorqtl.core import Residualizer, SimpleLogger, calculate_maf, get_t_pval, impute_mean

# Expose the stock cis methods used by the driver alongside the study method.
__all__ = ["map_cis", "map_nominal", "map_independent", "map_independent_5en8"]


def _cis_windows(generator, variant_index, verbose):
    """Yield cis windows with variant IDs on the study and stock tensorQTL readers."""
    if hasattr(generator, "generate_data_new"):
        yield from generator.generate_data_new(verbose=verbose)
        return

    for phenotype, genotypes, indices, phenotype_id in generator.generate_data(verbose=verbose):
        yield phenotype, variant_index[indices], genotypes, indices, phenotype_id


def map_independent_5en8(
    genotype_df, variant_df, cis_df, phenotype_df, phenotype_pos_df, covariates_df=None,
    group_s=None, maf_threshold=0, window=1000000, missing=-9,
    random_tiebreak=False, logger=None, logp=False, verbose=True,
):
    """Map conditionally independent cis signals using nominal P < 5e-8.

    cis_df is the per-phenotype output from tensorQTL map_cis. The study uses
    individual phenotypes, untransformed p-values and the first minimum in a tie.
    """
    if group_s is not None or logp or random_tiebreak:
        raise ValueError("The 5e-8 study analysis requires ungrouped phenotypes, "
                         "raw p-values and deterministic tie-breaking.")
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    if logger is None:
        logger = SimpleLogger()

    logger.write('cis-QTL mapping: conditionally independent variants')
    logger.write(f'  * {phenotype_df.shape[1]} samples')

    # subset significant phenotypes
    signif_df = cis_df[cis_df['pval_nominal']<5e-8].copy()
    if len(signif_df) == 0:
        raise ValueError(f"No significant phenotypes at p-value<5e-8.")
    cols = [ 'variant_id','start_distance', 'end_distance', 'ma_samples', 'ma_count', 'af',
            'pval_nominal', 'slope', 'slope_se']
    signif_df = signif_df[cols]
    ix = phenotype_df.index[phenotype_df.index.isin(signif_df.index)]
    logger.write(f'  * {signif_df.shape[0]}/{cis_df.shape[0]} significant phenotypes')

    assert phenotype_df.index.equals(phenotype_pos_df.index)

    if covariates_df is not None:
        assert covariates_df.index.equals(phenotype_df.columns), 'Sample names in phenotype matrix columns and covariate matrix rows do not match!'
        assert ~(covariates_df.isnull().any().any()), f'Missing or null values in covariates matrix, in columns {",".join(covariates_df.columns[covariates_df.isnull().any(axis=0)].astype(str))}'
        logger.write(f'  * {covariates_df.shape[1]} covariates')

    logger.write(f'  * {genotype_df.shape[0]} variants')
    if maf_threshold > 0:
        logger.write(f'  * applying in-sample {maf_threshold} MAF filter')
    logger.write(f'  * cis-window: ±{window:,}')
    phenotype_df = phenotype_df.loc[ix]
    phenotype_pos_df = phenotype_pos_df.loc[ix]

    genotype_ix = np.array([genotype_df.columns.tolist().index(i) for i in phenotype_df.columns])
    genotype_ix_t = torch.from_numpy(genotype_ix).to(device)
    ix_dict = {i:k for k,i in enumerate(genotype_df.index)}

    res_df = []
    igc = genotypeio.InputGeneratorCis(genotype_df, variant_df, phenotype_df, phenotype_pos_df, group_s=group_s, window=window)
    if igc.n_phenotypes == 0:
        raise ValueError('No valid phenotypes found.')
    logger.write('  * computing independent QTLs')
    start_time = time.time()

    for phenotype, variant_ids, genotypes, _, phenotype_id in _cis_windows(
        igc, genotype_df.index, verbose
    ):
        # Align the genotypes to the phenotype samples.
        phenotype_t = torch.tensor(phenotype, dtype=torch.float).to(device)
        genotypes_t = torch.tensor(genotypes, dtype=torch.float).to(device)
        genotypes_t = genotypes_t[:,genotype_ix_t]
        impute_mean(genotypes_t)

        if maf_threshold > 0:
            maf_t = calculate_maf(genotypes_t)
            mask_t = maf_t >= maf_threshold
            genotypes_t = genotypes_t[mask_t]
            mask = mask_t.cpu().numpy().astype(bool)
            variant_ids = variant_ids[mask]

        # 1) forward pass
        forward_df = [signif_df.loc[phenotype_id]]  # initialize results with top variant
        # initialize covariates
        if covariates_df is not None:
            covariates = covariates_df.values.copy()
        else:
            covariates = np.empty((phenotype_df.shape[1], 0))
        dosage_dict = {}
        while True:
            # add variant to covariates
            variant_id = forward_df[-1]['variant_id']
            ig = genotype_df.values[ix_dict[variant_id], genotype_ix].copy()
            m = ig == missing
            ig[m] = ig[~m].mean()
            dosage_dict[variant_id] = ig
            covariates = np.hstack([covariates, ig.reshape(-1,1)]).astype(np.float32)
            dof = phenotype_df.shape[1] - 2 - covariates.shape[1]
            residualizer = Residualizer(torch.tensor(covariates, dtype=torch.float32).to(device))

            res = calculate_cis_nominal(genotypes_t, phenotype_t, residualizer=residualizer)
            tstat, slope, slope_se, af, ma_samples, ma_count = [i.cpu().numpy() for i in res]
            df = pd.DataFrame({
                'pval_nominal' : get_t_pval(tstat, dof, log=logp),
                'slope': slope, 'slope_se': slope_se,
                'tstat': tstat, 'af': af, 'ma_samples': ma_samples, 'ma_count': ma_count,
            }, index=variant_ids)

            variant_id = df.pval_nominal.idxmin()
            if df.pval_nominal.min() < 5e-8:
                start_distance = variant_df.loc[variant_id,'pos']- igc.phenotype_start[phenotype_id]
                end_distance =  variant_df.loc[variant_id,'pos'] - igc.phenotype_end[phenotype_id]
                res_s = pd.Series(OrderedDict([
                    ('variant_id', variant_id),
                    ('start_distance',start_distance),
                    ('end_distance',end_distance),
                    ('ma_samples',np.nan),
                    ('ma_count',np.nan),
                    ('af',np.nan),
                    ('pval_nominal',np.nan),
                    ('slope',np.nan),
                    ('slope_se',np.nan)
                ]), name=phenotype_id)
                res_s[['ma_samples','ma_count','af','pval_nominal','slope','slope_se']] = df.loc[variant_id,['ma_samples','ma_count','af','pval_nominal','slope','slope_se']]

                forward_df.append(res_s)
            else:
                break
        forward_df = pd.concat(forward_df, axis=1, sort=False).T
        dosage_df = pd.DataFrame(dosage_dict)

        # 2) backward pass
        if forward_df.shape[0] > 1:
            back_df = []
            variant_set = set()
            for k,i in enumerate(forward_df['variant_id'], 1):
                if covariates_df is not None:
                    covariates = np.hstack([
                        covariates_df.values,
                        dosage_df[np.setdiff1d(forward_df['variant_id'], i)].values,
                    ])
                else:
                    covariates = dosage_df[np.setdiff1d(forward_df['variant_id'], i)].values
                dof = phenotype_df.shape[1] - 2 - covariates.shape[1]
                residualizer = Residualizer(torch.tensor(covariates, dtype=torch.float32).to(device))

                res = calculate_cis_nominal(genotypes_t, phenotype_t, residualizer=residualizer)
                tstat, slope, slope_se, af, ma_samples, ma_count = [i.cpu().numpy() for i in res]
                df = pd.DataFrame({
                    'pval_nominal' : get_t_pval(tstat, dof, log=logp),
                    'slope': slope, 'slope_se': slope_se,
                    'tstat': tstat, 'af': af, 'ma_samples': ma_samples, 'ma_count': ma_count,
                }, index=variant_ids)
                variant_id = df.pval_nominal.idxmin()
                if df.loc[variant_id, 'pval_nominal'] < 5e-8 and variant_id not in variant_set:
                    start_distance = variant_df.loc[variant_id,'pos']- igc.phenotype_start[phenotype_id]
                    end_distance =  variant_df.loc[variant_id,'pos'] - igc.phenotype_end[phenotype_id]
                    res_s = pd.Series(OrderedDict([
                        ('variant_id', variant_id),
                        ('start_distance',start_distance),
                        ('end_distance',end_distance),
                        ('ma_samples',np.nan),
                        ('ma_count',np.nan),
                        ('af',np.nan),
                        ('pval_nominal',np.nan),
                        ('slope',np.nan),
                        ('slope_se',np.nan)
                    ]), name=phenotype_id)
                    res_s[['ma_samples','ma_count','af','pval_nominal','slope','slope_se']] = df.loc[variant_id,['ma_samples','ma_count','af','pval_nominal','slope','slope_se']]
                    res_s['rank'] = k
                    back_df.append(res_s)
                    variant_set.add(variant_id)
            if len(back_df)>0:
                res_df.append(pd.concat(back_df, axis=1, sort=False).T)
        else:  # single independent variant
            forward_df['rank'] = 1
            res_df.append(forward_df)

    res_df = pd.concat(res_df, axis=0, sort=False)
    res_df.index.name = 'phenotype_id'
    logger.write(f'  Time elapsed: {(time.time()-start_time)/60:.2f} min')
    logger.write('done.')
    return res_df.reset_index()


# Adapted from tensorQTL (https://github.com/broadinstitute/tensorqtl).
# BSD 3-Clause License
#
# Copyright (c) 2018-2019, The Broad Institute, Inc. and The General Hospital Corporation.
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
#
# * Neither the name of the copyright holder nor the names of its
#   contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
