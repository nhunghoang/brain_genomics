'''
Single-gene stats for polygenic comparison. 

Run TWAS to get the r^2 values (p-value should 
match the PrediXcan Association version). Per 
regional phenotype, store matrix with dimensions 

(genes, [r^2, pval, FDR threshold, Bonferroni threshold])

The thresholds should be the same across genes of a 
given regional phenotype. 

- Nhung, Oct 2023
'''

import numpy as np
import pandas as pd
import h5py
import sys

from statsmodels.regression.linear_model import OLS as OLS
from statsmodels.tools.tools import add_constant
from statsmodels.stats.multitest import multipletests as sm

## params 
group = sys.argv[1] ## HCP, HCP/nonTwin, ...

regions = ['dlpfc', 'anterior-cingulate', 'amygdala', 'hippocampus', \
           'caudate', 'putamen', 'nucleus-accumbens', 'cerebellar-hemisphere']
phens = ['vol_mean', 'alff_mean', 'reho_noGS_mean', 'connmean_noGS_mean']
reg_phens = [(r,p) for r in regions for p in phens]

## paths
main_path = '/data1/rubinov_lab/brain_genomics/paper_twas'
twas_path = f'{main_path}/outputs_{group}/twas_JTI' ## /{phen}/{reg}.txt
grex_path = f'{main_path}/inputs_{group}/grex_JTI' ## /{reg}.hdf5
phen_path = f'{main_path}/inputs_{group}/phenotypes' ## /{phen}_regr.csv

outs_path = f'{main_path}/outputs_{group}/polygenic_models/single_stats.hdf5'

## load grex and phen data
grex_all = {} ## k: (reg, phen), v: (subj * gene)
phen_all = {} ## k: (reg, phen), v: (subj,)

## init outputs
sin_data = {} ## k: (reg, phen), v: ([gene], [twas pval])
sin_sigs = {} ## k: (reg, phen), v: {ben, bon}

for (reg, phen) in reg_phens: 

    ## get significant TWAS genes 
    tfile = f'{twas_path}/{phen}/{reg}.txt'
    df = pd.read_table(tfile, sep='\t', usecols=['gene', 'pvalue'])

    bon_keep, bons, _, bon_line = sm(df['pvalue'], method='bonferroni', alpha=0.05)
    ben_keep, bens, _, _ = sm(df['pvalue'], method='fdr_bh', alpha=0.05)

    sin_sigs[(reg,phen)] = {'bon_line': bon_line, 'ben_line': df['pvalue'][ben_keep].max()}

    df['bon'] = bons
    df['ben'] = bens

    keep = (df['pvalue'] <= 0.001)
    df_keep = df.loc[keep]
    twas_genes = df_keep['gene'].values
    twas_pvals = df_keep['pvalue'].values

    twas_bons = df_keep['bon'].values
    twas_bens = df_keep['ben'].values

    ## get grex matrices and slice 
    gpath = f'{grex_path}/{reg}.hdf5' 
    with h5py.File(gpath, 'r') as f:
        reg_gene = f['genes'][()].astype(str)
        grex_mat = f['pred_expr'][()]

    gset, idx1, idx2 = np.intersect1d(reg_gene, twas_genes, return_indices=True)
    grex_all[(reg, phen)] = grex_mat[idx1].T
    sin_data[(reg, phen)] = (gset, {'unc': twas_pvals[idx2], \
                                    'bon': twas_bons[idx2], \
                                    'ben': twas_bens[idx2]}) 

    ## get phenotype arrays 
    ppath = f'{phen_path}/{phen}_regr.csv' 
    df = pd.read_table(ppath, sep='\t', usecols=[reg])
    phen_all[(reg,phen)] = df[reg].values

## func: fit and report single-gene regression 
def regression(reg, phen, gidx):

    ## single-gene data 
    gg = grex_all[(reg,phen)][:,gidx] ## (subjs,)
    pp = phen_all[(reg,phen)] ## (subjs,)

    ## for OLS: add intercept constant
    gg = add_constant(gg)

    ## fit regression model
    model = OLS(pp, gg).fit()
    r2 = model.rsquared
    pv = model.pvalues[1]

    #print('tw: {}\nnh: {}\n'.format(,pv))
    n = sin_data[(reg,phen)][1]['unc'][gidx]
    h = pv 
    assert(n.round(6) == h.round(6))

    return r2, pv

## main loop 
sin_stats = {} ## k: (reg,phen), v: genes * [r2, pv, ben, bon]
for (reg, phen) in reg_phens: 

    nsig = sin_data[(reg,phen)][0].size
    r2s = np.zeros((nsig,1))
    pvs = np.zeros((nsig,1))
    bhs = np.zeros((nsig,1))
    bfs = np.zeros((nsig,1))

    for i in range(nsig): 
        r2, pv = regression(reg, phen, i)
        r2s[i][0] = r2
        pvs[i][0] = pv
        bhs[i][0] = sin_data[(reg,phen)][1]['ben'][i]
        bfs[i][0] = sin_data[(reg,phen)][1]['bon'][i]

    results = np.concatenate([r2s, pvs, bhs, bfs], axis=1)
    sin_stats[(reg,phen)] = results 

with h5py.File(outs_path, 'w') as f: 
    for (reg, phen), stats in sin_stats.items(): 
        key = '{}X{}'.format(reg,phen)
        f[key] = stats


