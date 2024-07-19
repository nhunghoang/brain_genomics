'''
'''

import numpy as np 
import pandas as pd 
import h5py 

from scipy.stats import pearsonr
from statsmodels.stats.multitest import multipletests as sm

## params 
group = 'HCP/nonTwin' #sys.argv[1] ## HCP, HCP/nonTwin, ...

regions = ['dlpfc', 'anterior-cingulate', 'amygdala', 'hippocampus', \
           'caudate', 'putamen', 'nucleus-accumbens', 'cerebellar-hemisphere']
phens = ['vol_mean', 'alff_mean', 'reho_noGS_mean', 'connmean_noGS_mean']
reg_phens = [(r,p) for r in regions for p in phens]

## paths 
main_path = '/data1/rubinov_lab/brain_genomics/paper_twas'

grex_path = f'{main_path}/inputs_{group}/grex_JTI'
phen_path = f'{main_path}/inputs_{group}/phenotypes'
twas_path = f'{main_path}/outputs_{group}/twas_JTI'
perm_path = f'{main_path}/outputs_{group}/poly_stats_perms.hdf5'

outs_path = f'{main_path}/outputs_{group}/poly_stats_observed.hdf5'

## load grex and phen data
gene_all = {} ## k: reg, v: (gene,)
grex_all = {} ## k: reg, v: (subj * gene)
phen_all = {} ## k: (reg, phen), v: (subj,)

for reg in regions:
    gpath = f'{grex_path}/{reg}.hdf5'
    with h5py.File(gpath, 'r') as f:
        reg_gene = f['genes'][()].astype(str)
        grex_mat = f['pred_expr'][()]

        gene_all[reg] = reg_gene
        grex_all[reg] = grex_mat.T

for phen in phens:
    ppath = f'{phen_path}/{phen}_regr.csv'
    df = pd.read_table(ppath, sep='\t', usecols=regions)
    for reg in regions:
        phen_all[(reg,phen)] = df[reg].values

## load perm corrs 
perm_corrs = {} ## k: (reg, phen), v: (perms,)
with h5py.File(perm_path, 'r') as f: 
    for key in f.keys(): 
        [reg, phen] = key.split('X')
        vals = f[key][()]
        perm_corrs[(reg,phen)] = np.abs(vals)  

############################################################

grex_sums = {} ## k: (reg, phen), v: (subj,)
phen_arrs = {} ## k: (reg, phen), v: (subj,)
corr_arrs = {} ## k: (reg, phen), v: ([pearson, pval])
num_genes = {} ## k: (reg, phen), v: (1,)

## strongest correlation between phen and a gene 
twa_corrs = {} ## k: (reg, phen), v: (1,) 

## k: (reg, phen), v: [best single-gene FDR, poly-gene FDR] 
fval_arrs = {(reg,phen): np.zeros(2) for (reg, phen) in reg_phens} 
pval_arrs = {(reg,phen): np.zeros(2) for (reg, phen) in reg_phens} 

for (reg, phen) in reg_phens: 

    ## load grex and phen data 
    grex_mat = grex_all[reg]
    phen_arr = phen_all[(reg,phen)] 

    ## select TWAS genes 
    tfile = f'{twas_path}/{phen}/{reg}.txt' 
    tdata = pd.read_table(tfile, usecols=['gene', 'pvalue'])
    tgene = tdata.loc[tdata['pvalue'] < 0.001]['gene'].values  

    ## store the best TWAS FDR  
    fval = sm(tdata['pvalue'], method='fdr_bh', alpha=0.05) 
    fval_arrs[(reg,phen)][0] = fval[1].min()
    pval_arrs[(reg,phen)][0] = tdata['pvalue'].min()

    ## slice grex accordingly 
    gmask = np.isin(gene_all[reg], tgene)
    grex_twa = grex_mat[:,gmask]

    ## flip grex signs based on corr 
    phen_norm = (phen_arr - phen_arr.mean()) / phen_arr.std() 
    grex_norm = (grex_twa - grex_twa.mean(axis=0)) / grex_twa.std(axis=0) 
    r_pearson = np.dot(phen_norm, grex_norm) / phen_norm.size     

    twa_corrs[(reg,phen)] = r_pearson.max()

    signs = np.sign(r_pearson) 
    grex_sign = grex_norm * signs 

    ## sum grex and correlate with phen 
    grex_sum = np.mean(grex_sign, axis=1) 
    rho, pval = pearsonr(phen_arr, grex_sum) 

    ## compute pval 
    pval = np.mean(perm_corrs[(reg,phen)] >= rho)
    if pval == 0: pval = np.array(1/10000)
    pval_arrs[(reg,phen)][1] = pval

    fval_arrs[(reg,phen)][1] = pval

    ## store
    grex_sums[(reg,phen)] = grex_sum 
    phen_arrs[(reg,phen)] = phen_arr
    corr_arrs[(reg,phen)] = [rho, pval]
    num_genes[(reg,phen)] = tgene.size 

    print(tgene.size, rho.round(4), pval.round(4), reg, phen)

## save outputs 
with h5py.File(outs_path, 'w') as f: 
    for (reg, phen) in reg_phens: 
        f[f'{reg}X{phen}Xgrex'] = grex_sums[(reg,phen)]
        f[f'{reg}X{phen}Xphen'] = phen_arrs[(reg,phen)]
        f[f'{reg}X{phen}Xcorr'] = corr_arrs[(reg,phen)]
        f[f'{reg}X{phen}Xgene'] = num_genes[(reg,phen)]

        f[f'{reg}X{phen}XpFDR'] = fval_arrs[(reg,phen)]
        f[f'{reg}X{phen}XpNOM'] = pval_arrs[(reg,phen)]
        f[f'{reg}X{phen}XrTWA'] = twa_corrs[(reg,phen)]

