'''
'''

import numpy as np 
import pandas as pd 
import h5py 

from scipy.stats import pearsonr
from sklearn.linear_model import LinearRegression

## params 
group = 'HCP/nonTwin' #sys.argv[1] ## HCP, HCP/nonTwin, ...
num_perms = 10000

regions = ['dlpfc', 'anterior-cingulate', 'amygdala', 'hippocampus', \
           'caudate', 'putamen', 'nucleus-accumbens', 'cerebellar-hemisphere']
phens = ['vol_mean', 'alff_mean', 'reho_noGS_mean', 'connmean_noGS_mean']
reg_phens = [(r,p) for r in regions for p in phens]

## paths 
main_path = '/data1/rubinov_lab/brain_genomics/paper_twas'

grex_path = f'{main_path}/inputs_{group}/grex_JTI'
phen_path = f'{main_path}/inputs_{group}/phenotypes'
covs_path = f'{main_path}/inputs_{group}/covariates.csv'
obsv_path = f'{main_path}/outputs_{group}/poly_stats_observed.hdf5'

twa0_path = f'{main_path}/outputs_{group}/twas_JTI/perms'
twa1_path = '/data/rubinov_lab/brain_genomics_project/polygenic/twas'

outs_path = f'{main_path}/outputs_{group}/poly_stats_perms.hdf5'

############################################################

## load number of genes from observed 
num_genes = {} ## k: (reg, phen), v: (1,) 
with h5py.File(obsv_path, 'r') as f:
    for (reg, phen) in reg_phens: 
        num_genes[(reg,phen)] = f[f'{reg}X{phen}Xgene'][()]

## load covariates 
covs_mat = pd.read_table(covs_path) \
             .drop(columns=['subject', 'sample']) \
             .values

## load grex 
gene_arr = {} ## k: reg, v: (gene,)
grex_mat = {} ## k: reg, v: (subj,gene)
for reg in regions: 
    gfile = f'{grex_path}/{reg}.hdf5' 
    with h5py.File(gfile, 'r') as f: 
        gene_arr[reg]= f['genes'][()].astype(str)
        grex_mat[reg] = f['pred_expr'][()].T 

###########################################################################

## func: compute (gene sum, phen resid) correlation
def compute_corr(p, phen_arr, reg, phen):

    ## load twas data
    if p < 3000: tfile = f'{twa0_path}/{phen}/{reg}_{p}.txt' 
    else: tfile = f'{twa1_path}/{phen}/{reg}_{p}.txt' 

    tdata = pd.read_table(tfile, usecols=['gene', 'pvalue'])

    ## take top N (from observed) genes 
    ngene = num_genes[(reg,phen)] 
    tgene = tdata[:ngene]['gene'].values

    ## slice grex accordingly 
    gmask = np.isin(gene_arr[reg], tgene) 
    assert(gmask.sum() > 0, f'{tfile} format is wrong')

    grex_twa = grex_mat[reg][:,gmask] 
    
    ## flip grex signs based on corr 
    phen_norm = (phen_arr - phen_arr.mean()) / phen_arr.std() 
    grex_norm = (grex_twa - grex_twa.mean(axis=0)) / grex_twa.std(axis=0) 
    r_pearson = np.dot(phen_norm, grex_norm) / phen_norm.size     

    signs = np.sign(r_pearson) 
    grex_sign = grex_norm * signs 

    ## sum grex and correlate with phen 
    grex_sum = np.mean(grex_sign, axis=1) 
    rho, pval = pearsonr(phen_arr, grex_sum) 
    return rho

###########################################################################

## main loop 
perm_corrs = {} ## k: (reg, phen), v: (perms,)
for phen in phens: 
    for reg in regions: 
        pcorrs = np.zeros(num_perms)

        ## load first 5k phen data 
        cols = [f'{reg}_{p}' for p in range(5000)]
        ppath0 = f'{phen_path}/perm_{phen}.csv'
        pdata0 = pd.read_table(ppath0, usecols=cols) 
        print('read pdata 0')

        for p in range(5000): 
            
            ## compute phen residuals 
            phen_old = pdata0[f'{reg}_{p}'].values
            lr = LinearRegression().fit(covs_mat, phen_old)
            ypred = lr.predict(covs_mat) 
            phen_arr = phen_old - ypred

            ## compute (gene sum, phen resid) corrs 
            pcorrs[p] = compute_corr(p, phen_arr, reg, phen)
            if (p%500) == 0: 
                txt = f'{p}/5000 {reg} {phen} ' + str(round(pcorrs[p], 3))
                print(txt)
        
        ## load second 5k phen data 
        cols = [f'{reg}_{p}' for p in range(5000, num_perms)]
        ppath1 = f'{phen_path}/perm_{phen}2.csv'
        pdata1 = pd.read_table(ppath1, usecols=cols) 
        print('read pdata 1')

        for p in np.arange(5000) + 5000: 
            
            ## compute phen residuals 
            phen_old = pdata1[f'{reg}_{p}'].values
            lr = LinearRegression().fit(covs_mat, phen_old)
            ypred = lr.predict(covs_mat) 
            phen_arr = phen_old - ypred

            ## compute (gene sum, phen resid) corrs 
            pcorrs[p] = compute_corr(p, phen_arr, reg, phen)
            if (p%500) == 0: 
                txt = f'{p}/{num_perms} {reg} {phen} ' + str(round(pcorrs[p], 3))
                print(txt)

        perm_corrs[(reg,phen)] = pcorrs

## write to file 
with h5py.File(outs_path, 'w') as f: 
    for (reg, phen) in reg_phens: 
        f[f'{reg}X{phen}'] = perm_corrs[(reg,phen)] 
