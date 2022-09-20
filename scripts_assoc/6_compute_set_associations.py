'''
For each regional phenotype, sort genes by 
correlation and then compute cumulative set 
associations. 

- Nhung, updated Sept 2022 
'''

import sys 
import os 
import numpy as np  
import h5py 
from statsmodels.stats.multitest import fdrcorrection 
from sklearn.metrics import auc as sk_AUC
from multiprocessing import Pool 
from time import time 
import warnings 

warnings.filterwarnings('ignore', 'invalid value encountered in true_divide')

## regional phenotypes 
phenotypes = ['gm_volume', 'myelination', 'alff', 'reho_noGS', 'connmean_noGS']
regions = ['hippocampus', 'amygdala', 'hypothalamus', 'substantia-nigra',\
        'caudate', 'putamen', 'nucleus-accumbens', 'anterior-cingulate',\
        'frontal-pole', 'cerebellar-hemisphere']

## paths 
phen_dir = '/data1/rubinov_lab/brain_genomics/analyses_HCP/DATA_OUTPUT/phen_regress'
expr_dir = '/data1/rubinov_lab/brain_genomics/analyses_HCP/DATA_OUTPUT/r0.3_p0.01/expr_regress'
rho_file = '/data1/rubinov_lab/brain_genomics/analyses_HCP/assoc_1M'
prm_file = '/data1/rubinov_lab/brain_genomics/analyses_HCP/DATA_OUTPUT/null_permutations_latest_100k.hdf5' 

main_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/assoc_1M/set_data'
if not os.path.exists(main_path): os.mkdir(main_path) 

## permutation data 
with h5py.File(prm_file, 'r') as f: 
    samp_order = np.array(f['samp_idx']) ## (subjs,) -- twin-sorted 
    subj_order = np.array(f['subj_idx']) ## (perms, subjs) 
    nsubjs = samp_order.size 

## phenotype (original order) 
phens = {} ## k: (reg, phen), v: subject array 
for phen in phenotypes: 
    with h5py.File('{}/{}.hdf5'.format(phen_dir, phen), 'r') as f: 
        for reg in regions: 
            phens[(reg,phen)] = np.array(f[reg]) ## (subjs,)      

## expression (twin-sorted order)  
genes = {}; exprs = {} ## k: reg, v: gene array / (genes, subjs) expression matrix  
for reg in regions:
    with h5py.File('{}/{}.hdf5'.format(expr_dir, reg), 'r') as f: 
        genes[reg] = np.array([g.decode("utf-8").split('.')[0] for g in np.array(f['genes'])])
        exprs[reg] = np.array(f['pred_expr'])[:,samp_order] ## (genes, subjs)

## function: vectorized correlation 
def compute_pearson(expr_mat, phen_mat):
    c_expr = expr_mat - expr_mat.mean(axis=1)[:,None] 
    c_phen = phen_mat - phen_mat.mean(axis=1)[:,None] 
    s_expr = (c_expr**2).sum(axis=1)[None] 
    s_phen = (c_phen**2).sum(axis=1)[:,None] 
    num = np.dot(c_phen, c_expr.T) 
    dom = np.sqrt(np.dot(s_phen, s_expr))
    return num / dom 

## function: compute AUC significance
def compute_auc(obsv_rhos, null_rhos):
    idxs = np.arange(obsv_rhos.size)
    nperms = null_rhos.shape[1]

    obsv_auc = sk_AUC(idxs, obsv_rhos)
    null_auc = np.empty(nperms)
    for i in range(nperms):
        null_auc[i] = sk_AUC(idxs, null_rhos[:,i])

    auc_pval = (null_auc >= obsv_auc).mean()
    return obsv_auc, null_auc, auc_pval

## function: sort genes and compute cumulative associations 
def compute_set_assocs(phen, reg): 

    ## read single-gene association results 
    ofile = '{}/pvals_{}/{}_1M.hdf5'.format(rho_file, phen, reg)
    with h5py.File(ofile, 'r') as f: 
        obsv_assoc = np.array(f['pearson']) ## (genes, [rho pval fdr]) 
    obsv_rhos = obsv_assoc[:,0][None] ## (1, genes) 
    obsv_pval = obsv_assoc[:,1][None] ## (1, genes) 
    num_sig = (obsv_pval <= 0.005).sum()

    nfile = '{}/nulls/pvals_{}/{}.hdf5'.format(rho_file, phen, reg)
    with h5py.File(nfile, 'r') as f: 
        null_assoc = np.array(f['null_pearson']) ## (perms, genes, [rho pval fdr]) 

    null_rhos = null_assoc[:,:,0] ## (perms, genes) 
    null_pval = null_assoc[:,:,1] ## (perms, genes) 

    ## sort genes/expression based on correlation 
    obsv_idx = np.argsort(obsv_rhos) ## (1, genes) 
    null_idx = np.argsort(null_rhos, axis=1) ## (perms, genes) 

    obsv_expr = exprs[reg][obsv_idx] ## (1, sorted genes, subjs) 
    null_expr = exprs[reg][null_idx] ## (perms, sorted genes, subjs) 

    ## change expression direction based on rho 
    obsv_rhos_sorted = np.take_along_axis(obsv_rhos, obsv_idx, 1) ## (1, sorted genes)
    null_rhos_sorted = np.take_along_axis(null_rhos, null_idx, 1) ## (perms, sorted genes) 

    obsv_rho_signs = np.where(obsv_rhos_sorted < 0, -1, 1) ## (1, sorted genes)
    null_rho_signs = np.where(null_rhos_sorted < 0, -1, 1) ## (perms, sorted genes)
    
    obsv_expr *= obsv_rho_signs[:,:,None] ## (1, genes, subjs) 
    null_expr *= null_rho_signs[:,:,None] ## (perms, genes, subjs) 

    ## compute cumulative expression 
    obsv_expr_sums = np.cumsum(obsv_expr, axis=1) ## (1, genes, subjs) 
    null_expr_sums = np.cumsum(null_expr, axis=1) ## (perms, genes, subjs) 

    ## compute pearsons 
    obsv_phen = phens[(reg,phen)][samp_order][None] ## (1, subjs) 
    null_phen = phens[(reg,phen)][subj_order] ## (perms, subjs)

    num_genes = null_assoc.shape[1]
    set_perms = null_assoc.shape[0]

    c_obsv_rhos = compute_pearson(obsv_expr_sums[0], obsv_phen).T ## (genes, 1)
    c_null_rhos = np.empty((num_genes, set_perms)) ## (genes, perms)
    for p in range(set_perms):
        c_null_rhos[:,p] = compute_pearson(null_expr_sums[p], null_phen[p][None])

    ## compute AUC significance 
    obsv_auc, null_auc, auc_pval = compute_auc(c_obsv_rhos, c_null_rhos)

    ## save pearsons, num sig genes, and auc pval for plotting
    out_path = '{}/{}_{}.hdf5'.format(main_path, phen, reg)
    with h5py.File(out_path, 'w') as f: 
        f['obsv_rhos'] = c_obsv_rhos ## (genes, 1) 
        f['null_rhos'] = c_null_rhos ## (genes, perms) 
        f['p0.005'] = num_sig ## (1,)
        f['auc_p'] = auc_pval ## (1,)
        f['obsv_auc'] = obsv_auc ## (1,)  
        f['null_auc'] = null_auc ## (perms,) 

##############################################################################
##############################################################################

start0 = time()
for phen in phenotypes: 
    for reg in regions: 
        start = time()
        compute_set_assocs(phen, reg)

        secs = time() - start 
        hr = int(secs//3600)
        mn = int((secs%3600)//60)
        sc = int((secs%3600)%60)
        print('{:d} hr, {:d} mn, {:d} sc | {} {}'.format(hr, mn, sc, phen, reg))

secs = time() - start0
hr = int(secs//3600)
mn = int((secs%3600)//60)
sc = int((secs%3600)%60)
print('total runtime: {:d} hr, {:d} mn, {:d} sc'.format(hr, mn, sc))


