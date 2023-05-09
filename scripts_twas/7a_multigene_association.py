'''
For each regional phenotype, rank genes by 
association p-value and then compute cumulative 
set associations. 

Repeat on permutations in order to compute the 
p-value of the observed average correlations. 

- Nhung, updated Feb 2023  
'''

import sys 
import numpy as np  
import h5py 
from statsmodels.stats.multitest import fdrcorrection 
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
base_path = '/data1/rubinov_lab/brain_genomics/scripts_assoc_clean'

phen_dir = base_path + '/inputs_HCP/phen_regress'
expr_dir = base_path + '/inputs_HCP/expr_regress'

prm_file = base_path + '/inputs_HCP/permutations_100k.hdf5' 
rho_file = base_path + '/outputs_HCP/assoc_1M'

out_path = base_path + '/outputs_HCP/multigene_ranked_sets'

## permutation indexing 
with h5py.File(prm_file, 'r') as f: 
    samp_order = f['samp_idx'][()] ## (subjs,) -- twin-sorted order  
    subj_order = f['subj_idx'][()] ## (perms, subjs) -- original unsorted order 
    nsubjs = samp_order.size 

## phenotype (original unsorted order) 
phens = {} ## k: (reg, phen), v: subject array 
for phen in phenotypes: 
    with h5py.File('{}/{}.hdf5'.format(phen_dir, phen), 'r') as f: 
        for reg in regions: 
            phens[(reg,phen)] = f[reg][()] ## (subjs,)      

## expression (twin-sorted order)  
genes = {} ## k: reg, v: gene array  
exprs = {} ## k: reg, v: (genes, subjs) expression matrix  
for reg in regions:
    with h5py.File('{}/{}.hdf5'.format(expr_dir, reg), 'r') as f: 
        genes[reg] = f['genes'][()]
        exprs[reg] = f['pred_expr'][()][:,samp_order] ## (genes, subjs)

## function: vectorized correlation 
def compute_pearson(expr_mat, phen_mat):
    c_expr = expr_mat - expr_mat.mean(axis=1)[:,None] 
    c_phen = phen_mat - phen_mat.mean(axis=1)[:,None] 
    s_expr = (c_expr**2).sum(axis=1)[None] 
    s_phen = (c_phen**2).sum(axis=1)[:,None] 
    num = np.dot(c_phen, c_expr.T) 
    dom = np.sqrt(np.dot(s_phen, s_expr))
    return num / dom 

## function: rank genes and compute cumulative associations 
def compute_set_assocs(phen, reg): 

    ## read single-gene association results 
    ofile = '{}/pvals_{}/{}.hdf5'.format(rho_file, phen, reg)
    with h5py.File(ofile, 'r') as f: 
        obsv_assoc = f['pearson'][()] ## (genes, [rho pval fdr]) 
    obsv_rhos = obsv_assoc[:,0][None] ## (1, genes) 
    obsv_pval = obsv_assoc[:,1][None] ## (1, genes) 
    num_sig = (obsv_pval <= 0.005).sum()

    nfile = '{}/nulls/pvals_{}/{}.hdf5'.format(rho_file, phen, reg)
    with h5py.File(nfile, 'r') as f: 
        null_assoc = f['null_pearson'][()] ## (perms, genes, [rho pval fdr]) 

    null_rhos = null_assoc[:,:,0] ## (perms, genes) 
    null_pval = null_assoc[:,:,1] ## (perms, genes) 

    ## rank genes/expression based on p-value 
    obsv_idx = np.argsort(obsv_pval) ## (1, genes) 
    null_idx = np.argsort(null_pval, axis=1) ## (perms, genes) 

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

    ## compute average gene-set rho 
    obsv_rho_avg = np.mean(c_obsv_rhos) ## (1,)
    null_rho_avg = np.mean(c_null_rhos, axis=0) ## (perms,)     
    pval_rho_avg = np.mean(null_rho_avg >= obsv_rho_avg) ## (1,)

    ## save results for plotting 
    out_file = '{}/{}_{}.hdf5'.format(out_path, phen, reg)
    with h5py.File(out_file, 'w') as f: 
        f['obsv_rhos'] = c_obsv_rhos ## (genes, 1) 
        f['null_rhos'] = c_null_rhos ## (genes, perms) 
        f['p0.005'] = num_sig ## (1,)

        f['avg_rho_pval'] = pval_rho_avg ## (1,)
        f['obsv_avg_rho'] = obsv_rho_avg ## (1,)  
        f['null_avg_rhos'] = null_rho_avg ## (perms,) 

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


