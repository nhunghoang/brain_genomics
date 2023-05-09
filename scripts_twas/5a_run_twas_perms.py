'''
Script for running permutations of single-gene associations. 
Each regional phenotype has 10k permutations, where the p-value 
of each permutation is based on 1M permutations.  

- Nhung, updated Sept 2022
'''

import sys 
import os 
import numpy as np  
import h5py 
from statsmodels.stats.multitest import fdrcorrection 
from multiprocessing import Pool  
from time import time 
import warnings 

warnings.filterwarnings('ignore', 'invalid value encountered in true_divide')

dset = sys.argv[1] ## HCP or UKB   
assoc_perms = int(1e6) 
set_perms = int(1e4) 
rng = np.random.default_rng() ## random number generator used for nontwin permutations 

regions = ['hippocampus', 'amygdala', 'hypothalamus', 'substantia-nigra',\
        'caudate', 'putamen', 'nucleus-accumbens', 'anterior-cingulate',\
        'frontal-pole', 'cerebellar-hemisphere']
phenotypes = ['gm_volume', 'myelination', 'alff', 'reho_noGS', 'connmean_noGS']

## paths 
phen_dir = '/data1/rubinov_lab/brain_genomics/scripts_assoc_clean/inputs_{}/phen_regress'.format(dset)
expr_dir = '/data1/rubinov_lab/brain_genomics/scripts_assoc_clean/inputs_{}/expr_regress'.format(dset)
prm_file = '/data1/rubinov_lab/brain_genomics/scripts_assoc_clean/inputs_{}/permutations_100k.hdf5'.format(dset)  

out_main = '/data1/rubinov_lab/brain_genomics/scripts_assoc_clean/outputs_{}/assoc_1M/nulls'.format(dset)
if not os.path.exists(out_main): 
    os.mkdir(out_main) 
    for phen in phenotypes: 
        os.mkdir('{}/pvals_{}'.format(out_main, phen))

## permutations 
with h5py.File(prm_file, 'r') as f: 
    samp_order = np.array(f['samp_idx'], dtype=int) ## (subjs,)  
    subj_order = np.array(f['subj_idx'], dtype=int)[:set_perms] ## (perms, subjs) 

    if dset == 'HCP':
        twin_idxs = np.array(f['twin_idx'], dtype=int) ## (twin pairs, 2) 
        nontwin_idxs = np.array(f['nontwin_idx'], dtype=int) ## (nontwins,) 

nsubjs = samp_order.size 
    
## phenotype (original order) 
phens = {} ## k: (reg, phen), v: subject array 
for phen in phenotypes: 
    with h5py.File('{}/{}.hdf5'.format(phen_dir, phen), 'r') as f: 
        for reg in regions: 
            phens[(reg,phen)] = np.array(f[reg])[:, None] ## (subjs, 1)      

## expression (twin-sorted order if HCP)  
genes = {}; exprs = {} ## k: reg, v: gene array / (genes, subjs) expression matrix  
for reg in regions:
    with h5py.File('{}/{}.hdf5'.format(expr_dir, reg), 'r') as f: 
        genes[reg] = np.array([g.split('.')[0] for g in \
                     np.array(f['genes']).astype(str)])  
        exprs[reg] = np.array(f['pred_expr'])[:,samp_order] ## (genes, subjs)

## function: single gene associations
def find_signif_genes(pool_data):

    tstart = time()

    ## load parameters 
    perm = pool_data['idx'] ## set perm 
    phen = pool_data['phen'] 
    reg = pool_data['reg']

    ## check that assoc results don't already exist 
    rfile = '{}/pvals_{}/{}_{}.hdf5'.format(out_main, phen, reg, perm) 
    #try: 
    #    with h5py.File(rfile, 'r') as f: pass 
    #    return 
    #except: 
    #    pass 
    if os.path.exists(rfile): return  

    ## input data     
    phen_array = phens[(reg,phen)] ## (subjs, 1)   
    expr_matrx = exprs[reg] ## (genes, subjs) 
    gene_array = genes[reg]
    ngenes = gene_array.size

    ## set permuted phen array 
    pidx = subj_order[perm]  
    phen_ordered = phen_array[pidx].T ## (1, subjs) 

    ## phenotype permutations 
    if dset == 'HCP':
        X = np.random.random((assoc_perms, twin_idxs.shape[0])) ## generate random, unsorted numbers (perms, twin pairs)
        X = X.argsort(axis=1) ## sort to get unique index perms (perms, twin pairs)

        twin_perms = twin_idxs[X] ## yields permuted indices (perms, twin pairs, 2)
        twin_perms = twin_perms.reshape((assoc_perms, -1)) ## flatten twin pairs (perms, twins)

        nontwin_perms = np.repeat(nontwin_idxs[None,:], assoc_perms, axis=0) ## copy nontwin indices (perms, nontwins)
        nontwin_perms = rng.permuted(nontwin_perms, axis=1) ## use RNG to permute indices per perm (perms, nontwins)

        phen_perms = np.concatenate((twin_perms, nontwin_perms), axis=1) ## concatenate twins and nontwins (perms, subjs)

    elif dset == 'UKB':
        phen_perms = np.repeat(samp_order[None,:], assoc_perms, axis=0) ## copy subject indices (perms, subjs)
        phen_perms = rng.permuted(phen_perms, axis=1) ## use RNG to permute indices per perm (perms, subjs)

    ## compute observed gene rhos 
    c_expr = expr_matrx - expr_matrx.mean(axis=1)[:,None] 
    s_expr = (c_expr**2).sum(axis=1) 

    c_phen = phen_ordered - phen_ordered.mean() 
    s_phen = (c_phen**2).sum() 
    
    num = np.dot(c_phen, c_expr.T) 
    dom = np.sqrt(np.dot(s_phen, s_expr))
    rhos = num / dom ## (1, genes)  

    ## compute null gene rhos 
    c_phen = phen_perms - phen_perms.mean(axis=1)[:,None] 
    s_phen = (c_phen**2).sum(axis=1) 

    num = np.dot(c_phen, c_expr.T) 
    dom = np.sqrt(np.dot(s_phen[:,None], s_expr[None]))
    null_rhos = num / dom ## (perms, genes)  

    ## consider assocs where rho != NaN 
    ## (due to zero-variance expression)  
    valid_mask = np.isfinite(rhos).flatten() 
    valid_idxs = np.arange(rhos.size)[valid_mask] 

    ## compute p-values 
    pvals = (np.abs(null_rhos) >= np.abs(rhos)).mean(axis=0) 
    pvals[~valid_mask] = np.nan

    ## compute FDR-correct pvals 
    rejected, corrected = fdrcorrection(pvals[valid_mask]) 

    ## format all gene information 
    data = np.empty((ngenes, 3)) ## rho, pval, fdr
    data[:,0] = rhos
    data[:,1] = pvals 
    np.put(data[:,2], valid_idxs, corrected) 
    data[~valid_mask,2] = np.nan

    ## save 
    rfile = '{}/pvals_{}/{}_{}.hdf5'.format(out_main, phen, reg, perm) 
    with h5py.File(rfile, 'w') as f: 
        f['pearson'] = data 

    ## time 
    secs = time() - tstart
    hr = int(secs//3600)
    mn = int((secs%3600)//60)
    sc = int((secs%3600)%60)
    print('runtime: {:d} hr, {:d} mn, {:d} sc | {} {} {}'.format(hr, mn, sc, reg, phen, perm))

## run using Pool (reg phen perm)  
perm_range = np.arange(2000, 5000)
perm_range = np.arange(7540, 7549)
pdata = [{'phen':ph, 'reg':rg, 'idx':pm} \
        for pm in perm_range \
        for rg in regions \
        for ph in phenotypes] 
start = time()
pool = Pool(processes=28) 
pool.map(find_signif_genes, pdata)

secs = time() - start
hr = int(secs//3600)
mn = int((secs%3600)//60)
sc = int((secs%3600)%60)
print('total runtime: {:d} hr, {:d} mn, {:d} sc'.format(hr, mn, sc))
