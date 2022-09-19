'''
Script for running single gene associations, 
p-values are based on 1M permutations. 

- Nhung, July 2022 
'''

import sys 
import os 
import numpy as np  
import h5py 
from statsmodels.stats.multitest import fdrcorrection 
from multiprocessing import Pool 
from time import time 
import warnings 

## ignore warnings b/c these values are expected and set to NaNs
warnings.filterwarnings('ignore', 'invalid value encountered in true_divide')

dataset = sys.argv[1] ## e.g., HCP or UKB 
nperms = int(1e6)   
rng = np.random.default_rng() ## random number generator used for nontwin permutations 

phenotypes = ['gm_volume', 'alff', 'reho_noGS', 'connmean_noGS', 'myelination']
regions = ['hippocampus', 'amygdala', 'hypothalamus', 'substantia-nigra',\
        'caudate', 'putamen', 'nucleus-accumbens', 'anterior-cingulate',\
        'frontal-pole', 'cerebellar-hemisphere']

## paths 
phen_dir = '/data1/rubinov_lab/brain_genomics/analyses_{}/DATA_OUTPUT/phen_regress'.format(dataset)
expr_dir = '/data1/rubinov_lab/brain_genomics/analyses_{}/DATA_OUTPUT/r0.3_p0.01/expr_regress'.format(dataset)
prm_file = '/data1/rubinov_lab/brain_genomics/analyses_{}/DATA_OUTPUT/null_permutations_latest_100k.hdf5'.format(dataset)  

out_main = '/data1/rubinov_lab/brain_genomics/analyses_{}/assoc_1M'.format(dataset)
if not os.path.exists('{}/pvals_{}'.format(out_main, phenotypes[0])):
    for phen in phenotypes: 
        os.mkdir('{}/pvals_{}'.format(out_main, phen))

## permutations 
with h5py.File(prm_file, 'r') as f: 
    samp_order = np.array(f['samp_idx'], dtype=int) ## (subjs,)

    if dataset == 'HCP':
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

    ## load parameters 
    reg = pool_data['reg']
    phen = pool_data['phen'] 
    nperms = pool_data['num_perms']

    ## check that assoc results don't already exist  
    out_path = '{}/pvals_{}/{}.hdf5'.format(out_main, phen, reg) 
    if os.path.exists(out_path): 
        print('done:', reg, phen) 
        return 

    ## input data 
    phen_array = phens[(reg,phen)] ## (subjs, 1)   
    expr_matrx = exprs[reg] ## (genes, subjs) 
    gene_array = genes[reg]
    ngenes = gene_array.size

    ## phenotype permutations 
    if dataset == 'HCP':
        X = np.random.random((nperms, twin_idxs.shape[0])) ## generate random, unsorted numbers (perms, twin pairs)  
        X = X.argsort(axis=1) ## sort to get unique index perms (perms, twin pairs)

        twin_perms = twin_idxs[X] ## yields permuted indices (perms, twin pairs, 2)
        twin_perms = twin_perms.reshape((nperms, -1)) ## flatten twin pairs (perms, twins)

        nontwin_perms = np.repeat(nontwin_idxs[None,:], nperms, axis=0) ## copy nontwin indices (perms, nontwins) 
        nontwin_perms = rng.permuted(nontwin_perms, axis=1) ## use RNG to permute indices per perm (perms, nontwins) 

        phen_perms = np.concatenate((twin_perms, nontwin_perms), axis=1) ## concatenate twins and nontwins (perms, subjs)

    elif dataset == 'UKB': 
        phen_perms = np.repeat(samp_order[None,:], nperms, axis=0) ## copy subject indices (perms, subjs) 
        phen_perms = rng.permuted(phen_perms, axis=1) ## use RNG to permute indices per perm (perms, subjs) 
        
    ## compute observed gene rhos 
    c_expr = expr_matrx - expr_matrx.mean(axis=1)[:,None] 
    s_expr = (c_expr**2).sum(axis=1) 

    phen_ordered = phen_array[samp_order].T ## (1, subjs) 
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

    ## save results 
    out_path = '{}/pvals_{}/{}_1M.hdf5'.format(out_main, phen, reg) 
    with h5py.File(out_path, 'w') as f: 
        f['genes'] = gene_array.astype(bytes)
        f['pearson'] = data ## rho, pval, fdr 

    ## print summary 
    p_sig = np.sum(data[:,1][valid_mask] <= 0.005)
    print('{:>3d} sig ({} {})'.format(p_sig, reg, phen))
    #f_sig = np.sum(data[:,2][valid_mask] <= 0.25)
    #print('> {:>3d} p, {:>3d} f - {} {} ({} perms)'.format(p_sig, f_sig, reg, phen, nperms))

pdata = [{'num_perms':nperms, 'reg':reg, 'phen':phen} for reg in regions for phen in phenotypes]
tstart = time() 
pool = Pool(processes=30)
pool.map(find_signif_genes, pdata)

secs = time() - tstart
hr = int(secs//3600)
mn = int((secs%3600)//60)
sc = int((secs%3600)%60)
print('runtime: {:d} hr, {:d} mn, {:d} sc'.format(hr, mn, sc))
