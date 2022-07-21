'''


null version of single gene assoc 

multithreading to update (perms * 3) array 

each reg phen has 10k permutations, whose sig is 
computed using another 10k permutations 

each phen has a dir in nulls/ and each reg has a file 

array is in regional loop (big reg phen loop) 
'''

import sys 
import os 
import numpy as np  
import h5py 
from scipy.stats import pearsonr, spearmanr 
from statsmodels.stats.multitest import fdrcorrection 
from multiprocessing import Pool  
from time import time 
import warnings 

warnings.filterwarnings('ignore', 'invalid value encountered in true_divide')

dset = sys.argv[1] ## HCP or UKB   
assoc_perms = 750 * 1000
set_perms = 10 * 1000

regions = ['hippocampus', 'amygdala', 'hypothalamus', 'substantia-nigra',\
        'caudate', 'putamen', 'nucleus-accumbens', 'anterior-cingulate',\
        'frontal-pole', 'cerebellar-hemisphere']
phenotypes = ['gm_volume', 'myelination', 'alff', 'reho_noGS', 'connmean_noGS']

## paths 
phen_dir = '/data1/rubinov_lab/brain_genomics/analyses_{}/DATA_OUTPUT/phen_regress'.format(dset)
expr_dir = '/data1/rubinov_lab/brain_genomics/analyses_{}/DATA_OUTPUT/r0.3_p0.01/expr_regress'.format(dset)
prm_file = '/data1/rubinov_lab/brain_genomics/analyses_{}/DATA_OUTPUT/null_permutations_latest_100k.hdf5'.format(dset)  

out_main = '/data1/rubinov_lab/brain_genomics/analyses_{}/assoc_750k/nulls'.format(dset)
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
        genes[reg] = np.array(f['genes']).astype(str) 
        exprs[reg] = np.array(f['pred_expr'])[:,samp_order] ## (genes, subjs)

## function: get permutation wrt original subj order
def get_permutation(): 
    if dset == 'HCP':
        twin_perm = np.random.permutation(twin_idxs) 
        nontwin_perm = np.random.permutation(nontwin_idxs) 
        perm_idx = np.concatenate((twin_perm.flatten(), nontwin_perm)) 
    elif dset == 'UKB': 
        perm_idx = np.random.permutation(samp_order) 
    return perm_idx 

## function: single gene associations
def find_signif_genes(pool_data):
    perm = pool_data['idx'] ## set perm 
    phen = pool_data['phen'] 
    reg = pool_data['reg']

    rfile = '{}/pvals_{}/{}_{}.hdf5'.format(out_main, phen, reg, perm) 
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
    phen_perms = np.zeros((assoc_perms, nsubjs)) ## (perms, subjs) 
    for p in range(assoc_perms): 
        phen_perms[p] = phen_array[get_permutation(),0]  

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
        f['res'] = data 
    print(phen, reg, perm) 

    #return data 

#for phen in ['regional_homogeneity']: 
#    for reg in regions: 
#        ngenes = genes[reg].size  
#        assoc_results = np.empty((nperms, ngenes, 3))
#
#        ## run associations  
#        for i in range(nperms): 
#            pdata = {'idx':i, 'reg':reg, 'phen':phen}
#            assoc_results[i] = find_signif_genes(pdata)
#            if (i%100)==0: print(i)
#
#        ## save assoc results 
#        out_file = '{}/pvals_{}/{}.hdf5'.format(out_main, phen, reg) 
#        with h5py.File(out_file, 'w') as f: 
#            f['assoc_results'] = assoc_results 
#        print('DONE: {} {}'.format(phen, reg))

## pool (reg phen perm) 
pdata = [{'phen':ph, 'reg':rg, 'idx':pm} \
    for ph in phenotypes for rg in regions for pm in np.arange(set_perms)] 
start = time()
pool = Pool(processes=15) 
pool.map(find_signif_genes, pdata)
duration = time() - start 

mn, sc = divmod(duration, 60) 
hr, mn = divmod(mn, 60) 
print('{:d} hr, {:d} mn, {:d} sc'.format(hr, mn, sc))
