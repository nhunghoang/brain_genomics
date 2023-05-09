'''
Script for running single gene associations, 
p-values are based on 1M permutations. 

- Nhung, updated April 2023
'''

import sys 
import os 
import numpy as np  
import h5py 
from statsmodels.stats.multitest import fdrcorrection 
from multiprocessing import Pool 
from time import time 
import warnings 

from scipy.stats import pearsonr

## ignore warnings b/c these values are expected and set to NaNs
warnings.filterwarnings('ignore', 'invalid value encountered in true_divide')

## random number generator used for permutations
rng = np.random.default_rng() 

## args
group = sys.argv[1] ## e.g., HCP or UKB 
model = sys.argv[2] ## e.g., PDX, JTI, UTM, FUS
nperms = int(1e6) 

## paths 
top_path = '/data1/rubinov_lab/brain_genomics/scripts_twas/'
phen_dir = top_path + 'inputs_{}/phen_regress_{}'.format(group, model)
expr_dir = top_path + 'inputs_{}/expr_regress_{}'.format(group, model)
outs_dir = top_path + 'outputs_{}/assoc_{}'.format(group, model)

if group == 'HCP': 
    perm_file = top_path + 'inputs_HCP/permutations_100k_euro.hdf5'

## regions and phenotypes 
regions = ['hippocampus', 'amygdala', 'caudate', 'putamen', \
           'nucleus-accumbens', 'cerebellar-hemisphere']
phenotypes = ['left_vol', 'right_vol']
#phenotypes = ['gm_volume', 'alff', 'reho_noGS', 'connmean_noGS', 'myelination']
#regions = ['hippocampus', 'amygdala', 'hypothalamus', 'substantia-nigra',\
#        'caudate', 'putamen', 'nucleus-accumbens', 'anterior-cingulate',\
#        'frontal-pole', 'cerebellar-hemisphere']

## create output folders if non-existent 
if not os.path.exists('{}/pvals_{}'.format(outs_dir, phenotypes[0])):
    for phen in phenotypes: 
        os.mkdir('{}/pvals_{}'.format(outs_dir, phen))

## BIG TMP 
#with h5py.File(perm_file, 'r') as f: 
#    euro_idx = f['euro_idx'][()].astype(bool) 

## phenotype (original order) 
phens = {} ## k: (reg, phen), v: subject array 
for phen in phenotypes: 
    with h5py.File('{}/{}.hdf5'.format(phen_dir, phen), 'r') as f: 
        for reg in regions: 
            #phens[(reg,phen)] = f[reg][()][:, None] ## (subjs, 1)      
            phens[(reg,phen)] = f[reg][()] ## (subjs,)      
            
            ## TMP 
            #phens[(reg,phen)] = f[reg][()][:, None][euro_idx] ## (subjs, 1)      

nsubjs = phens[(regions[0], phenotypes[0])].size 

## get sample indices
if group == 'HCP':
    with h5py.File(perm_file, 'r') as f: 
        samp_order = f['samp_idx'][()].astype(int) ## (subjs,)
        twin_idxs = f['twin_idx'][()].astype(int) ## (twin pairs, 2)
        nontwin_idxs = f['nontwin_idx'][()].astype(int) ## (nontwins,)

elif group == 'UKB': 
    samp_order = np.arange(nsubjs, dtype='int32') ## (subjs,)
    
## expression (twin-sorted order if HCP, else original order)  
genes = {}; exprs = {} ## k: reg, v: gene array / (genes, subjs) expression matrix  
for reg in regions:
    with h5py.File('{}/{}.hdf5'.format(expr_dir, reg), 'r') as f: 
        genes[reg] = f['genes'][()].astype(str)
        #exprs[reg] = f['pred_expr'][()][:,euro_idx][:,samp_order] ## (genes, subjs)
        exprs[reg] = f['pred_expr'][()][:,samp_order] ## (genes, subjs)

## TMP get valid UKB genes 
gpath = '/data1/rubinov_lab/brain_genomics/models_JTI/genes_r30_p01/'
for reg in regions: 
    jti_genes = np.loadtxt(gpath + reg + '.csv', dtype=str)
    mask = np.isin(genes[reg], jti_genes)
    genes[reg] = genes[reg][mask]
    exprs[reg] = exprs[reg][mask]
    print('{}: {} --> {}'.format(reg, mask.size, mask.sum())) 



## set phenotype permutation indices 
#if group == 'HCP':
#    X = np.random.random((nperms, twin_idxs.shape[0])) ## generate random, unsorted numbers (perms, twin pairs)  
#    X = X.argsort(axis=1) ## sort to get unique index perms (perms, twin pairs)
#
#    twin_perm_idxs = twin_idxs[X] ## yields permuted indices (perms, twin pairs, 2)
#    twin_perm_idxs = twin_perm_idxs.reshape((nperms, -1)) ## flatten twin pairs (perms, twins)
#
#    nontwin_perm_idxs = np.repeat(nontwin_idxs[None,:], nperms, axis=0) ## copy nontwin indices (perms, nontwins) 
#    nontwin_perm_idxs = rng.permuted(nontwin_perm_idxs, axis=1) ## use RNG to permute indices per perm (perms, nontwins) 
#
#    perm_idxs = np.concatenate((twin_perm_idxs, nontwin_perm_idxs), axis=1) ## concatenate twins and nontwins (perms, subjs)
#
#elif group == 'UKB': 
#    perm_idxs = np.repeat(samp_order[None,:], nperms, axis=0) ## copy subject indices (perms, subjs) 
#    perm_idxs = rng.permuted(perm_idxs, axis=1) ## use RNG to permute indices per perm (perms, subjs) 
        
## function: single gene associations
def find_signif_genes(pool_data):

    ## load parameters 
    reg = pool_data['reg']
    phen = pool_data['phen'] 

    ## check that assoc results don't already exist  
    #out_path = '{}/pvals_{}/{}.hdf5'.format(out_dir, phen, reg) 
    #if os.path.exists(out_path): 
    #    print('done:', reg, phen) 
    #    return 

    ## input data 
    phen_array = phens[(reg,phen)] ## (subjs, 1)   
    expr_matrx = exprs[reg] ## (genes, subjs) 
    gene_array = genes[reg]
    ngenes = gene_array.size
        
    ## compute observed gene rhos 
    #c_expr = expr_matrx - expr_matrx.mean(axis=1)[:,None] 
    #s_expr = (c_expr**2).sum(axis=1) 

    #phen_ordered = phen_array[samp_order].T ## (1, subjs) 
    #c_phen = phen_ordered - phen_ordered.mean() 
    #s_phen = (c_phen**2).sum() 
    #
    #num = np.dot(c_phen, c_expr.T) 
    #dom = np.sqrt(np.dot(s_phen, s_expr))
    #rhos = num / dom ## (1, genes)  

    ## compute null gene rhos 
    #phen_perms = phen_array[:,0][perm_idxs] 
    #c_phen = phen_perms - phen_perms.mean(axis=1)[:,None] 
    #s_phen = (c_phen**2).sum(axis=1) 

    #num = np.dot(c_phen, c_expr.T) 
    #dom = np.sqrt(np.dot(s_phen[:,None], s_expr[None]))
    #null_rhos = num / dom ## (perms, genes)  

    ## TMP P-VALS
    rhos = np.zeros(ngenes); pvals = np.zeros(ngenes) 
    for g in range(ngenes): 
        r, p = pearsonr(phen_array, expr_matrx[g])
        rhos[g] = r
        pvals[g] = p 

    ## consider assocs where rho != NaN 
    ## (due to zero-variance expression)  
    valid_mask = np.isfinite(rhos).flatten() 
    valid_idxs = np.arange(rhos.size)[valid_mask] 

    ## compute p-values 
    #pvals = (np.abs(null_rhos) >= np.abs(rhos)).mean(axis=0) 
    #pvals[~valid_mask] = np.nan

    ## compute FDR-correct pvals 
    rejected, corrected = fdrcorrection(pvals[valid_mask]) 

    ## format all gene information 
    data = np.empty((ngenes, 3)) ## rho, pval, fdr
    data[:,0] = rhos
    data[:,1] = pvals 
    np.put(data[:,2], valid_idxs, corrected) 
    data[~valid_mask,2] = np.nan

    ## save results 
    out_path = '{}/pvals_{}/{}.hdf5'.format(outs_dir, phen, reg) 
    with h5py.File(out_path, 'w') as f: 
        f['genes'] = gene_array.astype(bytes)
        f['pearson'] = data ## rho, pval, fdr 

    ## print summary 
    p_sig = np.sum(data[:,1][valid_mask] <= 0.005)
    print('{:>3d} sig ({} {})'.format(p_sig, reg, phen))
    #f_sig = np.sum(data[:,2][valid_mask] <= 0.25)
    #print('> {:>3d} p, {:>3d} f - {} {} ({} perms)'.format(p_sig, f_sig, reg, phen, nperms))

#pdata = [{'reg':reg, 'phen':phen} for reg in regions for phen in phenotypes]
#find_signif_genes(pdata[0])
#import sys; sys.exit()

pdata = [{'reg':reg, 'phen':phen} for reg in regions for phen in phenotypes]
tstart = time() 
pool = Pool(processes=10)
pool.map(find_signif_genes, pdata)

secs = time() - tstart
hr = int(secs//3600)
mn = int((secs%3600)//60)
sc = int((secs%3600)%60)
print('runtime: {:d} hr, {:d} mn, {:d} sc'.format(hr, mn, sc))
