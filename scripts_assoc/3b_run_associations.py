'''
https://stackoverflow.com/questions/30143417/computing-the-correlation-coefficient-between-two-multi-dimensional-arrays
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

dataset = sys.argv[1] ## e.g., HCP or UKB 

phenotypes = ['gm_volume', 'alff', 'reho_noGS', 'connmean_noGS', 'myelination']
regions = ['hippocampus', 'amygdala', 'hypothalamus', 'substantia-nigra',\
        'caudate', 'putamen', 'nucleus-accumbens', 'anterior-cingulate',\
        'frontal-pole', 'cerebellar-hemisphere']
reg_phens = [(reg, phen) for reg in regions for phen in phenotypes] 

## paths 
phen_dir = '/data1/rubinov_lab/brain_genomics/analyses_{}/DATA_OUTPUT/phen_regress'.format(dataset)
expr_dir = '/data1/rubinov_lab/brain_genomics/analyses_{}/DATA_OUTPUT/r0.3_p0.01/expr_regress'.format(dataset)
prm_file = '/data1/rubinov_lab/brain_genomics/analyses_{}/DATA_OUTPUT/null_permutations_latest_100k.hdf5'.format(dataset)  

out_main = '/data1/rubinov_lab/brain_genomics/analyses_{}/assoc_750k'.format(dataset)
if not os.path.exists('{}/pvals_{}'.format(out_main, phenotypes[0])):
    for phen in phenotypes: 
        os.mkdir('{}/pvals_{}'.format(out_main, phen))

## permutations 
with h5py.File(prm_file, 'r') as f: 
    samp_order = np.array(f['samp_idx'], dtype=int) ## (subjs,)  
    twin_idxs = np.array(f['twin_idx'], dtype=int) ## (twin pairs, 2) 
    nontwin_idxs = np.array(f['nontwin_idx'], dtype=int) ## (nontwins,) 

## phenotype (original order) 
phens = {} ## k: (reg, phen), v: subject array 
for phen in phenotypes: 
    with h5py.File('{}/{}.hdf5'.format(phen_dir, phen), 'r') as f: 
        for reg in regions: 
            phens[(reg,phen)] = np.array(f[reg])[:, None] ## (subjs, 1)      

## expression (twin-sorted order)  
genes = {}; exprs = {} ## k: reg, v: gene array / (genes, subjs) expression matrix  
for reg in regions:
    with h5py.File('{}/{}.hdf5'.format(expr_dir, reg), 'r') as f: 
        genes[reg] = np.array([g.split('.')[0] for g in \
                     np.array(f['genes']).astype(str)]) 
        exprs[reg] = np.array(f['pred_expr'])[:,samp_order] ## (genes, subjs)

nsubjs = phens[(regions[0], phenotypes[0])].size 

##
def get_permutation(): 
    twin_perm = np.random.permutation(twin_idxs) 
    nontwin_perm = np.random.permutation(nontwin_idxs) 
    perm_idx = np.concatenate((twin_perm.flatten(), nontwin_perm)) 
    return perm_idx 

## function: save gene-wise results 
## (gene name, rho, pval, fdr)  
def save_results(gene_array, data, filename):
    with h5py.File(filename, 'w') as f: 
        f['genes'] = gene_array.astype(bytes)
        f['pearson'] = data ## rho, pval, fdr 

    #with open(filename, 'w') as f:
    #    header = '\t'.join(['GENE', 'RHO', 'PVAL', 'FDR', '\n'])
    #    f.write(header)
    #    for g,d in zip(gene_array,data):
    #        line = '{}\t{:.8f}\t{:.8f}\t{:.8f}\n'.format(g.split('.')[0], d[0], d[1], d[2])
    #        f.write(line)

## function: single gene associations
def find_signif_genes(pool_data):
    reg = pool_data['reg']
    phen = pool_data['phen'] 
    nperms = pool_data['num_perms']

    ## input data 
    phen_array = phens[(reg,phen)] ## (subjs, 1)   
    expr_matrx = exprs[reg] ## (genes, subjs) 
    gene_array = genes[reg]
    ngenes = gene_array.size

    ## phenotype permutations 
    phen_perms = np.zeros((nperms, nsubjs)) ## (perms, subjs) 
    for p in range(nperms): 
        phen_perms[p] = phen_array[get_permutation(),0]  

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
    out_path = '{}/pvals_{}/{}.hdf5'.format(out_main, phen, reg) 
    save_results(gene_array, data, out_path)

    ## print summary 
    p_sig = np.sum(data[:,1][valid_mask] <= 0.005)
    print('{:>3d} sig ({} {})'.format(p_sig, reg, phen))
    #f_sig = np.sum(data[:,2][valid_mask] <= 0.25)
    #print('> {:>3d} p, {:>3d} f - {} {} ({} perms)'.format(p_sig, f_sig, reg, phen, nperms))

nperms = 750 * 1000  
pdata = [{'num_perms':nperms, 'reg':reg, 'phen':phen} for reg in regions for phen in phenotypes]
tstart = time() 
pool = Pool(processes=30)
pool.map(find_signif_genes, pdata)
tend = time() 
tpass = (tend - tstart)  
print('runtime: {:.2f} secs'.format(tpass))
