'''
Permutation testing:
  Compute single-gene and regional-phenotype correlations.

Usage: First argument is phenotype (all regional variations
  of this phenotype will be computed over). The current script
  uses pooling (rather than job arrays).

Notes:
- The gene models that are being considered have already been
  filtered based on completeness or quality thresholds.
- Permutation testing is used to compute p-values, and twin
  structures are maintained during the shuffling.
- Input expression and phenotype data should already be in their
  residual form (i.e. age, gender, PC1, PC2 confounders have been
  regressed out).

- Nhung, updated Jan 2022
'''

import sys 
import os 
import numpy as np  
import h5py 
from scipy.stats import pearsonr, spearmanr 
from statsmodels.stats.multitest import fdrcorrection 
from multiprocessing import Pool 

atlas = 'hoacer_sn_hth'
regions = ['caudate', 'nucleus-accumbens', 'hippocampus', 'amygdala', 'anterior-cingulate', \
    'putamen', 'hypothalamus', 'substantia-nigra', 'frontal-pole', 'cerebellar-hemisphere']
phenotypes = ['alff', 'regional_homogeneity', 'gm_volume']

dataset = sys.argv[1] ## e.g., HCP, UKB 

phen_reg_perm = int(sys.argv[2])
aphen = int(phen_reg_perm / 1000) % 3
areg = int((phen_reg_perm % 1000) / 100)
aperm = int((phen_reg_perm % 1000) % 100)

areg_name = regions[areg]
phenotype = phenotypes[aphen] 

## paths 
phn_file = '/data/rubinov_lab/brain_genomics_project/platypus6/phen_regress/{}.hdf5'.format(phenotype)
expr_dir = '/data/rubinov_lab/brain_genomics_project/platypus6/expr_regress'

obsv_perm_file = '/data/rubinov_lab/brain_genomics_project/platypus6/nulls_of_multi.hdf5' 
null_perm_file = '/data/rubinov_lab/brain_genomics_project/platypus6/nulls_of_nulls.hdf5' 

out_main = '/data/rubinov_lab/brain_genomics_project/platypus6/null_pvals_{}'.format(phenotype)

#if phen_reg_perm == 0:
#    for phen in phenotypes:
#        path = '/data/rubinov_lab/brain_genomics_project/platypus6/null_pvals_{}'.format(phen)
#        if not os.path.exists(path): os.mkdir(path) 
#        for reg in regions: 
#            os.mkdir(path + '/' + reg)

## permutations of the observed data 
with h5py.File(obsv_perm_file, 'r') as f: 
    obsv_subj_order = np.array(f['subj_idx'])
    obsv_samp_perms = {i:np.array(f['samp_idx_{}'.format(i)]) for i in range(100)}

## permutations of the null data 
with h5py.File(null_perm_file, 'r') as f: 
    null_subj_order = np.array(f['subj_idx'])
    null_samp_perms = {i:np.array(f['samp_idx_{}'.format(i)]) for i in range(1000)}
assert(np.array_equal(obsv_subj_order, null_subj_order))

## phenotype (sorted in twin order)  
with h5py.File(phn_file, 'r') as f:
    phens = {areg_name: np.array(f[areg_name])[obsv_subj_order]}

## expression  
genes = {}; exprs = {} 
with h5py.File('{}/{}.hdf5'.format(expr_dir, areg_name), 'r') as f:
    genes[areg_name] = np.array([g.decode("utf-8") for g in np.array(f['genes'])])
    exprs[areg_name] = np.array(f['pred_expr']) ## (genes * samps)

## function: compute p-value for specific gene association 
def compute_pvalue(expr, full_phen, rho): 
    N = 1000
    null = np.array([pearsonr(expr[null_samp_perms[n]], full_phen)[0] for n in range(N)])
    pval = np.mean(np.abs(null) >= np.abs(rho))
    return pval 

## function: save gene-wise results 
## (gene name, rho, pval, fdr)  
def save_results(gene_array, data, filename):
    with open(filename, 'w') as f:
        header = '\t'.join(['GENE', 'RHO', 'PVAL', 'FDR', '\n'])
        f.write(header) 
        for g,d in zip(gene_array,data): 
            line = '{}\t{:.5f}\t{:.5f}\t{:.5f}\n'.format(g, d[0], d[1], d[2])
            f.write(line) 

## function: single gene associations
def find_signif_genes():

    reg = areg_name
    obsv_perm = aperm
    obsv_perm_order = obsv_samp_perms[obsv_perm]

    ## input data 
    phen_array = phens[reg] ## in twin-sorted order  
    gene_array = genes[reg]

    expr_matrix = exprs[reg] ## in original order  
    expr_sorted = expr_matrix[:,obsv_perm_order] ## in twin-sorted order (for this obsv perm)  

    ngenes = gene_array.size

    ## gene loop starts here 
    data = np.zeros((gene_array.shape[0], 3)) ## rho, pval, fdr
    not_nan = [] 
    print('  {:>.2f}% [{}/{}] - {}'.format(0, 0, ngenes, reg))
    for g,gene in enumerate(gene_array): 
        
        ## return nan for genes with no expression variance  
        if np.var(expr_sorted[g]) == 0: 
            data[g,:] = np.nan
            continue  

        ## otherwise, compute correlation 
        not_nan.append(g) 
        rho, _ = pearsonr(phen_array, expr_sorted[g])
        pval = compute_pvalue(expr_matrix[g], phen_array, rho)  
        data[g,0] = rho
        data[g,1] = pval 

        ## status 
        if (g+1)%50 == 0: 
            perc = ((g+1)/ngenes) * 100 
            print('  {:>.2f}% [{}/{}] - {} ({})'.format(perc, g+1, ngenes, reg, obsv_perm))

    ## compute FDR-correct pvals 
    pvals = data[not_nan][:,1]
    rejected, corrected = fdrcorrection(pvals) 
    np.put(data[:,2], not_nan, corrected)

    ## save results 
    out_path = '{}/{}/{}.txt'.format(out_main, reg, obsv_perm) 
    save_results(gene_array, data, out_path)

    ## print summary 
    p_sig = np.sum(data[:,1] <= 0.05)
    f_sig = np.sum(data[:,2] <= 0.05)
    print('> {:>3d} p, {:>3d} f - {} ({})\n'.format(p_sig, f_sig, reg, obsv_perm))

find_signif_genes()
