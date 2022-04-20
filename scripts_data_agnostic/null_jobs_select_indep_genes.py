'''
Permutation testing:
  Compute single-gene and regional-phenotype correlations.

Usage: First argument is phenotype (all regional variations
  of this phenotype will be computed over). 

Notes:
- The gene models that are being considered have already been
  filtered based on completeness or quality thresholds.
- Permutation testing is used to compute p-values, and twin
  structures are maintained during the shuffling.
- Input expression and phenotype data should already be in their
  residual form (i.e. age, gender, PC1, PC2 confounders have been
  regressed out).

- Nhung, updated April 2022
'''

import sys 
import os 
import numpy as np  
import h5py 
from scipy.stats import pearsonr, spearmanr 
from statsmodels.stats.multitest import fdrcorrection 

## script arguments
dataset = sys.argv[1] ## e.g., HCP, UKB 

phen_reg_perm = int(sys.argv[2])
phen_idx = int(phen_reg_perm / 10000) % 4
reg_idx = int((phen_reg_perm % 10000) / 1000)
perm_idx = int((phen_reg_perm % 10000) % 1000)

## declare regional phenotype permutation  
atlas = 'hoacer_sn_hth'
regions = ['caudate', 'nucleus-accumbens', 'hippocampus', 'amygdala', 'anterior-cingulate', \
    'putamen', 'hypothalamus', 'substantia-nigra', 'frontal-pole', 'cerebellar-hemisphere']
phenotypes = ['alff', 'regional_homogeneity', 'gm_volume', 'connectivity_mean']

PHN = phenotypes[phen_idx] 
REG = regions[reg_idx] 
PRM = perm_idx 

## paths 
phn_file = '/data/rubinov_lab/brain_genomics_project/platypus7/phen_regress/{}.hdf5'.format(PHN)
expr_dir = '/data/rubinov_lab/brain_genomics_project/platypus7/expr_regress'
perm_file = '/data/rubinov_lab/brain_genomics_project/platypus7/null_permutations.hdf5' 
out_main = '/data/rubinov_lab/brain_genomics_project/platypus7/null_pvals_{}'.format(PHN)

#if phen_reg_perm == 0:
#    for phen in phenotypes:
#        path = '/data/rubinov_lab/brain_genomics_project/platypus6/null_pvals_{}'.format(phen)
#        if not os.path.exists(path): os.mkdir(path) 
#        for reg in regions: 
#            os.mkdir(path + '/' + reg)

## function: compute p-value for specific gene association 
def compute_pvalue(gene_expr, full_phen, rho, perms): 
    N = 1000
    null = np.array([pearsonr(gene_expr[perms[n]], full_phen)[0] for n in range(N)]) 
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

    ## permutations  
    with h5py.File(perm_file, 'r') as f: 
        obsv_subj_order = np.array(f['subj_idx'], dtype=int) ## (890,)
        obsv_samp_order = np.array(f['samp_idx'][PRM], dtype=int) ## (1000, 890) -> (890,) 
        null_samp_perms = np.array(f['null_idx'][PRM], dtype=int) ## (1000, 1000, 890) -> (1000, 890)  

    ## expression (original order)  
    genes = {}; exprs = {} 
    with h5py.File('{}/{}.hdf5'.format(expr_dir, REG), 'r') as f:
        genes[REG] = np.array([g.decode("utf-8") for g in np.array(f['genes'])])
        exprs[REG] = np.array(f['pred_expr']) ## (genes * samps)

    ## phenotype (original order)  
    with h5py.File(phn_file, 'r') as f:
        phens = {REG: np.array(f[REG])}

    ## input data 
    phen_array = phens[REG][obsv_subj_order]## in twin-sorted order  
    gene_array = genes[REG] ## as is 

    expr_matrix = exprs[REG] ## in original order (needed for this null's perms)  
    expr_sorted = expr_matrix[:,obsv_samp_order] ## in twin-sorted order (for this null)  

    ngenes = gene_array.size

    ## gene loop starts here 
    data = np.zeros((ngenes, 3)) ## per gene: rho, pval, fdr
    not_nan = [] 
    print('  {:>.2f}% [{}/{}] - {}'.format(0, 0, ngenes, REG))
    for g,gene in enumerate(gene_array): 
        
        ## return nan for genes with no expression variance  
        if np.var(expr_sorted[g]) == 0: 
            data[g,:] = np.nan
            continue  

        ## otherwise, compute correlation 
        not_nan.append(g) 
        rho, _ = pearsonr(expr_sorted[g], phen_array)
        pval = compute_pvalue(expr_matrix[g], phen_array, rho, null_samp_perms)  
        data[g,0] = rho
        data[g,1] = pval 

        ## status 
        if (g+1)%50 == 0: 
            perc = ((g+1)/ngenes) * 100 
            print('  {:>.2f}% [{}/{}] - {} ({})'.format(perc, g+1, ngenes, REG, PRM))

    ## compute FDR-correct pvals 
    pvals = data[not_nan][:,1]
    rejected, corrected = fdrcorrection(pvals) 
    np.put(data[:,2], not_nan, corrected)

    ## save results 
    out_path = '{}/{}/{}.txt'.format(out_main, REG, PRM) 
    save_results(gene_array, data, out_path)

    ## print summary 
    #p_sig = np.sum(data[:,1] <= 0.05)
    #f_sig = np.sum(data[:,2] <= 0.05)
    #print('> {:>3d} p, {:>3d} f - {} ({})\n'.format(p_sig, f_sig, reg, obsv_perm))

find_signif_genes()
