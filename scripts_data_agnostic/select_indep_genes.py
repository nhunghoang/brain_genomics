'''
Compute single-gene and regional-phenotype correlations. 

Usage: First argument is phenotype (all regional variations
  of this phenotype will be computed over). The current script
  uses pooling (rather than job arrays).  

Notes: 
- The gene models that are being considered have already been 
  filtered based on completeness or quality thresholds. 
- Permutation testing is used to compute p-values, and twin 
  structures are maintained during the shuffling (for HCP).  
- Input expression and phenotype data should already be in their 
  residual form (i.e. age, gender, PC1, PC2 confounders have been 
  regressed out). 

- Nhung, updated Jan 2022

No change for UKB

- Tim, Mar 15 2022
'''

import sys 
import os 
import numpy as np  
import h5py 
from scipy.stats import pearsonr, spearmanr 
from statsmodels.stats.multitest import fdrcorrection 
from multiprocessing import Pool 

atlas = 'hoacer_sn_hth'
phenotype = sys.argv[1] 
dataset = sys.argv[2] ## e.g., HCP or UKB 

## needed for HCP diffusion subjects only 
famd_flag = False 
nan_idx = None 
idx_path = '/data1/rubinov_lab/brain_genomics/data_HCP/hoacer_sn_hth/phenotypes/famd_isnan.txt'
if (dataset == 'HCP') and (phenotype in ['fa', 'md']):
    nan_idx = np.loadtxt(idx_path, dtype=int)
    nan_idx = nan_idx.astype(bool)
    famd_flag = True 

## paths 
phn_file = '/data1/rubinov_lab/brain_genomics/analyses_{}/DATA_OUTPUT/phen_regress/{}.hdf5'.format(dataset, phenotype)
expr_dir = '/data1/rubinov_lab/brain_genomics/analyses_{}/DATA_OUTPUT/expr_regress'.format(dataset)
prm_file = '/data1/rubinov_lab/brain_genomics/analyses_{}/DATA_OUTPUT/null_permutations.hdf5'.format(dataset)  
if famd_flag:
    prm_file = '/data1/rubinov_lab/brain_genomics/analyses_{}/DATA_OUTPUT/null_permutations_10k_famd.hdf5'.format(dataset)  

out_main = '/data1/rubinov_lab/brain_genomics/analyses_{}/assoc2/pvals_{}'.format(dataset, phenotype)
if not os.path.exists(out_main): os.mkdir(out_main) 

## gather permutations 
## permutations (apply to subjects, but not samples)  
## because sample permutation indices are wrt the original 890 indices
with h5py.File(prm_file, 'r') as f: 
    subj_order = np.array(f['subj_idx'], dtype=int) ## (890,)  
    samp_perms = np.array(f['samp_idx'], dtype=int) ## (1000, 890) 

## phenotype  
with h5py.File(phn_file, 'r') as f: 
    phens = {reg: np.array(f[reg])[subj_order] for reg in f.keys()}
regions = np.array(list(phens.keys()))

## expression  
genes = {}; exprs = {} 
for expr_file in os.listdir(expr_dir): 
    if expr_file[-5:] != '.hdf5': continue 
    with h5py.File('{}/{}'.format(expr_dir, expr_file), 'r') as f: 
        reg = expr_file.split('.')[0]
        genes[reg] = np.array([g.decode("utf-8") for g in np.array(f['genes'])])
        if famd_flag: exprs[reg] = np.array(f['pred_expr'])[:,~nan_idx] ## (genes * samps)
        else: exprs[reg] = np.array(f['pred_expr']) ## (genes * samps)

## function: compute p-value for specific gene association 
def compute_pvalue(expr, full_phen, rho): 
    N = 1000 ## Feb 2022 update: use 1k nulls, instead of 10k
    null = np.array([pearsonr(expr[samp_perms[n]], full_phen)[0] for n in range(N)])
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
def find_signif_genes(reg):

    ## input data 
    phen_array = phens[reg] ## in twin-sorted order  
    gene_array = genes[reg]

    expr_matrix = exprs[reg] ## in original order  
    expr_sorted = expr_matrix[:,subj_order] ## in twin-sorted order 

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
            print('  {:>.2f}% [{}/{}] - {}'.format(perc, g+1, ngenes, reg))

    ## compute FDR-correct pvals 
    pvals = data[not_nan][:,1]
    rejected, corrected = fdrcorrection(pvals) 
    np.put(data[:,2], not_nan, corrected)

    ## save results 
    out_path = '{}/{}.txt'.format(out_main, reg) 
    save_results(gene_array, data, out_path)

    ## print summary 
    p_sig = np.sum(data[:,1] <= 0.05)
    f_sig = np.sum(data[:,2] <= 0.05)
    print('> {:>3d} p, {:>3d} f - {}\n'.format(p_sig, f_sig, reg))

pool = Pool(processes=60)
pool.map(find_signif_genes, regions)

