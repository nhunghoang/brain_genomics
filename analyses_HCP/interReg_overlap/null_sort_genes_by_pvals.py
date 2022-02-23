'''
Off-diagonal permutation testing 

region: 100 * n_genes 
where gene names are sorted by best pvalue across all phenotypes 

'''

import numpy as np 
import os 
import h5py 

out_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/interReg_overlap/null_single_gene_sorted-by-pvals-FDR.hdf5'
expr_dir = '/data1/rubinov_lab/brain_genomics/analyses_HCP/DATA_OUTPUT/expr_regress'

regs = os.listdir('/data1/rubinov_lab/brain_genomics/analyses_HCP/assoc/null_pvals_alff')
phens = ['alff', 'regional_homogeneity', 'gm_volume'] 

## gather regional genes 
genes = {} ## k: reg, v: genes 
for REG in regs: 
    with h5py.File('{}/{}.hdf5'.format(expr_dir, REG), 'r') as f:
        genes[REG] = np.array(f['genes'])

genes_sorted = {} ## k: reg, v: sorted gene names (perms * genes) 
for REG in regs: 
    perm_pvals = np.zeros((100, genes[REG].size, 3)) ## (perms * genes * phenotypes) 

    for PRM in range(100): 
        for pidx, PHN in enumerate(phens):

            ## parse single-assoc for (p <= 0.05) gene indices
            pval_file = '/data1/rubinov_lab/brain_genomics/analyses_HCP/assoc/null_pvals_{}/{}/{}.txt'.format(PHN, REG, PRM)
            pval_data = np.loadtxt(pval_file, delimiter='\t', skiprows=1, usecols=[3])
            perm_pvals[PRM,:,pidx] = pval_data 

    ## for each perm: represent each gene by the lowest p-value across all phenotypes 
    best_pvals = np.min(perm_pvals, axis=2) ## (perms * genes) 

    ## for each perm: sort genes from lowest to highest p-value 
    sorted_idx = np.argsort(best_pvals, axis=1) ## (perms * genes) 
    perm_genes = [genes[REG][order] for order in sorted_idx]
    genes_sorted[REG] = np.array(perm_genes)  

## write sorted genes to file 
with h5py.File(out_path, 'w') as f: 
    for REG in regs: 
        f[REG] = genes_sorted[REG] ## (perms * genes) 
         
