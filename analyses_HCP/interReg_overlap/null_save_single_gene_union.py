'''
Diagonal permutation testing
'''

import numpy as np 
import os 
import h5py 

out_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/interReg_overlap/null_single_gene_unions_FDR.hdf5'

regs = os.listdir('/data1/rubinov_lab/brain_genomics/analyses_HCP/assoc/null_pvals_alff')
phens = ['alff', 'regional_homogeneity', 'gm_volume'] 

reg_perm_genes = {} ## k: reg, v: perm_genes
for REG in regs: 

    perm_genes = {} ## k: perm num, v: genes union  
    for PRM in range(100):

        ## final expression matrix (union of independently selected genes)
        reg_gene = None

        for PHN in phens:

            ## parse single-assoc for (p <= 0.05) gene indices
            pval_file = '/data1/rubinov_lab/brain_genomics/analyses_HCP/assoc/null_pvals_{}/{}/{}.txt'.format(PHN, REG, PRM)
            pval_data = np.loadtxt(pval_file, delimiter='\t', skiprows=1, usecols=[3])
            pval_mask = np.zeros_like(pval_data, dtype=bool)
            pval_mask[pval_data <= 0.05] = True

            ## starting gene array (single-assoc genes)
            expr_dir = '/data1/rubinov_lab/brain_genomics/analyses_HCP/DATA_OUTPUT/expr_regress'
            with h5py.File('{}/{}.hdf5'.format(expr_dir, REG), 'r') as f:
                genes_array = np.array(f['genes'])[pval_mask] 

            ## concatenate all independent genes
            if reg_gene is None:
                reg_gene = genes_array
            else:
                reg_gene = np.concatenate((reg_gene, genes_array))

        ## record the (unique) union of phenotype genes 
        uniq_genes, uniq_idx = np.unique(reg_gene, return_index=True)
        perm_genes[PRM] = uniq_genes 

    reg_perm_genes[REG] = perm_genes 

## save all permutation genes
with h5py.File(out_path, 'w') as f:
    for REG in regs:
        perm_genes = reg_perm_genes[REG]
        for p in range(100):
            f['{}-{}'.format(REG,p)] = perm_genes[p] 

