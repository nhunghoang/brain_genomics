'''
'''

import numpy as np 
import os 
import h5py 

out_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/DATA_OUTPUT/interreg_expr_regress' 
sng_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/interReg_overlap/single_gene_unions_FDR' 

regs = os.listdir('/data1/rubinov_lab/brain_genomics/analyses_HCP/assoc/null_pvals_alff')
phens = ['alff', 'regional_homogeneity', 'gm_volume'] 

for REG in regs: 

    ## final expression matrix (union of phenotype-selected genes) 
    reg_expr = None
    reg_gene = None

    ## final expression matrix (union of independently selected genes) 
    sng_expr = None
    sng_gene = None

    for PHN in phens:

        ## parse single-assoc for (p <= 0.05) gene indices
        pval_file = '/data1/rubinov_lab/brain_genomics/analyses_HCP/assoc/pvals_{}/{}.txt'.format(PHN, REG)
        pval_data = np.loadtxt(pval_file, delimiter='\t', skiprows=1, usecols=[3]) ## FDR
        pval_mask = np.zeros_like(pval_data, dtype=bool)
        pval_mask[pval_data <= 0.05] = True

        ## starting expression matrix (single-assoc genes)
        expr_dir = '/data1/rubinov_lab/brain_genomics/analyses_HCP/DATA_OUTPUT/expr_regress'
        with h5py.File('{}/{}.hdf5'.format(expr_dir, REG), 'r') as f:
            expr_matrix = np.array(f['pred_expr'])[pval_mask] ## (genes * subjects)
            genes_array = np.array(f['genes'])[pval_mask] 

        ## concatenate all independent genes 
        if sng_expr is None: 
            sng_expr = expr_matrix
            sng_gene = genes_array
        else: 
            sng_expr = np.concatenate((sng_expr, expr_matrix)) 
            sng_gene = np.concatenate((sng_gene, genes_array))

        ## parse simann results for weights
        #simann_file = '/data1/rubinov_lab/brain_genomics/analyses_HCP/assoc/simann_{}/{}.log'.format(PHN, REG)
        #weights = np.loadtxt(simann_file, skiprows=3, dtype=bool)

        ### middle expression matrix (multi-assoc genes) 
        #expr_matrix = expr_matrix[weights] 
        #genes_array = genes_array[weights] 

        ### concatenate all phenotype genes 
        #if reg_expr is None: 
        #    reg_expr = expr_matrix
        #    reg_gene = genes_array
        #else: 
        #    reg_expr = np.concatenate((reg_expr, expr_matrix)) 
        #    reg_gene = np.concatenate((reg_gene, genes_array))

    ## save the (unique) union of phenotype genes 
    #uniq_genes, uniq_idx = np.unique(reg_gene, return_index=True)
    #with h5py.File('{}/{}.hdf5'.format(out_path, REG), 'w') as f: 
    #    f['genes'] = uniq_genes 
    #    f['pred_expr'] = reg_expr[uniq_idx]
    #print('{}: {}/{} genes'.format(REG, uniq_genes.size, reg_gene.size)) 

    ## save the (unique) union of independent genes 
    uniq_genes, uniq_idx = np.unique(sng_gene, return_index=True)
    with h5py.File('{}/{}.hdf5'.format(sng_path, REG), 'w') as f: 
        f['genes'] = uniq_genes 
        f['pred_expr'] = sng_expr[uniq_idx]
    print('{}: {}/{} genes'.format(REG, uniq_genes.size, sng_gene.size)) 
