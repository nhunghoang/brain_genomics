'''
Regress confounders (PCs, age, gender) from 
each gene, as well as each phenotype. 

- Nhung (updated April 2023) 
'''

import numpy as np 
import pandas as pd 
import h5py
from sklearn.linear_model import LinearRegression

## paths 
top_path = '/data1/rubinov_lab/brain_genomics/'

hcp_paths = {}

ukb_paths = { \
    'phen': top_path + 'data_UKB/downloads/phenotypes_40k.csv', \
    'expr': top_path + 'scripts_twas/inputs_UKB/grex_JTI', \
    'covs': top_path + 'scripts_twas/inputs_UKB/covariates.csv', \
    'oute': top_path + 'scripts_twas/inputs_UKB/regress_expr_JTI', \
    'outp': top_path + 'scripts_twas/inputs_UKB/regress_phen', \
} 

## func: regress out all covariates 
## X: (n_samples, n_features)
## y: (n_samples) 
def compute_y_residuals(cov_data, y): 
    lr = LinearRegression().fit(cov_data, y) 
    coefs = lr.coef_ ## (n_features) 
    b = lr.intercept_ ## float 
    ccovs = np.matmul(cov_data, coefs) ## (n_samples) 
    y_residuals = y - ccovs - b ## (n_samples) 
    return y_residuals  

## func: format HCP phenotypes 

## func: get UKB phenotypes and covariates 
def get_UKB_data(paths): 

    regions = ['hippocampus', 'amygdala', 'caudate', 'putamen', \
               'cerebellar-hemisphere', 'nucleus-accumbens']

    leftt_ids = ['26562-2.0', '26563-2.0', '26559-2.0', '26560-2.0', '26557-2.0', '26564-2.0']
    right_ids = ['26593-2.0', '26594-2.0', '26590-2.0', '26591-2.0', '26588-2.0', '26595-2.0']

    keys = {} 
    for i in range(len(regions)): 
        keys[leftt_ids[i]] = 'left_vol_{}'.format(REGS[i])
        keys[right_ids[i]] = 'right_vol_{}'.format(REGS[i])
    
    covs_data = pd.read_csv(paths['covs'])
    eids = covs_data['eid'] 
    covs_vals = covs_data.values[:,1:]

    phen_data = pd.read_csv(paths['phen'], index_col='eid')
    phen_data = phen_data.reindex(eids)
    phen_data.rename(columns=keys, inplace=True)

    return phen_data, covs_vals 

## func: apply covariate regression to phenotypes  
def apply_phen_regr(phen_data, cov_mat, out_path): 

    for i in ['left_vol', 'right_vol']: 
        data = {} ## k: reg, v: subj vector  
        for reg in REGS: 
            col = '{}_{}'.format(i, reg) 
            yres = compute_y_residuals(cov_mat, phen_data[col].values)
            data[reg] = yres
            
        out_file = '{}/{}.hdf5'.format(out_path, i)
        with h5py.File(out_file, 'w') as f: 
            for reg, vals in data.items():
                f[reg] = vals

        print(i)

## func: apply covariate regression to expression 
def apply_expr_regr(expr_path, cov_mat, out_path): 
    
    for reg in REGS: 
        expr_file = '{}/{}.hdf5'.format(expr_path, reg) 
        with h5py.File(expr_file, 'r') as f: 
            genes = f['genes'][()] 
            expr = f['pred_expr'][()] 

        expr_regress = np.zeros_like(expr)
        for g in range(genes.size):
            yres = compute_y_residuals(cov_mat, expr[g])
            expr_regress[g] = yres 

        out_file = '{}/{}.hdf5'.format(out_path, reg)
        with h5py.File(out_file, 'w') as f: 
            (n_genes, n_samples) = expr_regress.shape
            n_gene_chunks = np.min([n_genes, 10])

            ds_preds = f.create_dataset('pred_expr', data=expr_regress, \
                                        chunks=(n_gene_chunks, n_samples), \
                                        dtype=np.dtype('float32'), \
                                        scaleoffset=4, \
                                        compression='gzip')

            ds_genes = f.create_dataset('genes', data=genes, \
                                        dtype=h5py.string_dtype())

        print(reg, 'expr')


phen_data, cov_mat = get_UKB_data(ukb_paths)
#apply_phen_regr(phen_data, cov_mat, ukb_paths['outp'])
apply_expr_regr(ukb_paths['expr'], cov_mat, ukb_paths['oute'])

