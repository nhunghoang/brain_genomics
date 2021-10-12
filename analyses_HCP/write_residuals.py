'''
- Regress confounders (PC1, PC2, age, gender) from 
  each gene, as well as each phenotype. 

Nhung, Oct 2021
'''

import numpy as np 
from sklearn.linear_model import LinearRegression
import os 
import h5py 

## paths 
cov_file = '/data1/rubinov_lab/brain_genomics/data_HCP/covariates.txt' 
phen_dir = '/data1/rubinov_lab/brain_genomics/data_HCP/hoacer_sn_hth/phenotypes'
expr_dir = '/data1/rubinov_lab/brain_genomics/data_HCP/expression/filtered_quality_r0.3_p0.01/cohort890'

phen_out = '/data1/rubinov_lab/brain_genomics/analyses_HCP/phen_regress'
expr_out = '/data1/rubinov_lab/brain_genomics/analyses_HCP/expr_regress'
if not os.path.exists(phen_out): os.mkdir(phen_out)
if not os.path.exists(expr_out): os.mkdir(expr_out)

## quickly grab regions 
regions = [f.split('.')[0] for f in os.listdir(expr_dir)]

## read in covariates 
## isNative isAsian isAfrican isEuropean isMulti isFemale age PC1 PC2
cov_text = ['isNative', 'isAsian', 'isAfrican', 'isEuropean', 'isMulti', \
            'isFemale', 'age', 'PC1', 'PC2']
cov_cols = [1,2,3,4,5,7,8,9,10]
cov_data = np.loadtxt(cov_file, delimiter='\t', skiprows=1, usecols=cov_cols)
#covs = {cov_text[i]:cov_data[:,i] for i in range(len(cov_text))}

## function: regress out all covariates 
## X: (n_samples, n_features)
## y: (n_samples) 
def compute_y_residuals(y): 
    lr = LinearRegression().fit(cov_data, y) 
    coefs = lr.coef_ ## (n_features) 
    b = lr.intercept_ ## float 
    ccovs = np.matmul(cov_data, coefs) ## (n_samples) 
    y_residuals = y - ccovs - b ## (n_samples) 
    return y_residuals  
    
## loop through the genes of all regions
## compute their residuals and save in the same format  
for reg in regions: 

    ## read original expression data (and genes) 
    with h5py.File('{}/{}.hdf5'.format(expr_dir, reg), 'r') as f: 
        expr_old = np.array(f['pred_expr']) ## (n_genes, n_samples)
        gene_old = np.array(f['genes']) ## (n_genes) 

    ## compute residuals per gene 
    expr_new = np.zeros_like(expr_old) 
    for g in range(gene_old.size): 
        expr_new[g] = compute_y_residuals(expr_old[g])

    ## save new expression data (and original genes) 
    with h5py.File('{}/{}.hdf5'.format(expr_out, reg), 'w') as f: 
        f['pred_expr'] = expr_new 
        f['genes'] = gene_old
        
## loop through the phenotypes of all regions 
## compute their residuals and save in the same format 
for phen in ['alff', 'regional_homogeneity', 'gm_volume']: 

    ## read original phenotype values 
    phen_old = {} ## k: region, v: (n_samples) 
    with h5py.File('{}/{}.hdf5'.format(phen_dir, phen), 'r') as f: 
        for reg in regions: 
            phen_old[reg] = np.array(f[reg])

    ## compute residuals per regional phenotype  
    phen_new = {} ## k: region, v: (n_samples) 
    for reg in regions: 
        phen_new[reg] = compute_y_residuals(phen_old[reg])

    ## save new phenotype values
    with h5py.File('{}/{}.hdf5'.format(phen_out, phen), 'w') as f: 
        for reg in regions: 
            f[reg] = phen_new[reg] 

