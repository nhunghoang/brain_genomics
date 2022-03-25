'''
- Regress confounders (PC1, PC2, age, gender) from 
  each gene, as well as each phenotype. 

Tim, updated Mar 15 2022 
'''

import numpy as np 
from sklearn.linear_model import LinearRegression
import os 
import h5py 
import sys 

dset = sys.argv[1] ## HCP, UKB

## paths 
cov_file = '/data1/rubinov_lab/brain_genomics/analyses_{}/DATA_OUTPUT/covariates.txt'.format(dset)
phen_dir = '/data1/rubinov_lab/brain_genomics/data_{}/hoacer_sn_hth/phenotypes'.format(dset) 
expr_dir = '/data1/rubinov_lab/brain_genomics/data_{}/expression/filtered_quality_r0.3_p0.01'.format(dset)

phen_out = '/data1/rubinov_lab/brain_genomics/analyses_{}/DATA_OUTPUT/phen_regress'.format(dset)
expr_out = '/data1/rubinov_lab/brain_genomics/analyses_{}/DATA_OUTPUT/expr_regress'.format(dset)
if not os.path.exists(phen_out): os.mkdir(phen_out)
if not os.path.exists(expr_out): os.mkdir(expr_out)

## quickly grab regions 
geno_regions = [f.split('.')[0] for f in os.listdir(expr_dir)]
phen_regions = ['amygdala', 'anterior-cingulate', 'caudate', 'cerebellar-hemisphere', 'frontal-pole', 'hippocampus', \
                'hypothalamus', 'nucleus-accumbens', 'putamen', 'substantia-nigra']

## grab genotype subjs for UKB
if dset == 'UKB':
    orig_subjs = np.loadtxt("/data1/rubinov_lab/Yiting/UKB_subject_lists/geno_subjs.txt", dtype=int)
    fin_subjs = np.loadtxt("/data1/rubinov_lab/brain_genomics/data_UKB/ordered_subject_list.txt", dtype=int)

## read in covariates 
## isNative isAsian isAfrican isEuropean isMulti isFemale age PC1 PC2
if dset == 'UKB':
    cov_text = ['isChinese', 'isSouth_Asian', 'isBlack', 'isWhite', 'isMulti', \
            'isFemale', 'age', 'PC1', 'PC2']
elif dset == 'HCP':
    cov_text = ['isNative', 'isAsian', 'isAfrican', 'isEuropean', 'isMulti', \
            'isFemale', 'age', 'PC1', 'PC2']

cov_cols = list(range(1,10))
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
for reg in geno_regions: 

    ## read original expression data (and genes) 
    with h5py.File('{}/{}.hdf5'.format(expr_dir, reg), 'r') as f: 
        expr_old = np.array(f['pred_expr']) ## (n_genes, n_samples)

        if dset == 'UKB':
            [_, idx, _] = np.intersect1d(orig_subjs, fin_subjs, return_indices=True)
            expr_old = expr_old[:, idx]
            
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
if dset == 'HCP':
    phens = ['alff', 'regional_homogeneity', 'gm_volume', \
            'connectivity_mean', 'connectivity_variance', \
            'falff', 'gradient', 'myelination', \
            'timeseries_variance', 'fa', 'md']
elif dset == 'UKB':
    phens = ['alff', 'regional_homogeneity', 'gm_volume', \
            'falff']

for phen in phens: 

    ## read original phenotype values 
    phen_old = {} ## k: region, v: (n_samples) 
    with h5py.File('{}/{}.hdf5'.format(phen_dir, phen), 'r') as f: 
        for reg in phen_regions: 
            phen_old[reg] = np.array(f[reg])

    ## compute residuals per regional phenotype  
    phen_new = {} ## k: region, v: (n_samples) 
    for reg in phen_regions: 
        phen_new[reg] = compute_y_residuals(phen_old[reg])

    ## save new phenotype values
    with h5py.File('{}/{}.hdf5'.format(phen_out, phen), 'w') as f: 
        for reg in phen_regions:
            f[reg] = phen_new[reg]

## diffusion-specific phenotypes (due to missing subjects)
if dset == 'HCP':
    idx_path = '/data1/rubinov_lab/brain_genomics/data_UKB/hoacer_sn_hth/phenotypes/famd_isnan.txt'
    nan_idx = np.loadtxt(idx_path, dtype=int)
    nan_idx = nan_idx.astype(bool)
    for phen in ['fa', 'md']: 

        ## read original phenotype values 
        phen_old = {} ## k: region, v: (n_samples) 
        with h5py.File('{}/{}.hdf5'.format(phen_dir, phen), 'r') as f: 
            for reg in regions: 
                phen_old[reg] = np.array(f[reg])

        ## compute residuals per regional phenotype  
        phen_new = {} ## k: region, v: (n_samples) 
        for reg in regions: 
            y = phen_old[reg] 
            lr = LinearRegression().fit(cov_data[~nan_idx], y)
            coefs = lr.coef_ ## (n_features)
            b = lr.intercept_ ## float
            ccovs = np.matmul(cov_data[~nan_idx], coefs) ## (n_samples)
            y_residuals = y - ccovs - b ## (n_samples)
            phen_new[reg] = y_residuals 

        ## save new phenotype values
        with h5py.File('{}/{}.hdf5'.format(phen_out, phen), 'w') as f: 
            for reg in regions: 
                f[reg] = phen_new[reg] 

