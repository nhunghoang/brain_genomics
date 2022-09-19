'''
Regress confounders (PC1, PC2, age, gender) from 
each gene, as well as each phenotype. 

- Nhung, Tim (updated June 2022) 
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
if dset == 'UKB': expr_dir += '/cohort1365'
if dset == 'HCP': expr_dir += '/cohort890'

phen_out = '/data1/rubinov_lab/brain_genomics/analyses_{}/DATA_OUTPUT/phen_regress'.format(dset)
expr_out = '/data1/rubinov_lab/brain_genomics/analyses_{}/DATA_OUTPUT/r0.3_p0.01/expr_regress'.format(dset)
if not os.path.exists(phen_out): os.mkdir(phen_out)
if not os.path.exists(expr_out): os.mkdir(expr_out)

if dset == 'UKB':
    subj_order1 = '/data1/rubinov_lab/brain_genomics/data_UKB/UKB_subjs_2732.txt' ## covariate order 
    subj_order2 = '/data1/rubinov_lab/brain_genomics/data_UKB/UKB_subjs_1365.txt' ## master order  

## region names 
phen_regions = ['amygdala', 'anterior-cingulate', 'caudate', 'cerebellar-hemisphere', 'frontal-pole', \
                'hippocampus', 'hypothalamus', 'nucleus-accumbens', 'putamen', 'substantia-nigra']

############################################################################################################

## read in covariates 
cov_cols = [-4, -3, -2, -1] ## isFemale, age, PC1, PC2 
COV_DATA = np.loadtxt(cov_file, delimiter='\t', skiprows=1, usecols=cov_cols)

if dset == 'UKB':
    ## subject mask 
    order1 = np.loadtxt(subj_order1) ## (2732,)
    order2 = np.loadtxt(subj_order2) ## (1365,)
    subj_idx = np.where(np.in1d(order1, order2))[0] ## indices of order2 wrt to order1
    assert(np.array_equal(order2, order1[subj_idx]))
    COV_DATA = COV_DATA[subj_idx] 

## function: regress out all covariates 
## X: (n_samples, n_features)
## y: (n_samples) 
def compute_y_residuals(y): 
    lr = LinearRegression().fit(COV_DATA, y) 
    coefs = lr.coef_ ## (n_features) 
    b = lr.intercept_ ## float 
    ccovs = np.matmul(COV_DATA, coefs) ## (n_samples) 
    y_residuals = y - ccovs - b ## (n_samples) 
    return y_residuals  

############################################################################################################

## function: compute and save expression residuals 
def save_expr_residuals():
    for reg in phen_regions: 

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

## function: compute and save phenotype residuals 
def save_phen_residuals():
    phens = ['gm_volume', 'myelination', 'alff', 'reho_noGS', 'connmean_noGS']
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

############################################################################################################

save_expr_residuals()
save_phen_residuals() 
