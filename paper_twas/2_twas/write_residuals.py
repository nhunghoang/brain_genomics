'''
Regress confounders (PCs, age, gender) from 
each gene, as well as each phenotype. 

- Nhung (updated Sept 2023) 
'''

import numpy as np 
import pandas as pd 
import h5py
import sys 
import os 
from sklearn.linear_model import LinearRegression

group = sys.argv[1]

## paths 
main_path = '/data1/rubinov_lab/brain_genomics/paper_twas/inputs_{}/'.format(group)
phen_path = main_path + 'phenotypes/{}.csv' ## phen
expr_path = main_path + 'grex_JTI/{}.hdf5' ## reg 
covs_path = main_path + 'covariates.csv' 

oute_path = main_path + 'grex_JTI_regr/{}.hdf5'
outp_path = main_path + 'phenotypes_regr/{}.csv' 

##################################################################

## func: regress out all covariates 
## X: (n_samples, n_features)
## y: (n_samples) 
def compute_y_residuals(cov_data, y): 
    lr = LinearRegression().fit(cov_data, y) 
    ypred = lr.predict(cov_data) 
    resid = y - ypred 
    return resid  

##################################################################

## func: regress covariates from grex  
def regress_grex(region, cov_data): 

    ## read original expression 
    rfile = expr_path.format(region) 
    with h5py.File(rfile, 'r') as f: 
        grexx = f['pred_expr'][()]
        genes = f['genes'][()]

    ## apply regression 
    grex_regress = np.zeros_like(grexx)
    for g in range(genes.size): 
        yres = compute_y_residuals(cov_data, grexx[g])
        grex_regress[g] = yres 

    ## write to file 
    ofile = oute_path.format(region)
    with h5py.File(ofile, 'w') as f:
        (n_genes, n_samples) = grex_regress.shape
        n_gene_chunks = np.min([n_genes, 10])

        ds_preds = f.create_dataset('pred_expr', data=grex_regress, \
                                    chunks=(n_gene_chunks, n_samples), \
                                    dtype=np.dtype('float32'), \
                                    scaleoffset=4, \
                                    compression='gzip')

        ds_genes = f.create_dataset('genes', data=genes, \
                                    dtype=h5py.string_dtype())

##################################################################

## func: regress covariates from phenotypes 
def regress_phen(phen, cov_data): 

    ## read original phenotypes 
    pfile = phen_path.format(phen)
    phen_old = pd.read_csv(pfile, sep='\t')

    ## apply regression 
    phen_new = pd.DataFrame().reindex_like(phen_old)

    if group == 'HCP':
        phen_new['sample_id'] = phen_old['sample_id']
        phen_new['subject_id'] = phen_old['subject_id']

        for pname in phen_new.columns[2:]:
            phen_vals = phen_old[pname].values
            phen_new[pname] = compute_y_residuals(cov_data, phen_vals) 

    elif group == 'UKB': 
        phen_new['eid'] = phen_old['eid']

        for pname in phen_new.columns[1:]:
            phen_vals = phen_old[pname].values
            phen_new[pname] = compute_y_residuals(cov_data, phen_vals) 

    ## write to file 
    ofile = outp_path.format(phen)
    phen_new.to_csv(ofile, sep='\t', index=False)

##################################################################

## main function 
def main(): 

    ## read covariates 
    if group == 'HCP': ids = ['sample_id', 'subject_id']
    elif group == 'UKB': ids = ['eid']
    else: ids = None

    covs = pd.read_csv(covs_path, sep='\t')
    covs = covs.drop(ids, axis=1).values

    ## loop through regions and regress expression 
    regions = ['hippocampus', 'amygdala', 'caudate', 'putamen', \
               'nucleus-accumbens', 'anterior-cingulate', 'dlpfc', \
               'dlpfc_psychencode', 'cerebellar-hemisphere']

    #for reg in regions: 
    #    regress_grex(reg, covs)

    ## loop through phenotypes and regress 
    phenotypes = ['vol_mean', 'alff_mean', \
                  'reho_noGS_mean', 'connmean_noGS_mean']
    if group == 'UKB': 
        phenotypes = ['vol_mean']
    
    for phen in phenotypes: 
        regress_phen(phen, covs)
        
main()

