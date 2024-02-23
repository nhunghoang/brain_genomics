'''
Regress covariates from phenotypes. 

Nhung, updated Feb 2024
'''

import numpy as np 
import pandas as pd 
from sklearn.linear_model import LinearRegression
import sys 

group = sys.argv[1] ## HCP/nonTwin 

## paths 
main_path = '/data1/rubinov_lab/brain_genomics/paper_twas/inputs_{}/'.format(group)
covs_path = main_path + 'covariates.csv'.format(group)
phen_path = main_path + 'phenotypes/{}.csv' ## phen 

outs_path = main_path + 'phenotypes/{}_regr.csv' ## phen 

## func: get residuals after regressing covariates 
def get_residuals(x_covs, y_phen): 
    lr = LinearRegression().fit(x_covs, y_phen)
    y_pred = lr.predict(x_covs) 
    resids = y_phen - y_pred
    return resids 

## read covariates 
ids = ['sample', 'subject']
covs = pd.read_csv(covs_path, sep='\t')
covs = covs.drop(ids, axis=1).values

## phen loop 
phens = ['vol_mean', 'alff_mean', 'reho_noGS_mean', 'connmean_noGS_mean']
for phen in phens: 

    ## get original data, init new data
    pfile = phen_path.format(phen) 
    phen_old = pd.read_csv(pfile, sep='\t')
    phen_new = phen_old[ids].copy()

    ## apply regression 
    regs = phen_old.columns[2:]
    for reg in regs: 
        reg_vals = phen_old[reg].values
        phen_new[reg] = get_residuals(covs, reg_vals)

    ## write to file 
    ofile = outs_path.format(phen) 
    phen_new.to_csv(ofile, index=False, sep='\t')

    
    
    
