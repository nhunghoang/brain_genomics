'''
Format phenotypes to adhere to the BGENIE style. 

source: https://jmarchini.org/bgenie/

- Nhung, April 2023 (updated Feb 2024)
''' 

import pandas as pd 
import numpy as np 

## paths 
main_path = '/data1/rubinov_lab/brain_genomics/paper_twas/' 
coho_path = main_path + 'inputs_UKB/cohort.txt'
phen_path = main_path + 'inputs_UKB/phenotypes/vol_mean.csv'
covs_path = main_path + 'inputs_UKB/covariates.csv'

phen_out_path = main_path + 'inputs_UKB/phenotypes/gwas_vol_mean.txt' 
covs_out_path = main_path + 'inputs_UKB/gwas_covariates.txt' 

## load cohort, phenotypes, and covariates
coho = pd.read_csv(coho_path, sep=' ', header=None)[0].values

phen_data = pd.read_csv(phen_path, sep='\t', index_col='eid')
covs_data = pd.read_csv(covs_path, sep='\t', index_col='eid')

## sanity check on subject order 
assert(np.array_equal(coho, phen_data.index.values))
assert(np.array_equal(coho, covs_data.index.values))

## format columns 
phen_data.index.rename('FID', inplace=True)
phen_data.insert(0, 'IID', coho)

covs_data.index.rename('FID', inplace=True)
covs_data.insert(0, 'IID', coho)

## save both tables 
phen_data.to_csv(phen_out_path, sep=' ')
covs_data.to_csv(covs_out_path, sep=' ')

