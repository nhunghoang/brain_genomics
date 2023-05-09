'''
Script for parsing HCP/UKB demographic files to 
generate a covariate file in the format needed 
for PrediXcanAssociation. 

Continuous variables are left as is, 
categorical variables are one-hot encoded. 

Columns: subject sample gender age PC1 PC2 ...  
(note that subject and sample IDs are the same IDs for UKB)

- Nhung (rewritten April 2023)  
'''

import csv
import numpy as np 
import pandas as pd 
import h5py
import sys 

## paths 
top_path = '/data1/rubinov_lab/brain_genomics/'

hcp_paths = { \
    'eur': top_path + 'TODO', \
    'dem': top_path + 'data_HCP/subject_demographics/COMPLETE_DEMOGRAPHICS.csv', \
    'gen': top_path + 'data_HCP/subject_demographics/HCP_unrestricted_4_23_2019_21_59_10.csv', \
    'pca': top_path + 'scripts_twas/inputs_HCP/eigen_output', \
    'sub': top_path + 'scripts_twas/inputs_HCP/subj_samp_assoc_order.hdf5', \
    'out': top_path + 'scripts_twas/inputs_HCP/covariates.csv' \
}

ukb_paths = { \
    'eur': top_path + 'scripts_twas/inputs_UKB/cohort_filtered.txt', \
    'dem': top_path + 'data_UKB/downloads/all_covariates.csv', \
    'out': top_path + 'scripts_twas/inputs_UKB/covariates.csv', \
}

## func: parse and format HCP data 
def get_HCP(paths): 
    return

## func: parse and format UKB covariates  
def get_UKB(paths): 

    eids = np.loadtxt(paths['eur'], dtype=int, usecols=[0])
    cols = ['eid', '21003-2.0', '22001-0.0'] + \
           ['22009-0.{}'.format(i) for i in range(1,41)]
    data = pd.read_csv(paths['dem'], index_col='eid', usecols=cols)
    data = data.reindex(eids) 
    data.columns = ['age', 'isMale']  + ['PC{}'.format(i) for i in range(1,41)]

    data.to_csv(paths['out'])

## function calls 
get_UKB(ukb_paths)

