'''
For a given cohort subset,  
slice the GREx, phenotype, and covariate
files accordingly, to fit the 
PrediXcanAssociation format. 

- Nhung, Feb 2024 
'''

import pandas as pd 
import numpy as np 
import h5py 

import sys 
import os 

from multiprocessing import Pool

group = sys.argv[1] ## HCP 
icoho = sys.argv[2] ## nonTwin
iphen = sys.argv[3] ## vol_mean

## regions 
regs = ['amygdala', 'hippocampus', 'putamen', 'nucleus-accumbens', \
        'caudate', 'cerebellar-hemisphere', 'anterior-cingulate', 'dlpfc']

## paths 
main_path = '/data1/rubinov_lab/brain_genomics/paper_twas'
coho_path = f'{main_path}/inputs_{group}/{icoho}/cohort.txt'
grex_path = f'{main_path}/inputs_{group}/grex_JTI'
phen_path = f'{main_path}/inputs_{group}/phenotypes/{iphen}.csv' 
covs_path = f'{main_path}/inputs_{group}/covariates.csv'

## out paths 
out_main_path = f'{main_path}/inputs_{group}/{icoho}'
out_grex_path = f'{out_main_path}/grex_JTI'
out_phen_path = f'{out_main_path}/phenotypes/{iphen}.csv' 
out_covs_path = f'{out_main_path}/covariates.csv'

if not os.path.exists(out_grex_path):
    os.mkdir(out_grex_path) 
    os.mkdir(f'{out_main_path}/phenotypes')

## load subset cohort 
coho = pd.read_csv(coho_path, sep=' ', header=None)[0].values

## load phen data, slice, and save 
full_phen = pd.read_csv(phen_path, sep='\t', index_col='subject')
coho_phen = full_phen.reindex(coho)
coho_phen.to_csv(out_phen_path, sep='\t')

if os.path.exists(covs_path): 
    sys.exit() 

## load covariate data and get full cohort 
full_covs = pd.read_csv(covs_path, sep='\t', index_col='subject') 
full_subj = full_covs.index.values

## get subset cohort indices wrt to full cohort 
coho_idxs = [np.where(full_subj == c)[0][0] for c in coho] 

## slice covs data and save 
coho_covs = full_covs.reindex(coho)
coho_covs.to_csv(out_covs_path, sep='\t')

## load grex data, slice, and save  
def slice_grex(reg):

    inn_path = f'{grex_path}/{reg}.hdf5'
    with h5py.File(inn_path, 'r') as f: 
        full_gene = f['genes'][()]
        full_grex = f['pred_expr'][()]

    coho_grex = full_grex[:, coho_idxs] 
    out_path = f'{out_grex_path}/{reg}.hdf5'
    with h5py.File(out_path, 'w') as f: 
        f['genes'] = full_gene 
        f['pred_expr'] = coho_grex 

    print(reg)

pool = Pool(processes=len(regs))
pool.map(slice_grex, regs)

    




