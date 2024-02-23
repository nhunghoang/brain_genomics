'''
Save the complement cohort of  
a given random cohort. 

- Nhung, Feb 2024 
'''

import pandas as pd 
import numpy as np 
import h5py 

import sys 
import os 

from multiprocessing import Pool

group = sys.argv[1] ## UKB 
compl = sys.argv[2] ## rand_400
iphen = sys.argv[3] ## vol_mean

## regions 
regs = ['amygdala', 'hippocampus', 'putamen', 'nucleus-accumbens', \
        'caudate', 'cerebellar-hemisphere', 'anterior-cingulate', 'dlpfc']

## paths 
main_path = '/data1/rubinov_lab/brain_genomics/paper_twas'
grex_path = f'{main_path}/inputs_{group}/grex_JTI'
phen_path = f'{main_path}/inputs_{group}/phenotypes/{iphen}.csv' 
covs_path = f'{main_path}/inputs_{group}/covariates.csv'

comp_path = f'{main_path}/inputs_{group}/{compl}/covariates.csv'

## out paths 
out_main_path = f'{main_path}/inputs_{group}/{compl}_compl'
out_grex_path = f'{out_main_path}/grex_JTI'
out_phen_path = f'{out_main_path}/phenotypes/{iphen}.csv' 
out_covs_path = f'{out_main_path}/covariates.csv'

if not os.path.exists(out_main_path):
    os.mkdir(out_main_path)
    os.mkdir(out_grex_path) 
    os.mkdir(f'{out_main_path}/phenotypes')

## load covariate data and get full cohort 
full_covs = pd.read_csv(covs_path, sep='\t', index_col='eid') 
full_subj = full_covs.index.values

## load existing random cohort and get the complement cohort  
comp_subj = pd.read_csv(comp_path, sep='\t')['eid'].values
rand_subj = np.setdiff1d(full_subj, comp_subj) 
rand_idxs = [np.where(full_subj == s)[0][0] for s in rand_subj]

## slice covs data and save 
rand_covs = full_covs.reindex(rand_subj)
rand_covs.to_csv(out_covs_path, sep='\t')

## load phen data, slice, and save 
full_phen = pd.read_csv(phen_path, sep='\t', index_col='eid')
rand_phen = full_phen.reindex(rand_subj)
rand_phen.to_csv(out_phen_path, sep='\t')

## load grex data, slice, and save  
def slice_grex(reg):

    inn_path = f'{grex_path}/{reg}.hdf5'
    with h5py.File(inn_path, 'r') as f: 
        full_gene = f['genes'][()]
        full_grex = f['pred_expr'][()]

    rand_grex = full_grex[:, rand_idxs] 
    out_path = f'{out_grex_path}/{reg}.hdf5'
    with h5py.File(out_path, 'w') as f: 
        f['genes'] = full_gene 
        f['pred_expr'] = rand_grex 

    print(reg)

pool = Pool(processes=len(regs))
pool.map(slice_grex, regs)

    




