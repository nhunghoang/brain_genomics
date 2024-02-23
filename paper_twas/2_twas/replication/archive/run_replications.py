'''
Run multiple UKB replications. 

- Nhung, Feb 2024 
'''

import pandas as pd 
import numpy as np 

import h5py
import sys 
import os

from multiprocessing import Pool 

## params 
disc_sizes = [772, 2000, 5500, 15000]
n_shuffles = np.arange(100)

regs = ['anterior-cingulate', 'dlpfc', 'amygdala', 'hippocampus', \
        'caudate', 'putamen', 'nucleus-accumbens', 'cerebellar-hemisphere'] 

## paths 
main_path = '/data1/rubinov_lab/brain_genomics/paper_twas'
covs_path = f'{main_path}/inputs_UKB/covariates.csv'
phen_path = f'{main_path}/inputs_UKB/phenotypes/vol_mean.csv'
grex_path = f'{main_path}/inputs_UKB/grex_JTI'

outs_path = f'{main_path}/inputs_UKB/replications'

## load covariates 
full_covs = pd.read_csv(covs_path, sep='\t', index_col='eid')
full_subj = full_covs.index.values
subj_arng = np.arange(full_subj.size) 

## load phenotypes 
full_phen = pd.read_csv(phen_path, sep='\t', index_col='eid')

## load grex
full_grex = {}; full_gene = {} 
for reg in regs: 
    gfile = f'{grex_path}/{reg}.hdf5'
    with h5py.File(gfile, 'r') as f: 
        full_grex[reg] = f['pred_expr'][()]
        full_gene[reg] = f['genes'][()]
    print(f'loaded full grex ({reg})')

## func: slice and save the data 
def slice_save(rand_idxs, out_main):

    rand_subj = full_subj[rand_idxs]

    ## out paths
    out_covs = f'{out_main}/covariates.csv'
    out_phen = f'{out_main}/phenotypes/vol_mean.csv'
    out_grex = f'{out_main}/grex_JTI'

    os.mkdir(out_main)
    os.mkdir(out_grex)
    os.mkdir(out_phen[:-13])

    ## slice covs data and save
    rand_covs = full_covs.reindex(rand_subj)
    rand_covs.to_csv(out_covs, sep='\t')

    ## slice phen data and save
    rand_phen = full_phen.reindex(rand_subj)
    rand_phen.to_csv(out_phen, sep='\t')

    ## slice grex data and save
    for reg, grex in full_grex.items(): 
        ofile = f'{out_grex}/{reg}.hdf5'
        rand_grex = grex[:, rand_idxs]
        with h5py.File(ofile, 'w') as f: 
            f['genes'] = full_gene[reg]
            f['pred_expr'] = rand_grex

## func: 
def rand_itr(params): 

    disc_size = params['ds']
    itr = params['it']

    ## discovery data 
    rand_idxs = np.random.choice(subj_arng, disc_size, replace=False)
    out_main = f'{outs_path}/{disc_size}_{itr}'

    if os.path.exists(out_main): return

    slice_save(rand_idxs, out_main)

    ## replication data 
    repl_idxs = np.setdiff1d(subj_arng, rand_idxs)
    out_main = f'{outs_path}/{disc_size}_{itr}c'
    slice_save(repl_idxs, out_main)

    print(disc_size, itr)

## parallel runs  
params = [{'ds': ds, 'it': it} \
           for ds in disc_sizes \
           for it in n_shuffles] 
pool = Pool(processes=20)
pool.map(rand_itr, params)


    

    




