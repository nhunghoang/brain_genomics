'''
Compile all phenotypes.

Nhung, Feb 2024
'''

import numpy as np 
import pandas as pd 
import h5py 

phens = ['Vol', 'ALFF', 'ReHo_GSout', 'FCmean_GSout']
names = ['vol', 'alff', 'reho_noGS', 'connmean_noGS']

## paths 
top_path = '/data1/rubinov_lab/brain_genomics/' 
mat_path = top_path + 'neuro_phenotypes/HCP/HCP_{}_FS.mat' ## Struct or Func
ids_path = top_path + 'paper_twas/inputs_HCP/phenotypes/all_phens_ids.txt'

out_path = top_path + 'paper_twas/inputs_HCP/phenotypes/all_phens.csv'

## fix this
reg_path = top_path + 'paper_twas/0inputs_HCP/archive/mat_files_reg_order.txt'
with open(reg_path, 'r') as f: regs = [r.strip() for r in f.readlines()]

## init dataframe with cohort IDs 
df = pd.read_csv(ids_path, sep='\t')

## read structural mat file
mat_file = mat_path.format('Struct')
with h5py.File(mat_file, 'r') as f: 
    data = f['structure']['Vol'][()] ## (1142 samps, 28 regs)

for i, reg in enumerate(regs): 
    df[f'vol_{reg}'] = data[:,i]

## read functional mat file 
mat_file = mat_path.format('Func')
with h5py.File(mat_file, 'r') as f: 
    data = f['functional']
    
    for phen, name in zip(phens[1:], names[1:]): 
        vals = data[phen][()] 
        for i, reg in enumerate(regs): 
            df[f'{name}_{reg}'] = vals[:,i]

## write to file 
df = df.rename(columns={'SUBJECT_ID': 'subject', 'SAMPLE_ID': 'sample'})
df.to_csv(out_path, index=False)

