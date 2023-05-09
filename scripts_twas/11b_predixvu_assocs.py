'''
Aggregate PrediXVU association results for relevant 
genes (i.e., passing the PrediXcan threshold). Save 
as a pandas dataframe for easy parsing. 

- Nhung, Feb 2023 
'''

import pandas as pd
import numpy as np 
import h5py
from statsmodels.stats.multitest import fdrcorrection as FDR

## paths 
pvu_path = '/dors/capra_lab/data/predixVU/byGene'
brn_file = '/data1/rubinov_lab/brain_genomics/scripts_assoc_clean/aux_files/predixvu_brain_phecodes.csv'
sym_file = '/data1/rubinov_lab/brain_genomics/models_PrediXcan_v8/predixvu_ens_sym_map.hdf5'

out_file = '/data1/rubinov_lab/brain_genomics/models_PrediXcan_v8/predixvu_assocs.csv'

## load phecode-to-trait mapping
brain_traits = ['mental disorders', 'neurological']
df = pd.read_csv(brn_file, sep=',')
clinical_map = df.set_index('phecode').to_dict()['phenotype'] ## k: float(phecode), v: clinical trait
brain_phecodes = list(clinical_map.keys())

## load ensembl-to-symbol mapping 
sym_map = {} ## k: ensembl, v: symbol 
with h5py.File(sym_file, 'r') as f: 
    for ens in f.keys(): 
        sym_map[ens] = f[ens][()].decode()

## output PrediXVU dataframe columns
df_data = {'ensembl': [], \
           'phecode': [], \
           'p-value': [], \
           'fdr-val': [], \
           'effsize': [], \
           'tissue': []}

## query PrediXVU for these PrediXcan genes 
for i, [ens, sym] in enumerate(sym_map.items()):
    if (i%20) == 0: print('{}/{} ({} records)'.format(i, len(sym_map), len(df_data['ensembl'])))
    pvu_file = '{}/{}_{}_predixVU.csv.gz'.format(pvu_path, sym, ens) 
    df = pd.read_csv(pvu_file) 
    df = df[df['gene'].notnull()]

    ## consider brain-related phecodes only
    phecodes = df['phecode'].str[1:].astype(float) ## example format: X008.52
    df = df.loc[phecodes.isin(brain_phecodes)]

    ## store this aggregate subset of PrediXVU
    df_data['ensembl'].extend(df['gene'])
    df_data['phecode'].extend(df['phecode'])
    df_data['p-value'].extend(df['pvalue'])
    df_data['effsize'].extend(df['effect_size'])
    df_data['tissue'].extend(df['tissue'])

print('{}/{}'.format(i+1, len(sym_map)))

## adjust p-values for multiple testing hypotheses 
_, fdr = FDR(df_data['p-value'])
df_data['fdr-val'] = fdr

## format phecodes 
codes = [float(i[1:]) for i in df_data['phecode']] ## example format: X008.52
df_data['phecode'] = codes  

## create PrediXVU dataframe 
pvu_df = pd.DataFrame.from_dict(df_data) 
pvu_df.to_csv(out_file)

print('done')
