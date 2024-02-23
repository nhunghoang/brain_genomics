'''
Save PrediXVU ontology based on 
nominal p < 0.001 for WebGestalt. 

- Nhung, Feb 2024
'''

import numpy as np 
import pandas as pd 

## phecodes to remove (non-specific) 
phe_rm = [290.2, 290.3, 306.0, 339.0, 346.0, 346.1, 346.2, 348.0, 348.9, 349.0]
brn_rm = ['Hypothalamus', 'Substantia_nigra', 'Spinal_cord_cervical_c-1']
brn_rm = [f'Brain_{b}_combo' for b in brn_rm] 

## paths 
top_path = '/data1/rubinov_lab/brain_genomics/paper_twas'
pdx_path = f'{top_path}/aux_files/predixvu_assocs.csv.gz'
phe_path = f'{top_path}/aux_files/predixvu_brain_phecodes.csv'

## out paths 
des_path = f'{top_path}/aux_files/pdx_nom0.001.des'
gmt_path = f'{top_path}/aux_files/pdx_nom0.001.gmt' 

## read table and filter 
cols = ['ensembl', 'phecode', 'p-value', 'tissue']
df = pd.read_csv(pdx_path, usecols=cols)

df = df.loc[df['p-value'] < 0.001] 
df = df.loc[df['tissue'].str.contains('Brain')]
df = df.loc[~df['tissue'].isin(brn_rm)]

df = df[['ensembl', 'phecode']].drop_duplicates()
df = df.loc[~df['phecode'].isin(phe_rm)]

## get phecode gene lists and save .gmt  
glists = df.groupby('phecode')['ensembl'].apply(list)

lines = [] 
for phecode, glist in glists.items(): 
    genes = '\t'.join(glist)
    line = f'{phecode}\tNA\t{genes}\n'
    lines.append(line)

with open(gmt_path, 'w') as f: 
    f.writelines(lines)

## save .des 
df = pd.read_csv(phe_path, usecols=['phecode', 'phenotype'])    
df = df.loc[df['phecode'].isin(glists.index.values)]

df.to_csv(des_path, sep='\t', index=False, header=False)
