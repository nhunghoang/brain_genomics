'''
HCP TWAS cohort - just Euro. 

Nhung, Feb 2024 
'''

import pandas as pd 
import numpy as np 

main_path = '/data1/rubinov_lab/brain_genomics/paper_twas/inputs_HCP'
phen_path = f'{main_path}/phenotypes/all_phens.csv' 
euro_path = f'{main_path}/ancestry.csv'
demo_path = f'{main_path}/demographics.csv'

outs_path = f'{main_path}/allEuro/cohort.txt'

## load phenotypes and get valid cohort  
phen = pd.read_csv(phen_path)
cols = phen.columns[2:] 
mask = phen[cols].isna().any(axis=1)
phen_co = phen[~mask][['subject', 'sample']]

## load ancestry and get euro cohort 
euro = pd.read_csv(euro_path)
euro = euro.loc[euro['cluster'] == 'White']
euro_co = euro[['subject', 'sample']]

## merge cohorts 
cohort = phen_co.merge(euro_co, on=['subject', 'sample'], how='inner')

## save cohort 
lines = [f'{s} {s}\n' for s in cohort['subject']]
with open(outs_path, 'w') as f: f.writelines(lines)
