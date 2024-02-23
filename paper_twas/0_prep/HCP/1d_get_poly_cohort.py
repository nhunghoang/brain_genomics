'''
HCP polygenic cohort - no twins. 

Nhung, Feb 2024 
'''

import pandas as pd 
import numpy as np 

main_path = '/data1/rubinov_lab/brain_genomics/paper_twas/inputs_HCP'
phen_path = f'{main_path}/phenotypes/all_phens.csv' 
euro_path = f'{main_path}/ancestry.csv'
demo_path = f'{main_path}/demographics.csv'

outs_path = f'{main_path}/nonTwin/cohort.txt'

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

## load twin info 
twin = pd.read_csv(demo_path, index_col='sample').to_dict()['twin_id']
cohort['twin_id'] = cohort['sample'].map(twin)

## keep one random twin per pair  
cohort['keep'] = True
for tid, tdf in cohort.groupby('twin_id'):
    if tdf.shape[0] == 2:
        idxs = tdf.index.values
        drop = np.random.choice(idxs, 1)
        cohort.loc[drop, 'keep'] = False

## sanity check 
mask = cohort['keep']
keep = cohort.loc[mask]
drop = cohort.loc[~mask]

keep_counts = keep['twin_id'].value_counts()
assert((keep_counts == 1).all())
assert(drop['twin_id'].isin(keep['twin_id']).all())

## save cohort 
lines = [f'{s} {s}\n' for s in keep['subject']]
with open(outs_path, 'w') as f: f.writelines(lines)
