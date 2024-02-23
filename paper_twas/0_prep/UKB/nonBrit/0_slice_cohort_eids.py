'''
From the non-White-British cohort (~5k people with vol and gen data), 
take the subset that self-reported as White. 

- Nhung, Feb 2024

'''

import numpy as np 
import pandas as pd 

top_path = '/data1/rubinov_lab/brain_genomics/paper_twas/inputs_UKB/'
anc_path = top_path + 'ancestry_500k.txt' 
nwb_path = top_path + 'allPops/cohort.txt' 

out_path = top_path + 'nonBrit/cohort.txt' 

####################################################################################

## ancestry of all 500k subjects 
anc = pd.read_csv(anc_path) 

## non White British cohort 
nwb = pd.read_csv(nwb_path, header=None, sep='\t')[0].values 

## get the subset 
mask = (anc['eid'].isin(nwb)) & (anc['top_label'] == 'White')
eids = anc.loc[mask]['eid']

## save subset eids 
lines = [f'{eid} {eid}' for eid in eids]
lines = '\n'.join(lines) + '\n' 

with open(out_path, 'w') as f: f.writelines(lines)

