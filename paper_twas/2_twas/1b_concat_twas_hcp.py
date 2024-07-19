'''
Concat HCP TWAS results (all regional phens). 

- Nhung, March 2024
'''

import pandas as pd 
import numpy as np 
from statsmodels.stats.multitest import multipletests as sm

regs = ['dlpfc', 'anterior-cingulate', 'amygdala', 'hippocampus', \
        'putamen', 'caudate', 'nucleus-accumbens', 'cerebellar-hemisphere']
phens = ['vol_mean', 'alff_mean', 'reho_noGS_mean', 'connmean_noGS_mean']

regphens = [(reg, phen) for reg in regs for phen in phens]

df = None 

main_path = '/data1/rubinov_lab/brain_genomics/paper_twas/outputs_HCP/allEuro/twas_JTI'
outs_path = f'{main_path}/summary.csv'

for (reg, phen) in regphens: 
    tfile = f'{main_path}/{phen}/{reg}.txt'

    tdf = pd.read_table(tfile)
    tdf['reg'] = reg
    tdf['phen'] = phen 
    tdf['FDR'] = sm(tdf['pvalue'], method='fdr_bh', alpha=0.05)[1]

    tdf = tdf[['phen', 'reg', 'gene', 'effect', 'pvalue', 'FDR']]
    if df is None: 
        df = tdf 
    else: 
        df = pd.concat([df, tdf])
    print(reg, phen)

print('')
df.to_csv(outs_path, index=False)
