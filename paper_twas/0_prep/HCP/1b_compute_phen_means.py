'''
Compute phenotype means for 
all subjects for all regions. 

 Nhung, Feb 2024
'''

import pandas as pd
import numpy as np 

## paths 
main_path = '/data1/rubinov_lab/brain_genomics/paper_twas/inputs_HCP'
phen_path = f'{main_path}/phenotypes/all_phens.csv'

outs_path = f'{main_path}/phenotypes'

## reg phens 
phens = ['vol', 'alff', 'reho_noGS', 'connmean_noGS']
regs = ['anterior-cingulate', 'dlpfc', 'amygdala', 'hippocampus', \
        'caudate', 'putamen', 'nucleus-accumbens', 'cerebellar-hemisphere']
regphens = [(reg, phen) for reg in regs for phen in phens]

## load phenotypes 
df = pd.read_csv(phen_path)

## reg-phen loop 
for (reg, phen) in regphens: 

    if reg == 'anterior-cingulate': 
        rlist = ['caudal-anterior-cingulate', 'rostral-anterior-cingulate']
    elif reg == 'dlpfc': 
        rlist = ['caudal-middle-frontal', 'rostral-middle-frontal']
    else: 
        rlist = [reg] 

    mean = 0 
    for reg2 in rlist: 
        left = f'{phen}_L_{reg2}'
        righ = f'{phen}_R_{reg2}'
        cols = [left, righ]
        mean += np.mean(df[cols], axis=1)
    mean /= len(rlist)
    df[f'{phen}_{reg}'] = mean

## save file per phen 
for phen in phens: 
    cols = {f'{phen}_{reg}': reg for reg in regs}
    pdf = df[cols.keys()] 
    pdf = pdf.rename(columns=cols)

    pdf.insert(0, 'sample', df['sample'])
    pdf.insert(1, 'subject', df['subject'])
        
    pfile = f'{outs_path}/{phen}_mean.csv'
    pdf.to_csv(pfile, sep='\t', index=False)
        
    
