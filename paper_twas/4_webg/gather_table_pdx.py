'''
Generate S5 table. 
Columns: phecode, phename, vol, FDR(eratio, FDR), BON(eratio, FDR)

- Nhung, Feb 2024
'''

import pandas as pd 

top_path = '/data1/rubinov_lab/brain_genomics/paper_twas/outputs_UKB' 
FDR_path = f'{top_path}/enrich_pdx_nom0.001/FDR_vol_mean_interreg/enrichment_results' ## _{reg}.txt 
BON_path = f'{top_path}/enrich_pdx_nom0.001/BON_vol_mean_interreg/enrichment_results' 

out_path = f'{top_path}/tables/S5_BioVU.csv'

cols = {'geneSet': 'phecode', 'description': 'phename', \
        'enrichmentRatio': 'enrichment ratio', 'FDR': 'FDR'}

regs = ['dlpfc', 'anterior_cingulate', 'amygdala', 'hippocampus', \
        'caudate', 'putamen', 'nucleus_accumbens', 'cerebellar_hemisphere']

data = None
for reg in regs:
    ff = pd.read_table(f'{FDR_path}_{reg}.txt', usecols=cols.keys()).rename(columns=cols)
    bb = pd.read_table(f'{BON_path}_{reg}.txt', usecols=cols.keys()).rename(columns=cols)

    ff = ff.loc[ff['FDR'] < 0.05]
    bb = bb.loc[bb['FDR'] < 0.05]

    df = ff.merge(bb, on=['phecode', 'phename'], how='outer', suffixes=['_FDR', '_Bonf'])
    df.insert(2, 'region volume', reg)

    if data is None: data = df
    else: data = pd.concat([data, df])

data.to_csv(out_path, index=False)
