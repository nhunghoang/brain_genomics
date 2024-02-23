'''
Generate S3 table. 
Columns: annotation, vol, FDR(eratio, FDR), BON(eratio, FDR)

Also generate a full summary (no FDR filter). 

- Nhung, Feb 2024
'''

import pandas as pd 

top_path = '/data1/rubinov_lab/brain_genomics/paper_twas' 
FDR_path = f'{top_path}/outputs_UKB/enrich_gwas_catalog/FDR_vol_mean/enrichment_results' ## _{reg}.txt 
BON_path = f'{top_path}/outputs_UKB/enrich_gwas_catalog/BON_vol_mean/enrichment_results' 

brn_path = f'{top_path}/aux_files/brain_relatedness_gwas_catalog.csv'

out_path = f'{top_path}/outputs_UKB/tables/S3_GWAS_Catalog.csv'
sum_path = f'{top_path}/outputs_UKB/enrich_summary_GWAS_Catalog.csv'

cols = {'geneSet': 'annotation', 'enrichmentRatio': 'enrichment ratio', 'FDR': 'FDR'}

regs = ['dlpfc', 'anterior_cingulate', 'amygdala', 'hippocampus', \
        'caudate', 'putamen', 'nucleus_accumbens', 'cerebellar_hemisphere']

brain_rel = pd.read_csv(brn_path, index_col='annotation').to_dict()['brain_related']

data = None
summ = None
for reg in regs:
    ff = pd.read_table(f'{FDR_path}_{reg}.txt', usecols=cols.keys()).rename(columns=cols)
    bb = pd.read_table(f'{BON_path}_{reg}.txt', usecols=cols.keys()).rename(columns=cols)

    df = ff.merge(bb, on='annotation', how='outer', suffixes=['_FDR', '_Bonf'])
    df['annotation'] = df['annotation'].apply(lambda x: x.split(' | ')[0])
    df.insert(1, 'region volume', reg)

    br = df['annotation'].map(brain_rel)
    df.insert(0, 'brain_related', br)

    if summ is None: summ = df
    else: summ = pd.concat([summ, df])

    mask = (df['FDR_FDR'] < 0.05) | (df['FDR_Bonf'] < 0.05)
    df = df.loc[mask]

    if data is None: data = df
    else: data = pd.concat([data, df])

summ.to_csv(sum_path, index=False)
data.to_csv(out_path, index=False)
