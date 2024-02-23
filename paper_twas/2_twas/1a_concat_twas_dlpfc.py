'''
Concatenate just the regional DLPFC (GTEx and PsychEncode) TWAS. 

- Nhung, Feb 2024
'''

import pandas as pd 
import sys

from statsmodels.stats.multitest import multipletests as sm

## params 
group = 'UKB' 
phens = 'vol_mean' 

## paths 
main_path = '/data1/rubinov_lab/brain_genomics'
twas_path = f'{main_path}/paper_twas/outputs_{group}/twas_JTI/{phens}'

gene_path = f'{main_path}/models_JTI/syms_by_tissue/jti_ens_to_sym.csv'

outs_path = f'{main_path}/paper_twas/outputs_{group}/twas_JTI/summary_dlpfc.csv'

## read ens-to-sym map 
ens2sym = pd.read_csv(gene_path, index_col='ens').to_dict()['sym']

## reg loop 
regs = {'dlpfc': 'GTEx', 'dlpfc_psychencode': 'PsychEncode'}
tcols = {'gene': 'ens', 'pvalue': 'pval', 'effect': 'beta'}

dfs = [] 
for rname, rmodel in regs.items():

    tfile = f'{twas_path}/{rname}.txt'

    ## gene names
    df = pd.read_table(tfile, usecols=tcols.keys())
    df = df.rename(columns=tcols) 
    df['sym'] = df['ens'].map(ens2sym)

    ## pvalues 
    df['FDR'] = sm(df['pval'], method='fdr_bh', alpha=0.05)[1]
    df['BON'] = sm(df['pval'], method='bonferroni', alpha=0.05)[1]

    dfs.append(df)

## save 
df = dfs[0].merge(dfs[1], on=['sym', 'ens'], how='outer', suffixes=['_GTEx', '_PsychE'])

cols = ['sym', 'ens', 'beta_GTEx', 'pval_GTEx', 'FDR_GTEx', 'BON_GTEx', \
        'beta_PsychE', 'pval_PsychE', 'FDR_PsychE', 'BON_PsychE']

df = df[cols]
df.to_csv(outs_path, index=False)

