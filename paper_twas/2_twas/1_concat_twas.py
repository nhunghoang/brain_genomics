'''
Concatenate TWAS results (regional and interregional). 

- Nhung, Feb 2024
'''

import pandas as pd 
import sys

from statsmodels.stats.multitest import multipletests as sm

## params 
group = sys.argv[1] ## UKB 
phens = sys.argv[2] ## vol_mean 

regs = ['dlpfc', 'anterior-cingulate', 'caudate', 'nucleus-accumbens', \
        'putamen', 'amygdala', 'hippocampus', 'cerebellar-hemisphere']

## paths 
main_path = '/data1/rubinov_lab/brain_genomics'
twas_path = f'{main_path}/paper_twas/outputs_{group}/twas_JTI/{phens}'
itwa_path = f'{main_path}/paper_twas/outputs_{group}/twas_JTI/cross_regs/{phens}'

gene_path = f'{main_path}/models_JTI/syms_by_tissue/jti_ens_to_sym.csv'

outs_path = f'{main_path}/paper_twas/outputs_{group}/twas_JTI/summary.csv'

## read ens-to-sym map 
ens2sym = pd.read_csv(gene_path, index_col='ens').to_dict()['sym']

## reg loop 
rlist = [(gr, pr) for gr in regs for pr in regs]
tcols = {'gene': 'ens', 'pvalue': 'pval', 'effect': 'beta'}

dfs = [] 
for (gr, pr) in rlist: 

    if gr == pr: tfile = f'{twas_path}/{gr}.txt'
    else: tfile = f'{itwa_path}/grex_{gr}_phen_{pr}.txt'

    ## gene names
    df = pd.read_table(tfile, usecols=tcols.keys())
    df = df.rename(columns=tcols) 
    df['sym'] = df['ens'].map(ens2sym)

    ## pvalues 
    df['FDR'] = sm(df['pval'], method='fdr_bh', alpha=0.05)[1]
    df['BON'] = sm(df['pval'], method='bonferroni', alpha=0.05)[1]

    ## grex and phen 
    df['grex'] = gr
    df['phen'] = pr

    dfs.append(df)

## save 
df = pd.concat(dfs)
cols = ['sym', 'ens', 'grex', 'phen', 'beta', 'pval', 'FDR', 'BON']
df = df[cols]
df.to_csv(outs_path, index=False)

