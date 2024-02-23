'''
Compare a discovery TWAS against a replication TWAS.

- Nhung, Feb 2024
'''

import pandas as pd 
import numpy as np 
import sys 

from statsmodels.stats.multitest import multipletests as sm
from scipy.stats import pearsonr, spearmanr

repl_twas = sys.argv[1] ## 8020 

## paths 
main_path = '/data1/rubinov_lab/brain_genomics'
#ukbb_path = f'{main_path}/paper_twas/outputs_UKB/rand{repl_twas}/twas_JTI/vol_mean'
#repl_path = f'{main_path}/paper_twas/outputs_UKB/rand{repl_twas}_compl/twas_JTI/vol_mean'

repl_path = f'{main_path}/paper_twas/outputs_UKB/rand{repl_twas}/twas_JTI/vol_mean'
ukbb_path = f'{main_path}/paper_twas/outputs_UKB/rand{repl_twas}_compl/twas_JTI/vol_mean'

## regions 
regs = ['amygdala', 'hippocampus', 'caudate', 'putamen', 'nucleus-accumbens', \
        'cerebellar-hemisphere', 'dlpfc', 'anterior-cingulate']
tcols = ['gene', 'pvalue']

## func: read UKB TWAS 
def read_ukb(reg): 
    upath = f'{ukbb_path}/{reg}.txt'
    df = pd.read_csv(upath, sep='\t', usecols=tcols)
    df['FDR'] = sm(df['pvalue'], method='fdr_bh', alpha=0.05)[1] 

    return df.loc[df['FDR'] < 0.05] 

ukb_twas = {reg: read_ukb(reg) for reg in regs}

## stats 
## repl sig based on FDR and nom (indep) 
## corr b/n nom of uk sig and repl 
## frac of uk sig that is nom < 0.05 in repl 
## (no) frac of uk sig that is FDR < 0.05 in repl 

## loop through UKB TWAS 
for reg, utwas in ukb_twas.items(): 

    print('UKB: {} genes, FDR < 0.05 ({})'.format(utwas.shape[0], reg))

    ## read replication TWAS 
    rpath = f'{repl_path}/{reg}.txt'
    rtwas = pd.read_csv(rpath, sep='\t', usecols=tcols)

    ## get its indep FDR 
    #iFDR = sm(rtwas['pvalue'], method='fdr_bh', alpha=0.05)[0].sum()
    #print('     repl: {} genes, iFDR < 0.05'.format(iFDR))

    ## get its FDR based on UKB 
    df = utwas.merge(rtwas, on='gene', how='inner', suffixes=['_UKBB', '_repl'])
    #uFDR = sm(df['pvalue_repl'], method='fdr_bh', alpha=0.05)[0].sum()
    #print('     repl: {} genes, uFDR < 0.05'.format(uFDR))

    ## get correlation based on UKB genes 
    ukbb_p = -np.log10(df['pvalue_UKBB'].values) 
    repl_p = -np.log10(df['pvalue_repl'].values)
    pr = pearsonr(ukbb_p, repl_p)[0]
    ps = spearmanr(ukbb_p, repl_p)[0]
    print('     rp: {:.2f}, rs: {:.2f}'.format(pr, ps))

    ## get fractions 
    gsize = (df['pvalue_repl'] < 0.05).sum()
    gperc = gsize / utwas.shape[0]
    print('     {} / {} genes ({:.2f}) based on repl nom < 0.05\n'\
          .format(gsize, utwas.shape[0], gperc))


    
    
