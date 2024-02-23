'''
For the UKB TWAS genes that pass FDR < 0.05, 
build a regression model and predict the UKB-repl 
volumes - one model per gene. 

- Nhung, Feb 2024
'''

import numpy as np 
import pandas as pd 
import h5py 

from statsmodels.regression.linear_model import OLS as OLS
from statsmodels.tools.tools import add_constant
from statsmodels.stats.multitest import multipletests as sm

from scipy.stats import pearsonr

regions = ['dlpfc', 'anterior-cingulate', 'amygdala', 'hippocampus', \
           'caudate', 'putamen', 'nucleus-accumbens', 'cerebellar-hemisphere']

## paths 
main_path = '/data1/rubinov_lab/brain_genomics/paper_twas'

utwa_path = f'{main_path}/outputs_UKB/twas_JTI/vol_mean'
ugre_path = f'{main_path}/inputs_UKB/grex_JTI'
uvol_path = f'{main_path}/inputs_UKB/phenotypes/vol_mean.csv'

rtwa_path = f'{main_path}/outputs_UKB/nonBrit/twas_JTI/vol_mean'
rgre_path = f'{main_path}/inputs_UKB/nonBrit/grex_JTI'
rvol_path = f'{main_path}/inputs_UKB/nonBrit/phenotypes/vol_mean.csv'

## get UKB FDR genes 
ugenes = {} ## k: reg, v: [genes] 
for reg in regs: 
    df = pd.read_csv(f'{utwa_path}/{reg}.txt', sep='\t')
    mask = sm(df['pvalue'], method='fdr_bh', alpha=0.05)[0]
    ugenes = df[mask]['gene'].values

## load UKB phen and grex data 
uphen = {} ## k: reg, v: (subj,)
ugrex = {} ## k: reg, v: (FDR_gene, subj)

pdf = pd.read_csv(uvol_path, sep='\t')
for reg in regs: 

    uphen[reg] = pdf[reg].values 
    with h5py.File(f'{ugre_path}/{reg}.hdf5', 'r') as f: 
        glist = f['gene'][()].astype(str)
        grexs = f['pred_expr'][()]

        mask = [np.where(glist == g)[0][0] for g in ugenes[reg]]
        ugrex[reg] = grexs[mask]

## load replication phen and grex data 
rphen = {} ## k: reg, v: (subj,)
rgrex = {} ## k: reg, v: (FDR_gene, subj)

pdf = pd.read_csv(rvol_path, sep='\t')
for reg in regs: 

    rphen[reg] = pdf[reg].values 
    with h5py.File(f'{rgre_path}/{reg}.hdf5', 'r') as f: 
        glist = f['gene'][()].astype(str)
        grexs = f['pred_expr'][()]

        mask = [np.where(glist == g)[0][0] for g in ugenes[reg]]
        rgrex[reg] = grexs[mask]

## apply regression per gene 
rdf = pd.DataFrame([], columns=['reg', 'gene', 'ukb_r2', 'rep_r2', 'ukb_rp', 'rep_rp'])
for reg, glist in ugenes.items(): 
    pdata = uphen[reg]

    for g, gene in enumerate(glist): 
        gdata = ugrex[reg][g] 
    
        ## build regression on UKB 
        gdata = add_constant(gdata) 
        model = OLS(pdata, gdata).fit() 

        ukb_r2 = model.rsquared 
        ppred = model.predict(gdata) 
        ukb_rp = pearsonr(pdata, pphen)[0]

        ## predict repl data 
        rp_data = rphen[reg]
        rg_data = rgrex[reg][g]
        rep_r2 = model.score(rg)Dat
        

        
    

