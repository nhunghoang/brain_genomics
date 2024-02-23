'''
We have HCP nonTwin TWAS results and UKB TWAS results for N-nonTwin discovery and replication. 

For every UKB discovery TWAS gene set, get the fraction of HCP and UKB replication genes that match. 

- Nhung, Feb 2024
'''

import pandas as pd 
import numpy as np 

from statsmodels.stats.multitest import multipletests as sm

## paths 
main_path = '/data1/rubinov_lab/brain_genomics/paper_twas/outputs'

disc_path = f'{main_path}_UKB/repl_HCPnonTwin' ## /{itr}_{reg}.txt 
repl_path = f'{main_path}_UKB/repl_HCPnonTwin' ## /repl_{itr}_{reg}.txt

ntwn_path = f'{main_path}_HCP/nonTwin/twas_JTI/vol_mean' ## /{reg}.txt 

outs_path = f'{main_path}_UKB/HCPnonTwin_summary.csv'

## params 
regs = ['dlpfc', 'anterior-cingulate', 'caudate', 'putamen', 'amygdala', \
        'hippocampus', 'nucleus-accumbens', 'cerebellar-hemisphere']

n_itrs = 300

## parse twas results
disc = {} ## k: (reg, itr), v: list of FDR genes 
repl = {} ## k: (reg, itr), v: list of p 0.05 genes 
ntwn = {} ## k: reg, v: list of p 0.05 genes

for reg in regs: 

    nfile = f'{ntwn_path}/{reg}.txt'
    df = pd.read_table(nfile, usecols=['gene', 'pvalue'])
    mask = df['pvalue'] < 0.05
    ntwn[reg] = df['gene'][mask].values

    for itr in range(n_itrs):
        dfile = f'{disc_path}/{itr}_{reg}.txt'
        df = pd.read_table(dfile, usecols=['gene', 'pvalue'])
        mask = sm(df['pvalue'], method='fdr_bh', alpha=0.05)[0]
        disc[(reg, itr)] = df['gene'][mask].values

        rfile = f'{repl_path}/repl_{itr}_{reg}.txt'
        df = pd.read_table(rfile, usecols=['gene', 'pvalue'])
        mask = df['pvalue'] < 0.05
        repl[(reg, itr)] = df['gene'][mask].values

## get replication stats
stats = pd.DataFrame([], columns=['reg', 'itr', 'n_disc', 'ukb_repl', 'hcp_nontwin'])

for reg in regs: 
    ngenes = ntwn[reg] 

    for i in range(n_itrs): 
        dgenes = disc[(reg, i)]
        rgenes = repl[(reg, i)]

        ds = dgenes.size
        hcp = np.intersect1d(dgenes, ngenes).size / ds
        ukb = np.intersect1d(dgenes, rgenes).size / ds

        data = {'reg': reg, 'itr': i, 'n_disc': ds, 'ukb_repl': ukb, 'hcp_nontwin': hcp}
        stats = stats.append(data, ignore_index=True)

stats.to_csv(outs_path, index=False)
        
        
 




