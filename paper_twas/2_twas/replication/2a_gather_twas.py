'''
Gather all the replication TWAS into a readable format (FOR HCP)

- Nhung, Feb 2024
'''

import pandas as pd 
import numpy as np 

from statsmodels.stats.multitest import multipletests as sm

## paths 
main_path = '/data1/rubinov_lab/brain_genomics/paper_twas'
twas_path = f'{main_path}/outputs_HCP/allEuro/subsampling/nonTwin' ## _{i}_{reg}.txt

outs_path = f'{main_path}/outputs_HCP/allEuro/subsampling_twas.csv'

## params 
R = 100 ## num reps

regs = ['dlpfc', 'anterior-cingulate', 'caudate', 'putamen', 'amygdala', \
        'nucleus-accumbens', 'hippocampus', 'cerebellar-hemisphere']

## func: read TWAS and return sig genes  
def read_twas(tpath, disc): 
    tdf = pd.read_table(tpath, sep='\t', usecols=['gene', 'pvalue'])
    if disc: 
        keep = sm(tdf['pvalue'], method='fdr_bh')[0]
    else: 
        keep = tdf['pvalue'] < 0.05
    return tdf['gene'][keep]

## output data: repl, itr, reg, num FDR-sig, num repl
cols = ['repl_set', 'itr', 'vol', 'n_FDR', 'perc_repl']
#nrow = len(repl_sizes) * R * len(regs)
#ncol = len(cols)
#data = np.zeros((nrow, ncol))
data = pd.DataFrame([], columns=cols)

## loop 
loop = [(rs, reg) for rs in repl_sizes for reg in regs]
for (rs0, reg) in loop:
    rs = str(rs0)
    for i in range(R): 
        dd = {'repl_set': rs, 'itr': i, 'vol': reg}

        ## discovery twas 
        dfile = f'{twas_path}/{rs}_{i}_{reg}.txt' 
        dgene = read_twas(dfile, True)

        ## replication twas 
        rfile = f'{twas_path}/{rs}c_{i}_{reg}.txt' 
        rgene = read_twas(rfile, False)

        ## data 
        dd['n_FDR'] = dgene.size 
        if dgene.size == 0: 
            dd['perc_repl'] = 0 
        else:
            dd['perc_repl'] = np.intersect1d(dgene, rgene).size / dgene.size

        data = data.append(dd, ignore_index=True)

        ## HCP if 772
        if rs == '772':
            rfile = f'{hcpp_path}/{reg}.txt' 
            rgene = read_twas(rfile, False)

            dd['repl_set'] = 'HCP'
            dd['perc_repl'] = np.intersect1d(dgene, rgene).size / dgene.size

            data = data.append(dd, ignore_index=True)

    print(rs, reg)

## save 
data.to_csv(outs_path, index=False)


