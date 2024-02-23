'''
Save interest gene sets to pass into WebGestalt. 
Note: save both regional and interregional sets.

- Nhung, Feb 2024
'''

import numpy as np
import pandas as pd 
from statsmodels.stats.multitest import multipletests as sm

import sys 

## params 
group = sys.argv[1] ## UKB
phens = sys.argv[2] ## vol_mean 
ptype = sys.argv[3] ## FDR, BON

## regs 
regs = ['dlpfc', 'anterior-cingulate', 'caudate', 'nucleus-accumbens', \
        'putamen', 'amygdala', 'hippocampus', 'cerebellar-hemisphere']

## paths
main_path = '/data1/rubinov_lab/brain_genomics'
twas_path = f'{main_path}/paper_twas/outputs_{group}/twas_JTI/{phens}'
itwa_path = f'{main_path}/paper_twas/outputs_{group}/twas_JTI/cross_regs/{phens}'

outs_path = f'{main_path}/paper_twas/inputs_{group}/enrich_sets'

################################################################################

## regional loop 
if ptype == 'BON': alg = 'bonferroni'
else: alg = 'fdr_bh'

for reg in regs: 
    tfile = f'{twas_path}/{reg}.txt'
    tdata = pd.read_table(tfile, usecols=['gene', 'pvalue'])
    
    fkeep = sm(tdata['pvalue'], method=alg, alpha=0.05)[0]
    genes = tdata['gene'][fkeep]

    ofile = f'{outs_path}/{ptype}_{phens}_{reg}.txt'
    genes.to_csv(ofile, index=False, header=False)

################################################################################

## interregional loop (incl. regional)
glist = {reg: [] for reg in regs} 
rlist = [(gr, pr) for gr in regs for pr in regs]

for (gr, pr) in rlist: 

    if gr == pr: tfile = f'{twas_path}/{gr}.txt'
    else: tfile = f'{itwa_path}/grex_{gr}_phen_{pr}.txt'

    tdata = pd.read_table(tfile, usecols=['gene', 'pvalue'])
    
    fkeep = sm(tdata['pvalue'], method=alg, alpha=0.05)[0]
    genes = tdata['gene'][fkeep].values
 
    glist[pr] = np.union1d(glist[pr], genes)

for reg, genes in glist.items(): 
    ofile = f'{outs_path}/{ptype}_{phens}_interreg_{reg}.txt'
    with open(ofile, 'w') as f: 
        lines = '\n'.join(genes) + '\n'
        f.writelines(lines)
    

        
