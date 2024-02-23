'''
Compare enrichments for TWAS genes passing FDR and Bonferroni. 

- Nhung, Feb 2024
'''

import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt

import pandas as pd
import numpy as np
from scipy.stats import spearmanr

import sys 

## params 
enric = sys.argv[1] ## gwas_catalog, pdx_nom0.001
group = 'UKB'
phens = 'vol_mean_interreg'

## paths 
top_path = f'/Users/nhunghoang/Desktop/remote_platypus/paper_twas/outputs_{group}'
FDR_path = f'{top_path}/enrich_{enric}/FDR_{phens}/enrichment_results'
BON_path = f'{top_path}/enrich_{enric}/BON_{phens}/enrichment_results'

out_path = f'{top_path}/plots/s2b_gwas_catalog.pdf'

## regs 
regs = ['dlpfc', 'anterior_cingulate', 'putamen', 'caudate', 'amygdala', \
        'nucleus_accumbens', 'hippocampus', 'cerebellar_hemisphere']

## load enrichment results 
fs = [] 
bs = [] 
for reg in regs: 
    FDR = pd.read_csv(f'{FDR_path}_{reg}.txt', sep='\t', usecols=['geneSet', 'FDR'])
    BON = pd.read_csv(f'{BON_path}_{reg}.txt', sep='\t', usecols=['geneSet', 'FDR'])
    

    FDR['reg'] = reg; BON['reg'] = reg
    fs.append(FDR); bs.append(BON)

## concat 
fs = pd.concat(fs)
bs = pd.concat(bs)
df = fs.merge(bs, on=['reg', 'geneSet'], how='inner', suffixes=['_tFDR', '_tBON'])

## log pvals 
df['logf'] = -np.log10(df['FDR_tFDR'])
df['logb'] = -np.log10(df['FDR_tBON'])

## correlation 
rho = spearmanr(df['FDR_tFDR'], df['FDR_tBON'])[0]

## plot 
plt.ion()
fig, ax = plt.subplots(1, 1, figsize=(3,3))

sc = ax.scatter(df['logf'], df['logb'], c='k', s=10, alpha=0.1, edgecolor='none', marker='s')
ec = ax.scatter(df['logf'], df['logb'], c='none', s=10, edgecolor='k', marker='s', linewidth=0.5)

if enric == 'gwas_catalog': end = 14
else: end = 11
ax.set_xlim([-0.2, end])
ax.set_ylim([-0.2, end])

ax.set_xticks([0, 5, 10])
ax.set_yticks([0, 5, 10])

ax.vlines(-np.log10(0.05), -0.2, end, color='k', linewidth=0.5)
ax.vlines(-np.log10(0.01), -0.2, end, color='k', linewidth=0.8)

ax.hlines(-np.log10(0.05), -0.2, end, color='k', linewidth=0.5)
ax.hlines(-np.log10(0.01), -0.2, end, color='k', linewidth=0.8)

ax.set_xlabel('-log10(p) (FDR-corrected associations,\n FDR-corrected enrichment)', size=8)
ax.set_ylabel('-log10(p) (Bonf-corrected associations,\n FDR-corrected enrichment)', size=8)

xleft, xright = ax.get_xlim()
ybottom, ytop = ax.get_ylim()
ax.set_aspect(abs((xright-xleft)/(ybottom-ytop)))

txt = ax.text(9, 12.5, '$r_s = {:.3f}$'.format(rho))

plt.tight_layout()
plt.savefig(out_path, format='pdf')
