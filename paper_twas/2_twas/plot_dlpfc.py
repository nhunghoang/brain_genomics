'''
Scatter plot the two DLPFC TWAS. 

- Nhung, Feb 2024
'''

import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

## paths
main_path = '/Users/nhunghoang/Desktop/remote_platypus/paper_twas'
gtex_path = f'{main_path}/outputs_UKB/twas_JTI/vol_mean/dlpfc.txt'
psyc_path = f'{main_path}/outputs_UKB/twas_JTI/vol_mean/dlpfc_psychencode.txt'

gsig_path = f'{main_path}/outputs_UKB/TWAS_GWAS_sigbars.csv'  
psig_path = f'{main_path}/outputs_UKB/TWAS_GWAS_dlpfcPE_sigbars.csv' 

plot_path = f'{main_path}/outputs_UKB/plots/fig_s5b.pdf'

## load sig data 
gs = pd.read_csv(gsig_path, usecols=['TWAS_FDR', 'TWAS_BON', 'phenotype'])
gs = gs.loc[gs['phenotype'] == 'dlpfc']

ps = pd.read_csv(psig_path, usecols=['TWAS_FDR', 'TWAS_BON', 'phenotype'])
ps = ps.loc[ps['phenotype'] == 'dlpfc_psychencode']

## load twas results 
gg = pd.read_table(gtex_path, usecols=['gene', 'pvalue'])
pp = pd.read_table(psyc_path, usecols=['gene', 'pvalue'])

df = gg.merge(pp, on='gene', how='inner', suffixes=['_gtex', '_psyc'])
df['glog'] = -np.log10(df['pvalue_gtex'])
df['plog'] = -np.log10(df['pvalue_psyc'])

rho = spearmanr(df['glog'], df['plog'])[0]

## plot
plt.ion()
fig, ax = plt.subplots(1, 1, figsize=(3,3))

colr = '#b29fd1'

sc = ax.scatter(df['glog'], df['plog'], c=colr, s=10, alpha=0.4, edgecolor='none')
ec = ax.scatter(df['glog'], df['plog'], c='none', s=10, edgecolor=colr, linewidth=0.5)

ax.set_xlim([0, 16])
ax.set_ylim([0, 16])

ax.set_xticks([0, 4, 8, 12, 16])
ax.set_yticks([0, 4, 8, 12, 16])

ax.vlines(-np.log10(gs['TWAS_FDR']), 0, 16, color='k', linewidth=0.5)
ax.vlines(-np.log10(gs['TWAS_BON']), 0, 16, color='k', linewidth=0.5, linestyle='dashed')

ax.hlines(-np.log10(ps['TWAS_FDR']), 0, 16, color='k', linewidth=0.5)
ax.hlines(-np.log10(ps['TWAS_BON']), 0, 16, color='k', linewidth=0.5, linestyle='dashed')

ax.set_xlabel('DLPFC GTEx: $log_{10}p$ (TWAS)', size=9)
ax.set_ylabel('DLPFC PsychEncode: $log_{10}p$ (TWAS)', size=9)

xleft, xright = ax.get_xlim()
ybottom, ytop = ax.get_ylim()
ax.set_aspect(abs((xright-xleft)/(ybottom-ytop)))

txt = ax.text(10, 14, '$r_s = {:.3f}$'.format(rho))

plt.tight_layout()
plt.savefig(plot_path, format='pdf')


