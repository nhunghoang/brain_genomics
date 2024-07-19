'''
Plot S4 (observed vs predicted phenotypes for all reg phens).

- Nhung, July 2024
'''

import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

import numpy as np
import pandas as pd
import h5py
import sys

from scipy.stats import pearsonr
from statsmodels.stats.multitest import multipletests as sm

## params 
group = 'HCP/nonTwin' #sys.argv[1] ## HCP, HCP/nonTwin, ... 

regs = ['dlpfc', 'anterior-cingulate', 'amygdala', 'hippocampus', \
        'caudate', 'putamen', 'nucleus-accumbens', 'cerebellar-hemisphere']

phens = ['vol_mean', 'alff_mean', 'reho_noGS_mean', 'connmean_noGS_mean']
pnames = {p:n for p,n in zip(phens, ['volume', 'amplitude', 'homogeneity', 'coactivity'])}

reg_phens = [(r,p) for p in phens for r in regs]

## paths 
main_path = f'/Users/nhunghoang/Desktop/remote_platypus/paper_twas'
colr_path = f'{main_path}/aux_files/color_dict.txt'
tick_path = f'{main_path}/aux_files/fig_s4_ticks.txt'

data_path = f'{main_path}/outputs_{group}/poly_stats_observed.hdf5'
plot_path = f'{main_path}/outputs_HCP/plots/fig_s4_nonTwin_updated.pdf'

## regional colors
reg_color = pd.read_table(colr_path, index_col='label').to_dict()['hex']

## axis ticks 
ticks = pd.read_table(tick_path, index_col='regXphen')

##############################################################################

## load data 
grex_sums = {} ## k: (reg, phen), v: (subj,)
phen_arrs = {} ## k: (reg, phen), v: (subj,)
corr_arrs = {} ## k: (reg, phen), v: ([pearson, pval])
num_genes = {} ## k: (reg, phen), v: (1,)

with h5py.File(data_path, 'r') as f:
    for (reg,phen) in reg_phens:
        grex_sums[(reg,phen)] = f[f'{reg}X{phen}Xgrex'][()]
        phen_arrs[(reg,phen)] = f[f'{reg}X{phen}Xphen'][()]
        corr_arrs[(reg,phen)] = f[f'{reg}X{phen}Xcorr'][()]
        num_genes[(reg,phen)] = f[f'{reg}X{phen}Xgene'][()]

##############################################################################

## func: square plots
def square(ax):
    xleft, xright = ax.get_xlim()
    ybottom, ytop = ax.get_ylim()
    ax.set_aspect(abs((xright-xleft)/(ybottom-ytop)))
    return

## init plot
plt.ion()
fig, axes_ = plt.subplots(4, 8, figsize=(13.5, 6))
axes = {rp: ax for rp, ax in zip(reg_phens, axes_.flatten())}

for (reg, phen) in reg_phens: 
    ax = axes[(reg,phen)]

    ## get data 
    obsv = grex_sums[(reg,phen)]
    pred = phen_arrs[(reg,phen)]

    sc = ax.scatter(obsv, pred, c=reg_color[reg], s=3.5, alpha=0.1, edgecolor='none')
    ec = ax.scatter(obsv, pred, c='none', s=3.5, alpha=0.5, edgecolor=reg_color[reg], linewidth=0.5)

    ## format ticks 
    key = f'{reg}X{phen}'
    yt = ticks['yt'][key]

    print(ax.get_xlim())

    if (reg == 'anterior-cingulate') and (phen == 'alff_mean'): 
        ax.set_xticks([-1, 0, 5])
    elif (reg == 'amygdala') and (phen == 'alff_mean'): 
        ax.set_xticks([-2, 0, 1])
    elif (reg == 'nucleus-accumbens') and (phen == 'connmean_noGS_mean'): 
        ax.set_xticks([-2, 0, 1])
    else: 
        ax.set_xticks([-1, 0, 1])

    if (reg == 'hippocampus') and (phen == 'vol_mean'): 
        ax.set_yticks([-2000, 0, 1200])
    elif (reg == 'putamen') and (phen == 'vol_mean'): 
        ax.set_yticks([-1600, 0, 2200])
    elif (reg == 'dlpfc') and (phen == 'alff_mean'): 
        ax.set_yticks([-60, 0, 80])
    elif (reg == 'dlpfc') and (phen == 'reho_noGS_mean'): 
        ax.set_yticks([-0.01, 0, 0.03])
    elif (reg == 'hippocampus') and (phen == 'reho_noGS_mean'): 
        ax.set_yticks([-0.002, 0, 0.014])
    elif (reg == 'caudate') and (phen == 'reho_noGS_mean'): 
        ax.set_yticks([-0.005, 0, 0.01])
    elif (reg == 'cerebellar-hemisphere') and (phen == 'reho_noGS_mean'): 
        ax.set_yticks([-0.002, 0, 0.01])
    else:
        ax.set_yticks([-yt, 0, yt])

    ax.tick_params(size=2, labelsize=8)
    square(ax)

    ## labels 
    [rho, pval] = corr_arrs[(reg,phen)] 
    ngene = num_genes[(reg,phen)]

    if pval == 0: pval = 1/10000

    title = '$r_p$ = {:.2f} ($p$ = {:.4f})\n{} genes'.format(rho, pval, ngene)
    ax.set_title(title, size=8)

    if reg == 'dlpfc': 
        ax.set_ylabel(pnames[phen], size=8)

plt.tight_layout()
plt.savefig(plot_path, format='pdf')
