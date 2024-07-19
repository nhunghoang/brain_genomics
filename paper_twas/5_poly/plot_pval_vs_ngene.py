'''

- Nhung, updated July 2024
'''

import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd
import h5py
import sys

from scipy.stats import pearsonr

## params 
group = 'HCP/nonTwin' #sys.argv[1] ## HCP, HCP/nonTwin, ...

regs = ['dlpfc', 'anterior-cingulate', 'amygdala', 'hippocampus', 
        'caudate', 'putamen', 'nucleus-accumbens', 'cerebellar-hemisphere']

phens = ['vol_mean', 'alff_mean', 'reho_noGS_mean', 'connmean_noGS_mean']

reg_phens = [(r,p) for r in regs for p in phens]

## paths 
main_path = '/Users/nhunghoang/Desktop/remote_platypus/paper_twas'
colr_path = f'{main_path}/aux_files/color_dict.txt'

obsv_path = f'{main_path}/outputs_{group}/poly_stats_observed.hdf5'

outs_path = f'{main_path}/outputs_HCP/plots/polygenic_pvals_vs_ngenes.pdf'

## func: square plots
def square(ax):
    xleft, xright = ax.get_xlim()
    ybottom, ytop = ax.get_ylim()
    ax.set_aspect(abs((xright-xleft)/(ybottom-ytop)))
    return

## regional colors
reg_color = pd.read_table(colr_path, index_col='label').to_dict()['hex']

########################################################################

## load observed data 
num_genes = {} ## k: (reg, phen), v: (1,)

## k: (reg, phen), v: (1,)
pval_arrs = {}

with h5py.File(obsv_path, 'r') as f: 
    for (reg,phen) in reg_phens: 
        num_genes[(reg,phen)] = f[f'{reg}X{phen}Xgene'][()]
        pval_arrs[(reg,phen)] = f[f'{reg}X{phen}XpNOM'][()][1]

########################################################################

## init plot
plt.ion()
fig, axes_ = plt.subplots(1, 4, figsize=(6,2))
axes = {phen: ax for phen, ax in zip(phens, axes_.flatten())}

colors = [reg_color[reg] for reg in regs]

## main loop
for phen, ax in axes.items(): 

    ngenes = [num_genes[(reg,phen)] for reg in regs]
    ppvals = [pval_arrs[(reg,phen)] for reg in regs]

    ppvals = -np.log10(ppvals)

    sc = ax.scatter(ngenes, ppvals, c=colors, s=14, alpha=0.3, edgecolor='none')
    ec = ax.scatter(ngenes, ppvals, c='none', s=14, alpha=1.0, edgecolor=colors, linewidth=0.6)

    if phen == 'vol_mean':
        xl = [5, 45]
        yl = [0, 3.7]
        xt = np.arange(10, 45, 10)
        yt = np.arange(0, 4, 1)

    if phen == 'alff_mean':
        xl = [-1, 35]
        yl = [-0.2, 2]
        xt = np.arange(0, 35, 10)
        yt = np.arange(0, 2.1, 1)

    if phen == 'reho_noGS_mean':
        xl = [5, 50]
        yl = [0.4, 4.25]
        xt = np.arange(10, 51, 10)
        yt = np.arange(1, 5, 1)

    if phen == 'connmean_noGS_mean':
        xl = [4, 17]
        yl = [-0.1, 3.6]
        xt = np.arange(5, 17, 5)
        yt = np.arange(0, 4, 1)

    ax.set_xlim([-2, 50])
    ax.set_ylim([-0.2, 4.3])
    ax.set_xticks(np.arange(0, 50, 20))
    ax.set_yticks(np.arange(0, 5, 2))

    ax.tick_params(labelsize=8)

    ax.set_xlabel('num. genes', size=8)
    ax.set_ylabel('-log10(p)', size=8)
    square(ax)

    rho, pval = pearsonr(ngenes, ppvals)
    if pval.round(3) == 0: pval = 0.001
    title = '{}\n$r_p$ = {:.3f} (p < {:.3f})'.format(phen, rho, pval)
    ax.set_title(title, size=8)

plt.tight_layout()
plt.savefig(outs_path, format='pdf')
