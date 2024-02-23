'''
Plot S4 (observed vs predicted phenotypes for all reg phens).

- Nhung, Feb 2024
'''

import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

import numpy as np
import pandas as pd
import h5py
import sys

from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests as sm

## params 
group = sys.argv[1] ## HCP, HCP/nonTwin, ... 

regs = ['dlpfc', 'anterior-cingulate', \
        'amygdala', 'hippocampus', 'caudate', 'putamen', \
        'nucleus-accumbens', 'cerebellar-hemisphere']

phens = ['vol_mean', 'alff_mean', 'reho_noGS_mean', 'connmean_noGS_mean']
pnames = {p:n for p,n in zip(phens, ['volume', 'amplitude', 'homogeneity', 'coactivity'])}

reg_phens = [(r,p) for p in phens for r in regs]

## paths 
main_path = f'/Users/nhunghoang/Desktop/remote_platypus/paper_twas'
colr_path = f'{main_path}/aux_files/color_dict.txt'
tick_path = f'{main_path}/aux_files/fig_s4_ticks.txt'

data_path = f'{main_path}/outputs_{group}/polygenic_models/median_ytrue_ypred.hdf5'
stat_path = f'{main_path}/outputs_{group}/polygenic_models/split_r2s.hdf5'
plot_path = f'{main_path}/outputs_HCP/plots/fig_s4.pdf'

## regional colors
reg_color = pd.read_table(colr_path, index_col='label').to_dict()['hex']

## axis ticks 
ticks = pd.read_table(tick_path, index_col='regXphen')

## func: square plots
def square(ax):
    xleft, xright = ax.get_xlim()
    ybottom, ytop = ax.get_ylim()
    ax.set_aspect(abs((xright-xleft)/(ybottom-ytop)))
    return

## load subject data for median multi-models
median_data = {} ## k: (reg,phen), v: [observed, predicted] * subjs
with h5py.File(data_path, 'r') as f: 
    for key in f.keys(): 
        data = f[key][()]
        [reg, phen] = key.split('X')
        median_data[(reg,phen)] = data 

## load median r2s 
med_r2s = {} 
with h5py.File(stat_path, 'r') as f:
    for key in f.keys(): 
        [reg, phen] = key.split('X')
        stats = f[key][()]

        midx = np.argsort(stats)[321]
        med_r2s[(reg,phen)] = stats[midx]

## init plot
plt.ion()
fig, axes_ = plt.subplots(4, 8, figsize=(13.5, 6))
axes = {rp: ax for rp, ax in zip(reg_phens, axes_.flatten())}

for (reg, phen) in reg_phens: 
    ax = axes[(reg,phen)]

    ## get data 
    obsv = median_data[(reg,phen)][:,0]
    pred = median_data[(reg,phen)][:,1]

    sc = ax.scatter(obsv, pred, c=reg_color[reg], s=4, alpha=0.3, edgecolor='none')
    ec = ax.scatter(obsv, pred, c='none', s=4, edgecolor=reg_color[reg], linewidth=0.5)

    ## format ticks 
    key = f'{reg}X{phen}'
    xt = ticks['xt'][key]
    yt = ticks['yt'][key]

    ax.set_xticks([-xt, 0, xt])
    ax.set_yticks([-yt, 0, yt])

    flag = (reg == 'anterior-cingulate') and (phen == 'alff_mean')
    if flag: 
        ax.set_ylim([-1, 5.5])
        ax.set_yticks([0, 2, 4])

    ax.tick_params(size=2, labelsize=8)
    #ax.set_xlabel('observed', size=8)
    #ax.set_ylabel('predicted', size=8)

    square(ax)

    ## labels 
    r2 = med_r2s[(reg, phen)]
    rs = spearmanr(obsv, pred)[0]

    title = '$r^2$ = {:.2f} ($r_s$ = {:.2f})'.format(r2, rs)
    if phen == 'vol_mean': 
        title = f'{reg}\n{title}'
    ax.set_title(title, size=8)

    if reg == 'dlpfc': 
        ax.set_ylabel(pnames[phen], size=8)

plt.tight_layout()
plt.savefig(plot_path, format='pdf')
