'''
Generate Figure 5. 

- Nhung, Dec 2023
'''

import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd
import h5py
import sys

from statsmodels.stats.multitest import multipletests as sm
from scipy.stats import spearmanr

## params 
group = sys.argv[1] ## HCP, HCP/nonTwin, ...

regs = ['dlpfc', 'anterior-cingulate', 'amygdala', 'hippocampus', 
        'caudate', 'putamen', 'nucleus-accumbens', 'cerebellar-hemisphere']

phens = ['vol_mean', 'alff_mean', 'reho_noGS_mean', 'connmean_noGS_mean']
pnames = {'vol_mean': 'volume', 'alff_mean': 'ALFF', \
          'reho_noGS_mean': 'ReHo', 'connmean_noGS_mean': 'mean co-activity'}

regphens = [(r,p) for r in regs for p in phens]

## paths 
main_path = '/Users/nhunghoang/Desktop/remote_platypus/paper_twas'
colr_path = f'{main_path}/aux_files/color_dict.txt'

midd_path = f'{main_path}/outputs_{group}/polygenic_models/median_ytrue_ypred.hdf5'
sing_path = f'{main_path}/outputs_{group}/polygenic_models/single_stats.hdf5'
rsqs_path = f'{main_path}/outputs_{group}/polygenic_models/split_r2s.hdf5'
pval_path = f'{main_path}/outputs_{group}/polygenic_models/split_pvs.hdf5'
null_path = f'{main_path}/outputs_{group}/polygenic_models/nulls.hdf5'

outs_path = f'{main_path}/outputs_{group}/plots/fig5.pdf'

## func: square plots
def square(ax):
    xleft, xright = ax.get_xlim()
    ybottom, ytop = ax.get_ylim()
    ax.set_aspect(abs((xright-xleft)/(ybottom-ytop)))
    return

## regional colors
reg_color = pd.read_table(colr_path, index_col='label').to_dict()['hex']

########################################################################

## load median polygenic stats (subject obsv vs pred phens) 
median_data = {} ## k: (reg,phen), v: [observed, predicted] * subjs
with h5py.File(midd_path, 'r') as f: 
    for key in f.keys(): 
        data = f[key][()]
        [reg, phen] = key.split('X')
        median_data[(reg,phen)] = data 

## load TWAS stats (keep best based on r2)  
sin_r2s = {} ## k: (reg, phen), v: best r2
sin_pvs = {} ## k: (reg, phen), v: FDR of best r2
with h5py.File(sing_path, 'r') as f: 
    for key in f.keys():
        data = f[key][()]
        [reg, phen] = key.split('X')

        ## get best r2
        ibest = np.argmax(data[:,0])
        sin_r2s[(reg,phen)] = data[:,0][ibest]
        sin_pvs[(reg,phen)] = data[:,2][ibest]

## load polygenic stats (r2s and pvals for all splits) 
mul_data = {(r,p): np.zeros((642,2)) for r in regs for p in phens} ## k: (reg, phen), v: split * [r2, pv]
with h5py.File(rsqs_path, 'r') as f: 
    for key in f.keys():  
        [reg, phen] = key.split('X')
        mul_data[(reg,phen)][:,0] = f[key][()]

with h5py.File(pval_path, 'r') as f: 
    for key in f.keys():  
        [reg, phen] = key.split('X')
        mul_data[(reg,phen)][:,1] = f[key][()]

## get the median stats across splits
mul_r2s = {} ## k: (reg, phen), v: median r2
mul_pvs = {} ## k: (reg, phen), v: FDR of median r2

for (reg, phen), data in mul_data.items(): 
    midx = np.argsort(data[:,0])[321]
    mul_r2s[(reg,phen)] = data[:,0][midx]
    mul_pvs[(reg,phen)] = data[:,1][midx]

    print(reg, phen, '{:2f}'.format(mul_r2s[(reg,phen)]))

## FDR correct the (polygenic) median pvals
med_pvs = [mul_pvs[(r,p)] for (r,p) in regphens]
ben = sm(med_pvs, method='fdr_bh', alpha=0.05)[1]
for i, (reg, phen) in enumerate(regphens): 
    mul_pvs[(reg,phen)] = ben[i]

## load polygenic null r2s for all splits  
null_data = {} ## k: (reg, phen), v: flatten(10k nulls, 642 splits)
with h5py.File(null_path, 'r') as f: 
    for key in f.keys(): 
        [reg, phen] = key.split('X')
        null_data[(reg,phen)] = f[key][()].flatten()

########################################################################

## init plot
plt.ion()
fig = plt.figure(figsize=(7,8))
grid = plt.GridSpec(4, 3, figure=fig, width_ratios=[1/5, 3/5, 1/5])

axes = {} 
for a, phen in enumerate(phens):
    axes['multi_scatter_' + phen] = plt.subplot(grid[a,0])
    axes['r2_boxes_' + phen] = plt.subplot(grid[a,1])
    axes['pval_scatter_' + phen] = plt.subplot(grid[a,2])

## func: customize boxplot 
def customize(bp, which):
    if which == 'face':
        _ = [box.set_facecolor(rc) for box,rc in zip(bp['boxes'], rcolors)]
        _ = [box.set_alpha(0.5) for box in bp['boxes']]
        _ = [box.set_edgecolor('none') for box in bp['boxes']]
        _ = [whk.set_color('none') for whk in bp['whiskers']]
        _ = [med.set_color('none') for med in bp['medians']]
    if which == 'edge':
        _ = [box.set_facecolor('none') for box in bp['boxes']]
        _ = [med.set_color('k') for med in bp['medians']]

## main loop
for phen in phens:

    ########################################################################

    ## panel A (median polygenic scatters)
    rselect = {p:r for p,r in zip(phens, \
        ['dlpfc', 'nucleus-accumbens', 'putamen', 'cerebellar-hemisphere'])}
    reg = rselect[phen]

    obsv = median_data[(reg,phen)][:,0]
    pred = median_data[(reg,phen)][:,1]

    ax = axes['multi_scatter_' + phen] 
    sc = ax.scatter(obsv, pred, c=reg_color[reg], s=4, alpha=0.3, edgecolor='none')
    ec = ax.scatter(obsv, pred, c='none', s=4, edgecolor=reg_color[reg], linewidth=0.5)

    if phen == 'vol_mean': 
        xr = [-3000, 0, 3000]
        yr = [-1600, 0, 1600]
    if phen == 'alff_mean': 
        xr = [-50, 0, 50]
        yr = [-20, 0, 20]
    if phen == 'reho_noGS_mean': 
        xr = [-0.003, 0, 0.003]
        yr = [-0.002, 0, 0.002]
    if phen == 'connmean_noGS_mean': 
        xr = [-0.1, 0, 0.1]
        yr = [-0.04, 0, 0.04]

    ax.tick_params(labelsize=8)
    ax.set_xticks(xr)
    ax.set_yticks(yr)

    rr = spearmanr(obsv, pred)[0]
    r2 = mul_r2s[(reg, phen)]
    title = '$r^2$ = {:.2f} ($r$ = {:.2f})'.format(r2, rr)
    ax.set_title(title, size=8)

    square(ax)
    ax.set_xlabel('observed {}'.format(pnames[phen]), size=8)
    ax.set_ylabel('predicted {}'.format(pnames[phen]), size=8)

    ########################################################################

    ## panel B (boxes) 
    positions = np.arange(1, 20, 2.5)
    fly_props = {'marker': '.', 'markersize': 3, 'markeredgecolor': 'none', 'alpha': 1}
    null_colr = '#5e5e5e'

    ax = axes['r2_boxes_' + phen]
    ax.set_xlim([0.5, 20])
    rcolors = [reg_color[reg] for reg in regs]
    
    ## single 
    sin_box = [sin_r2s[(reg,phen)] for reg in regs]
    sc = ax.scatter(positions, sin_box, c=rcolors, s=8, alpha=1, linewidth=1) 

    ## nulls 
    null_box = [null_data[(reg,phen)] for reg in regs] 
    eb = ax.boxplot(null_box, positions=positions+0.5, patch_artist=True, \
                    showcaps=False, showmeans=False, showfliers=False) 
    customize(eb, 'face')

    nb = ax.boxplot(null_box, positions=positions+0.5, patch_artist=True, showfliers=False) 
    customize(nb, 'edge')

    ## multi 
    mul_box = [mul_data[(reg,phen)][:,0] for reg in regs]
    bp = ax.boxplot(mul_box, positions=positions+1, patch_artist=True, \
                    showcaps=False, showmeans=False, showfliers=False) 
    customize(bp, 'face')

    be = ax.boxplot(mul_box, positions=positions+1, patch_artist=True, showfliers=False) 
    customize(be, 'edge')

    ## axis formatting 
    yt = np.arange(0, 0.5, 0.2)
    ax.set_yticks(yt)

    ax.set_xticks(positions + 0.5)
    ax.set_xticklabels('')
    ax.set_ylabel('$r^2$', size=8)

    ########################################################################

    ## panel C (multi vs single pvals)
    for reg in regs: 
        sin_ben = sin_pvs[(reg,phen)]
        sin_ben = -np.log10(sin_ben)
        mul_ben = -np.log10(mul_pvs[(reg,phen)]) 

        ax = axes['pval_scatter_' + phen]
        sc = ax.scatter(sin_ben, mul_ben, c=reg_color[reg], s=8, alpha=0.3, edgecolor='none')
        ec = ax.scatter(sin_ben, mul_ben, c='none', s=8, edgecolor=reg_color[reg], linewidth=0.5)

    ax = axes['pval_scatter_' + phen]
    lp = -np.log10(0.05) - 0.5

    xl = [-0.1, 2]

    if phen == 'vol_mean': 
        yl = [min(0, lp), 11]
        xt = np.arange(0, 3, 1)
        yt = np.arange(0, 11, 5)

    if phen == 'alff_mean': 
        yl = [min(0.5, lp), 8.5]
        xt = np.arange(0, 3, 1)
        yt = np.arange(0, 9, 4)

    if phen == 'reho_noGS_mean': 
        yl = [min(4, lp), 11.5]
        xt = np.arange(0, 3, 1)
        yt = np.arange(0, 11, 5)

    if phen == 'connmean_noGS_mean': 
        yl = [min(2, lp), 8.5]
        xt = np.arange(0, 3, 1)
        yt = np.arange(0, 9, 4)

    ax.set_xlim(xl)
    ax.set_ylim(yl)
    ax.set_xticks(xt)
    ax.set_yticks(yt)

    ax.vlines(-np.log10(0.05), ax.get_ylim()[0], ax.get_ylim()[1], color='k', linewidth=0.3)
    ax.hlines(-np.log10(0.05), ax.get_xlim()[0], ax.get_xlim()[1], color='k', linewidth=0.3)

    ax.margins(0)
    square(ax)

    ax.set_xlabel('best single pval', size=8)
    ax.set_ylabel('median multi pval', size=8)

plt.tight_layout()

## standardize height 
for phen in phens:
    sx, sy, sw, sh = axes['multi_scatter_' + phen].get_position().bounds
    rx, ry, rw, rh = axes['r2_boxes_' + phen].get_position().bounds 
    rc_pos = [rx, sy, rw, sh]
    axes['r2_boxes_' + phen].set_position(rc_pos) 

plt.savefig(outs_path, format='pdf')
