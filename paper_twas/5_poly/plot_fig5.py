'''
Generate Figure 5. 

- Nhung, updated July 2024
'''

import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd
import h5py
import sys

from statsmodels.stats.multitest import multipletests as sm
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
perm_path = f'{main_path}/outputs_{group}/poly_stats_perms.hdf5'

outs_path = f'{main_path}/outputs_HCP/plots/fig5_nonTwin_updated.pdf'

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
grex_sums = {} ## k: (reg, phen), v: (subj,)
phen_arrs = {} ## k: (reg, phen), v: (subj,)
corr_arrs = {} ## k: (reg, phen), v: ([pearson, pval])
num_genes = {} ## k: (reg, phen), v: (1,)
twa_corrs = {} ## k: (reg, phen), v: (1,)

## k: (reg, phen), v: [best single-gene FDR, poly-gene FDR]
fval_arrs = {(reg,phen): np.zeros(2) for (reg, phen) in reg_phens}
pval_arrs = {(reg,phen): np.zeros(2) for (reg, phen) in reg_phens}

with h5py.File(obsv_path, 'r') as f: 
    for (reg,phen) in reg_phens: 
        grex_sums[(reg,phen)] = f[f'{reg}X{phen}Xgrex'][()]
        phen_arrs[(reg,phen)] = f[f'{reg}X{phen}Xphen'][()]
        corr_arrs[(reg,phen)] = f[f'{reg}X{phen}Xcorr'][()]
        num_genes[(reg,phen)] = f[f'{reg}X{phen}Xgene'][()]
        fval_arrs[(reg,phen)] = f[f'{reg}X{phen}XpFDR'][()]
        pval_arrs[(reg,phen)] = f[f'{reg}X{phen}XpNOM'][()]
        twa_corrs[(reg,phen)] = f[f'{reg}X{phen}XrTWA'][()]

## load perm data 
perm_corrs = {} ## k: (reg, phen), v: (perms,)
with h5py.File(perm_path, 'r') as f: 
    for (reg,phen) in reg_phens: 
        perm_corrs[(reg,phen)] = f[f'{reg}X{phen}'][()]

########################################################################

## init plot
plt.ion()
fig = plt.figure(figsize=(7,8))
grid = plt.GridSpec(4, 3, figure=fig, width_ratios=[1/5, 3/5, 1/5])

axes = {} 
for a, phen in enumerate(phens):
    axes['poly_scatter_' + phen] = plt.subplot(grid[a,0])
    axes['corr_boxes_' + phen] = plt.subplot(grid[a,1])
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
        _ = [box.set_linewidth(0.4) for box in bp['boxes']]
        _ = [med.set_color('k') for med in bp['medians']]
        _ = [med.set_linewidth(0.4) for med in bp['medians']]
        _ = [whk.set_linewidth(0.4) for whk in bp['whiskers']]
        _ = [cap.set_linewidth(0.4) for cap in bp['caps']]

## main loop
for phen in phens:

    ########################################################################

    ## panel A (polygenic scatters)
    rselect = {p:r for p,r in zip(phens, \
        ['dlpfc', 'nucleus-accumbens', 'caudate', 'cerebellar-hemisphere'])}
    reg = rselect[phen]

    obsv = grex_sums[(reg,phen)]
    pred = phen_arrs[(reg,phen)]

    ax = axes['poly_scatter_' + phen] 
    sc = ax.scatter(obsv, pred, c=reg_color[reg], s=3.5, alpha=0.1, edgecolor='none')
    ec = ax.scatter(obsv, pred, c='none', s=3.5, alpha=0.5, edgecolor=reg_color[reg], linewidth=0.5)

    if phen == 'vol_mean': 
        yr = [-4500, 0, 4500]
    if phen == 'alff_mean': 
        yr = [-70, 0, 70]
    if phen == 'reho_noGS_mean': 
        yr = [-0.005, 0, 0.010]
    if phen == 'connmean_noGS_mean': 
        yr = [-0.12, 0, 0.12]

    xr = [-1, 0, 1]

    ax.tick_params(labelsize=8)
    ax.set_xticks(xr)
    ax.set_yticks(yr)

    [rho, pval] = corr_arrs[(reg,phen)]
    ngene = num_genes[(reg,phen)]
    if pval == 0: pval = np.array(1/10000)

    title = '$r_p$ = {:.2f} ($p$ = {:.4f})\n{} genes'.format(rho, pval, ngene)
    ax.set_title(title, size=8)

    square(ax)
    #ax.set_xlabel('observed {}'.format(phen), size=8)
    #ax.set_ylabel('predicted {}'.format(phen), size=8)

    ########################################################################

    ## panel B (boxes) 
    positions = np.arange(1, 20, 2.5)
    fly_props = {'marker': '.', 'markersize': 3, 'markeredgecolor': 'none', 'alpha': 1}
    null_colr = '#5e5e5e'

    ax = axes['corr_boxes_' + phen]
    ax.set_xlim([0.5, 20])
    rcolors = [reg_color[reg] for reg in regs]

    ## single 
    sin_box = np.abs([twa_corrs[(reg,phen)] for reg in regs])
    sc = ax.scatter(positions, sin_box, c=rcolors, s=7, alpha=1, linewidth=1) 

    ## nulls 
    null_box = [np.abs(perm_corrs[(reg,phen)]) for reg in regs] 
    eb = ax.boxplot(null_box, positions=positions+0.5, widths=0.4, patch_artist=True, \
                    showcaps=False, showmeans=False, showfliers=False) 
    customize(eb, 'face')

    nb = ax.boxplot(null_box, positions=positions+0.5, widths=0.4, patch_artist=True, showfliers=False) 
    customize(nb, 'edge')

    ## multi 
    mul_box = np.abs([corr_arrs[(reg,phen)][0] for reg in regs])
    mc = ax.scatter(positions+1.1, mul_box, c=rcolors, s=16, marker='_', alpha=1, linewidth=2.5) 
    m2 = ax.scatter(positions+1.3, mul_box, c=rcolors, s=16, marker='_', alpha=1, linewidth=2.5) 

    ## axis formatting 
    ax.set_ylim([0, 0.7]) 
    if phen == phens[2]: ax.set_ylim([0, 0.8])
    yt = np.arange(0, 0.61, 0.2)
    ax.set_yticks(yt)

    ax.set_xlim([0.5, 20.25])
    ax.set_xticks(positions + 0.5)
    ax.set_xticklabels('')
    ax.set_ylabel('abs($r_p$)', size=8)

    ## pval legend
    def p2txt(p): 
        if p < 0.0005: return '***'
        if p < 0.005: return '**' 
        if p < 0.05: return '*'
        return ''

    xpos = [2, 4.5, 7, 9.5, 12, 14.5, 17, 19.5] 
    ypos = mul_box + 0.035

    for r, reg in enumerate(regs): 
        pval = pval_arrs[(reg,phen)][1]
        astr = p2txt(pval) 
        if astr == '': continue 

        ptxt = ax.text(xpos[r], ypos[r], astr, size=7, rotation=90)
        
    
    ########################################################################

    ## panel C (multi vs single pvals)
    colors = [reg_color[reg] for reg in regs]
    sin_fvals = [fval_arrs[(reg,phen)][0] for reg in regs]
    mul_fvals = [fval_arrs[(reg,phen)][1] for reg in regs]

    sin_fvals = -np.log10(sin_fvals)
    mul_fvals = -np.log10(mul_fvals)

    ax = axes['pval_scatter_' + phen]
    sc = ax.scatter(sin_fvals, mul_fvals, c=colors, s=8, alpha=0.3, edgecolor='none')
    ec = ax.scatter(sin_fvals, mul_fvals, c='none', s=8, edgecolor=colors, linewidth=0.5)

    if phen == 'vol_mean': 
        xl = [0, 2]
        yl = [0, 3.7]
        xt = np.arange(4, 6.1, 1)
        yt = np.arange(1, 4, 1)

    if phen == 'alff_mean': 
        xl = [-0.2, 2]
        yl = [-0.2, 2]
        xt = np.arange(3, 6, 1)
        yt = np.arange(0, 2.1, 1)

    if phen == 'reho_noGS_mean': 
        xl = [0, 2]
        yl = [0.4, 4.25]
        xt = np.arange(4, 6, 1)
        yt = np.arange(1, 5, 1)

    if phen == 'connmean_noGS_mean': 
        xl = [0, 1.5]
        yl = [-0.1, 3.6]
        xt = np.arange(4, 6, 1)
        yt = np.arange(0, 4, 1)

    print(ax.get_xlim())
    print(ax.get_ylim())
    ax.set_xlim(xl)
    ax.set_ylim(yl)
    #ax.set_xticks(xt)
    #ax.set_yticks(yt)

    ax.vlines(-np.log10(0.05), ax.get_ylim()[0], ax.get_ylim()[1], color='k', linewidth=0.3)
    ax.hlines(-np.log10(0.05), ax.get_xlim()[0], ax.get_xlim()[1], color='k', linewidth=0.3)

    ax.margins(0)
    square(ax)

    ax.set_xlabel('best single pval (FDR)', size=8)
    ax.set_ylabel('polygenic pval', size=8)

    ########################################################################

plt.tight_layout()

## standardize height 
for phen in phens:
    sx, sy, sw, sh = axes['poly_scatter_' + phen].get_position().bounds
    rx, ry, rw, rh = axes['corr_boxes_' + phen].get_position().bounds 
    rc_pos = [rx, sy, rw, sh]
    axes['corr_boxes_' + phen].set_position(rc_pos) 

plt.savefig(outs_path, format='pdf')
