'''
Heatmaps (reg phen * reg phen) of the co-expression. 

- Nhung, Feb 2022
'''

import numpy as np
from scipy.stats import pearsonr, spearmanr
import h5py
from time import time
import sys
import os
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

out_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/interReg_overlap/data_matrices'
png_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/interReg_overlap/regphen_coexpression_p05.png'

phens = ['alff', 'regional_homogeneity', 'gm_volume'] 
phens_short = ['ALFF', 'ReHo', 'GMVol']

reg_idx = np.array([5,0,6,9,2,8,7,1,4,3])
regs = os.listdir('/data1/rubinov_lab/brain_genomics/analyses_HCP/assoc/null_pvals_alff')
regs = np.array(sorted(regs)) ## abc order
regs = regs[reg_idx] 

phen_regs = [(phen,reg) for phen in phens for reg in regs] 

## gather data 
data = {} ## k: (phen, reg), v: gene array
expr_data = {} ## k: (phen, reg), v: (genes * subjs)  
for phen in phens:
    for reg in regs:

        ## parse single-assoc for (p <= 0.05) gene indices
        pval_file = '/data1/rubinov_lab/brain_genomics/analyses_HCP/assoc/pvals_{}/{}.txt'.format(phen,reg)
        pval_data = np.loadtxt(pval_file, delimiter='\t', skiprows=1, usecols=[2])
        pval_mask = np.zeros_like(pval_data, dtype=bool)
        pval_mask[pval_data <= 0.05] = True

        ## starting expression matrix (single-assoc genes)
        expr_dir = '/data1/rubinov_lab/brain_genomics/analyses_HCP/DATA_OUTPUT/expr_regress'
        with h5py.File('{}/{}.hdf5'.format(expr_dir, reg), 'r') as f:
            expr_matrix = np.array(f['pred_expr'])[pval_mask] ## (genes * subjects)
            data[(phen,reg)] = np.array(f['genes'])[pval_mask]
            expr_data[(phen,reg)] = expr_matrix 

## format data into matrix 
mat_pearson = np.zeros((30,30), dtype=float) 
mat_spearman = np.zeros((30,30), dtype=float) 

for i, (phen1, reg1) in enumerate(phen_regs):

    mat_pearson[i][i] = 1
    mat_spearman[i][i] = 1

    for ii, (phen2, reg2) in enumerate(phen_regs[i+1:]):
        j = i + ii + 1

        genes1 = data[(phen1,reg1)]
        genes2 = data[(phen2,reg2)]

        ## compute correlation on gene intersection 
        common, idx1, idx2 = np.intersect1d(genes1, genes2, return_indices=True) 
        expr1 = expr_data[(phen1,reg1)][idx1].flatten() 
        expr2 = expr_data[(phen2,reg2)][idx2].flatten()

        rho, _ = pearsonr(expr1, expr2)
        mat_pearson[i][j] = rho; mat_pearson[j][i] = rho 

        if rho == 1: 
            print(reg1, phen1, phen2)

        rho, _ = spearmanr(expr1, expr2)
        mat_spearman[i][j] = rho; mat_spearman[j][i] = rho 
        
with h5py.File('{}/regphen_coexpression_matrices.hdf5'.format(out_path), 'w') as f: 
    f['regphen_pearson'] = mat_pearson 
    f['regphen_spearman'] = mat_spearman 

## load data 
with h5py.File('{}/regphen_coexpression_matrices.hdf5'.format(out_path), 'r') as f: 
    mat_pearson = np.array(f['regphen_pearson']) 
    mat_spearman = np.array(f['regphen_spearman']) 

## plot data as heatmaps 
fig, axes = plt.subplots(1, 2, figsize=(60,20))
kwargs = {'fmt':'.2f', 'cmap':'Blues', 'annot_kws':{'size':15}, \
        'linewidth':1, 'linecolor':'w', 'vmin':0.6, 'vmax':None}

## plot Pearson matrix
sns.heatmap(mat_pearson, annot=True, **kwargs, ax=axes[0])
axes[0].collections[0].colorbar.ax.set_ylabel('pearson', size=40)
axes[0].collections[0].colorbar.ax.tick_params(labelsize=30)
axes[0].collections[0].colorbar.ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

## plot Spearman matrix 
sns.heatmap(mat_spearman, annot=True, **kwargs, ax=axes[1])
axes[1].collections[0].colorbar.ax.set_ylabel('spearman', size=40)
axes[1].collections[0].colorbar.ax.tick_params(labelsize=30)
axes[1].collections[0].colorbar.ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

ax_labels = ['{}-{}'.format(phen,reg[:3]) for phen in phens_short for reg in regs]
for ax in axes:
    ax.set_xticks(np.arange(30) + 0.5)
    ax.set_yticks(np.arange(30) + 0.5)
    ax.set_xticklabels(ax_labels, rotation=20, ha='right', size=15)
    ax.set_yticklabels(ax_labels, rotation=0, size=15)
    xleft, xright = ax.get_xlim()
    ybottom, ytop = ax.get_ylim()
    ax.set_aspect(abs((xright-xleft)/(ybottom-ytop)))

plt.tight_layout()
plt.savefig(png_path)
plt.close('all')
