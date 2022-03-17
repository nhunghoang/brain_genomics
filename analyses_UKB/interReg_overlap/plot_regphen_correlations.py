'''
Heatmaps (reg phen * reg phen) of the phenotype correlations. 

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

out_path = '/data1/rubinov_lab/brain_genomics/analyses_UKB/interReg_overlap'
png_path = '/data1/rubinov_lab/brain_genomics/analyses_UKB/interReg_overlap/regphen_correlations.png'

phens = ['alff', 'regional_homogeneity', 'gm_volume', \
        'connectivity_mean', 'connectivity_variance', \
        'falff', 'gradient', 'myelination', \
        'timeseries_variance', 'fa', 'md']

phens_short = ['ALFF', 'ReHo', 'GMVol', 'ConnMean', 'ConnVar', \
                'FALFF', 'Grad', 'Myel', 'TimeVar', 'FA', 'MD']

phens = phens[:3] 
phens_short = phens_short[:3] 

#idx = [5, 1, 3, 4, 8, 9, 10, 2, 7, 0, 6] 
#phens = np.array(phens)[idx] 
#phens_short = np.array(phens_short)[idx] 

reg_idx = np.array([5,0,6,9,2,8,7,1,4,3])
regs = os.listdir('/data1/rubinov_lab/brain_genomics/analyses_HCP/assoc/null_pvals_alff')
regs = np.array(sorted(regs)) ## abc order
regs = regs[reg_idx] 

phen_regs = [(phen,reg) for phen in phens for reg in regs] 
n_pairs = len(phen_regs) 

## gather data 
data = {} ## k: (phen, reg), v: subject phenotypes  
for phen in phens: 
    phen_path = '/data1/rubinov_lab/brain_genomics/analyses_UKB/DATA_OUTPUT/phen_regress'
    with h5py.File('{}/{}.hdf5'.format(phen_path, phen), 'r') as f: 
        for reg in regs:
            data[(phen,reg)] = np.array(f[reg])

## load subject indices for FA/MD 
#idx_path = '/data1/rubinov_lab/brain_genomics/data_HCP/hoacer_sn_hth/phenotypes/famd_isnan.txt'
#nan_idx = np.loadtxt(idx_path, dtype=int)
#nan_idx = nan_idx.astype(bool)

## format data into matrix 
mat_pearson = np.zeros((n_pairs, n_pairs), dtype=float) 
mat_spearman = np.zeros((n_pairs, n_pairs), dtype=float) 
for i, (phen1, reg1) in enumerate(phen_regs): 
    mat_pearson[i][i] = 0 
    mat_spearman[i][i] = 0 
    data1 = data[(phen1,reg1)]

    for ii, (phen2, reg2) in enumerate(phen_regs[i+1:]):
        j = i + ii + 1
        data2 = data[(phen2,reg2)]

        ## slice subjects in FA/MD case 
        famd = ['fa', 'md'] 
        if (phen1 in famd) and (phen2 not in famd): 
            p_rho, _ = pearsonr(data1, data2[~nan_idx])
            s_rho, _ = spearmanr(data1, data2[~nan_idx])
        elif (phen1 not in famd) and (phen2 in famd): 
            p_rho, _ = pearsonr(data1[~nan_idx], data2)
            s_rho, _ = spearmanr(data1[~nan_idx], data2)
        else:
            p_rho, _ = pearsonr(data1, data2)
            s_rho, _ = spearmanr(data1, data2)

        mat_pearson[i][j] = np.abs(p_rho); mat_pearson[j][i] = np.abs(p_rho) 
        mat_spearman[i][j] = np.abs(s_rho); mat_spearman[j][i] = np.abs(s_rho) 
        
## save data 
with h5py.File('{}/regphen_correlation_matrices.hdf5'.format(out_path), 'w') as f: 
    f['regphen_pearson'] = mat_pearson 
    f['regphen_spearman'] = mat_spearman 

## load data 
with h5py.File('{}/regphen_correlation_matrices.hdf5'.format(out_path), 'r') as f: 
    mat_pearson = np.array(f['regphen_pearson']) 
    mat_spearman = np.array(f['regphen_spearman']) 

#mask = (1 - np.eye(n_pairs)).astype(bool)
#plt.hist(mat_pearson[mask], bins=50) 
#plt.savefig('tmp.png') 

## plot data as heatmaps 
fig, axes = plt.subplots(1, 2, figsize=(60,20))
title_size = 40 
label_size = 30 
ticks_size = 6 
xmin = 0.1
xmax = 0.8
kwargs = {'fmt':'.2f', 'cmap':'Blues', 'annot_kws':{'size':15}, \
         'linewidth':1, 'linecolor':'w', 'vmin':xmin, 'vmax':xmax}

## plot Pearson matrix
sns.heatmap(mat_pearson, annot=False, **kwargs, ax=axes[0])
axes[0].collections[0].colorbar.ax.set_ylabel('pearson', size=title_size)
axes[0].collections[0].colorbar.ax.tick_params(labelsize=label_size)
axes[0].collections[0].colorbar.ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

## plot Spearman matrix
#sns.heatmap(mat_spearman, annot=False, **kwargs, ax=axes[1])
#axes[1].collections[0].colorbar.ax.set_ylabel('spearman', size=title_size)
#axes[1].collections[0].colorbar.ax.tick_params(labelsize=label_size)
#axes[1].collections[0].colorbar.ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

ax_labels = ['{}-{}'.format(phen,reg[:3]) for phen in phens_short for reg in regs]
for ax in axes:
    ax.set_xticks(np.arange(n_pairs) + 0.5)
    ax.set_yticks(np.arange(n_pairs) + 0.5)
    ax.set_xticklabels(ax_labels, rotation=30, ha='right', size=ticks_size)
    ax.set_yticklabels(ax_labels, rotation=0, size=ticks_size)
    xleft, xright = ax.get_xlim()
    ybottom, ytop = ax.get_ylim()
    ax.set_aspect(abs((xright-xleft)/(ybottom-ytop)))

plt.tight_layout()
plt.savefig(png_path)
plt.close('all')
