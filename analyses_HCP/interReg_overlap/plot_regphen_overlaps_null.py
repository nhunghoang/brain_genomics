'''
Heatmaps (reg phen * reg phen) of significant 
genes from single-gene association. 

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

out_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/interReg_overlap/data_matrices/null_regphen_matrices_p05-FDR.hdf5'
png_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/interReg_overlap/null_regphen_overlap_single_p05-FDR.png'

reg_idx = np.array([5,0,6,9,2,8,7,1,4,3])

phens = ['alff', 'regional_homogeneity', 'gm_volume'] 
phens_short = ['ALFF', 'ReHo', 'GMVol']

regs = os.listdir('/data1/rubinov_lab/brain_genomics/analyses_HCP/assoc/null_pvals_alff')
regs = np.array(sorted(regs)) ## abc order
regs = regs[reg_idx] 

def gather_nulls(perm):
    mat_int = np.zeros((30,30), dtype=int) 
    mat_jac = np.zeros((30,30), dtype=float) 

    ## gather data 
    data = {} ## k: (phen, reg), v: gene array  
    for phen in phens: 
        for reg in regs: 

            ## parse single-assoc for (p <= 0.05) gene indices
            pval_file = '/data1/rubinov_lab/brain_genomics/analyses_HCP/' + \
                'assoc/null_pvals_{}/{}/{}.txt'.format(phen,reg,perm)
            pval_data = np.loadtxt(pval_file, delimiter='\t', skiprows=1, usecols=[2])
            pval_mask = np.zeros_like(pval_data, dtype=bool)
            pval_mask[pval_data <= 0.01] = True

            ## starting expression matrix (single-assoc genes)
            expr_dir = '/data1/rubinov_lab/brain_genomics/analyses_HCP/DATA_OUTPUT/expr_regress'
            with h5py.File('{}/{}.hdf5'.format(expr_dir, reg), 'r') as f:
                #expr_matrix = np.array(f['pred_expr'])[pval_mask] ## (genes * subjects)
                data[(phen,reg)] = np.array(f['genes'])[pval_mask] 

    ## format data into matrices 
    phen_regs = [(phen,reg) for phen in phens for reg in regs] 
    for i, (phen1, reg1) in enumerate(phen_regs): 

        mat_int[i][i] = data[(phen1,reg1)].size
        mat_jac[i][i] = 1 

        for ii, (phen2, reg2) in enumerate(phen_regs[i+1:]):
            j = i + ii + 1

            genes1 = data[(phen1,reg1)]
            genes2 = data[(phen2,reg2)]

            ## compute gene intersection and normalized jaccard
            intersect = np.intersect1d(genes1, genes2).size
            union = np.union1d(genes1, genes2).size

            if union == 0: union = 1e-6

            jaccard = intersect / union
            best_jac = min(genes1.size, genes2.size) / union

            if best_jac == 0: best_jac = 1e-6

            norm_jac = jaccard / best_jac

            mat_int[i][j] = intersect; mat_int[j][i] = intersect 
            mat_jac[i][j] = norm_jac; mat_jac[j][i] = norm_jac 

    return mat_int, mat_jac

#null_int = np.zeros((100,30,30))
#null_jac = np.zeros((100,30,30))
#for n in range(100): 
#    mat_int, mat_jac = gather_nulls(n)
#    null_int[n] = mat_int 
#    null_jac[n] = mat_jac
#    print('{}/100 nulls'.format(n+1))
#
#with h5py.File(out_path, 'w') as f: 
#    f['regphen_intersect'] = null_int 
#    f['regphen_jaccard'] = null_jac

## load data
with h5py.File(out_path, 'r') as f:
    null_int = np.array(f['regphen_intersect'])
    null_jac = np.array(f['regphen_jaccard'])
mat_int = null_int.mean(axis=0)
mat_jac = null_jac.mean(axis=0)

## plot data as heatmaps 
fig, axes = plt.subplots(1, 2, figsize=(60,20))
kwargs_int = {'fmt':'.2f', 'cmap':'Greens', 'annot_kws':{'size':15}, 'linewidth':1, 'linecolor':'w'}
kwargs_jac = {'fmt':'.2f', 'cmap':'Oranges', 'annot_kws':{'size':15}, 'linewidth':1, 'linecolor':'w'}

## plot Jaccardian matrix
sns.heatmap(mat_jac, annot=True, **kwargs_jac, ax=axes[0], vmin=0, vmax=0.05)
axes[0].collections[0].colorbar.ax.set_ylabel('normalized Jaccard similarity', size=40)
axes[0].collections[0].colorbar.ax.tick_params(labelsize=30)
axes[0].collections[0].colorbar.ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

## plot intersection and Jaccard off diagonals
int_mask = np.triu(np.ones((30,30)), k=1).astype(bool) ## lower triangle and diagonal
jac_mask = np.tril(np.ones((30,30)), k=0).astype(bool) ## upper triangle

sns.heatmap(mat_int, annot=True, **kwargs_int, vmin=0, vmax=1)
#sns.heatmap(mat_int, mask=int_mask, annot=True, **kwargs_int, vmin=0, vmax=50)
#sns.heatmap(mat_jac, mask=jac_mask, annot=True, **kwargs_jac, cbar=False, vmin=0, vmax=0.3)
axes[1].collections[0].colorbar.ax.set_ylabel('intersect count', size=40)
axes[1].collections[0].colorbar.ax.tick_params(labelsize=30)

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
