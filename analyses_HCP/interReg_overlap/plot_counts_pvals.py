'''
'''

import numpy as np
from scipy.stats import pearsonr, spearmanr
import h5py
from time import time
import sys
import os
import seaborn as sns 
import matplotlib.pyplot as plt  
from matplotlib.patches import Rectangle
from matplotlib.ticker import FormatStrFormatter
import bct 

regs = os.listdir('/data1/rubinov_lab/brain_genomics/analyses_HCP/lofo_assoc/train_pvals_alff')
regs = np.array(sorted(regs)) ## abc order 

## gather observed data 
expr_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/interReg_overlap/single_gene_unions_p01'
genes = {} ## k: reg, v: gene list
for reg in regs:
   with h5py.File('{}/{}.hdf5'.format(expr_path, reg), 'r') as f:
       genes[reg] = np.array(f['genes'])

## count observed interregional gene overlap 
count_obsv = np.zeros((10,10), dtype=int) ## intersection matrix
for i, reg1 in enumerate(regs):
    count_obsv[i][i] = genes[reg1].size
    for ii, reg2 in enumerate(regs[i+1:]):
        j = i + ii + 1
        genes1 = genes[reg1]
        genes2 = genes[reg2]
        intersect = np.intersect1d(genes1, genes2).size
        count_obsv[i][j] = intersect; count_obsv[j][i] = intersect

## gather permutation data - regional gene counts   
perm_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/interReg_overlap/null_single_gene_unions_p01.hdf5'
perm_reg_counts = {} ## k: reg, v: (100 counts) 
with h5py.File(perm_path, 'r') as f: 
    for reg in regs: 
        perm_reg_counts[reg] = np.array([np.array(f['{}-{}'.format(reg,p)]).size for p in range(100)]) 

## gather permutation genes in pvalue-sorted order
perm_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/interReg_overlap/null_single_gene_sorted-by-pvals.hdf5'
perm_genes_sorted = {} ## k: reg, v: sorted gene names (perms * genes) 
with h5py.File(perm_path, 'r') as f:
    for reg in regs: 
        perm_genes_sorted[reg] = np.array(f[reg]) 

## permutations: take top N genes based on observed regional counts  
perm_genes = {} ## k: reg, v: (perms * num observed genes) 
for reg in regs: 
    n_obsv = genes[reg].size 
    perm_genes[reg] = perm_genes_sorted[reg][:,:n_obsv]

## count number of interregional gene overlaps per permutation (for top N genes)  
count_perm = np.zeros((100,10,10), dtype=int) ## intersection matrix
for p in range(100):
    for i, reg1 in enumerate(regs):
        count_perm[p][i][i] = perm_genes[reg1][p].size
        for ii, reg2 in enumerate(regs[i+1:]):
            j = i + ii + 1
            genes1 = perm_genes[reg1][p]
            genes2 = perm_genes[reg2][p]
            intersect = np.intersect1d(genes1, genes2).size
            count_perm[p][i][j] = intersect; count_perm[p][j][i] = intersect

## heatmap data: effect size 
effect_size = count_obsv / count_perm.mean(axis=0) 

## heatmap data: p-value for regional counts on diagonal, zero otherwise  
diag_obsv = np.diag(count_obsv) 
pvals_raw_count = np.ones((10,10)) * 1e-5
for r in range(10): 
    pval = (perm_reg_counts[regs[r]] >= diag_obsv[r]).sum() / 100
    #if pval < 1e-6: 
    #    pvals_raw_count[r][r] = 0
    #else: 
    #    pvals_raw_count[r][r] = -1 * np.log10(pval)
    pvals_raw_count[r][r] = pval 

## heatmap data: p-value for interregional counts on off-diagonal, zero otherwise
pvals_bootstrap = np.ones((10,10)) * 1e-5
for i, reg1 in enumerate(regs):
    for ii, reg2 in enumerate(regs[i+1:]):
        j = i + ii + 1
        pval = (count_perm[:,i,j] >= count_obsv[i,j]).sum() / 100 
        #if pval < 1e-6: pval = 0
        #else: pval = -1 * np.log10(pval)
        pvals_bootstrap[i][j] = pval; pvals_bootstrap[j][i] = pval

## manual reordering of regions  
idx = np.array([5,0,6,9,2,8,7,1,4,3])
effect_size = effect_size[idx][:,idx]
pvals_raw_count = pvals_raw_count[idx][:,idx]
pvals_bootstrap = pvals_bootstrap[idx][:,idx]
regs = regs[idx]

## create figure  
fig, axes = plt.subplots(1, 2, figsize=(60,20))
kwargs = {'annot_kws':{'size':35}, 'linewidth':1, 'linecolor':'w'}

## plot effect size 
sns.heatmap(effect_size, ax=axes[0], fmt='.2f', annot=True, cmap='Purples', vmin=1, **kwargs) 
cbar = axes[0].collections[0].colorbar
cbar.ax.tick_params(labelsize=30)
cbar.ax.set_ylabel('effect size', size=40)
cbar.ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

## plot regional p-values (diagonal) and interregional p-values (off-diagonal) 
median_val = np.median(pvals_bootstrap)
sns.heatmap(pvals_raw_count, mask=1-np.eye(10), ax=axes[1], vmax=median_val, \
            fmt='.2f', **kwargs, cbar=False, cmap='Blues_r', annot=False)

#notsig_labels = np.around(pvals_bootstrap, 2).astype(str)
#notsig_labels[pvals_bootstrap < (-1*np.log10(0.05))] = ''
notsig_labels = np.around(pvals_bootstrap, 2).astype(str)
notsig_labels[pvals_bootstrap <= 0.05] = ''
sns.heatmap(pvals_bootstrap, mask=np.eye(10), ax=axes[1], vmax=median_val, \
            fmt='', annot=notsig_labels, cmap='Blues_r', cbar=False, \
            annot_kws={'size':35, 'color':'gray'}) 

#sig_labels = np.around(pvals_bootstrap, 2).astype(str)
#sig_labels[pvals_bootstrap >= (-1*np.log10(0.05))] = ''
sig_labels = np.around(pvals_bootstrap, 2).astype(str)
sig_labels[pvals_bootstrap > 0.05] = ''
sns.heatmap(pvals_bootstrap, mask=np.eye(10), ax=axes[1], vmax=median_val, \
            fmt='', annot=sig_labels, cmap='Blues_r', \
            annot_kws={'size':45, 'weight':'bold'}) 

cbar = axes[1].collections[2].colorbar
cbar.ax.set_ylabel('p-value', size=40)
cbar.ax.tick_params(labelsize=30)
cbar.ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

#for r in range(10):
#    axes[1].add_patch(Rectangle((r, r), 1, 1, fill=False, edgecolor='k', lw=6))

for ax in axes:
    ax.set_xticks(np.arange(10) + 0.5)
    ax.set_yticks(np.arange(10) + 0.5)
    ax.set_xticklabels(regs, rotation=20, ha='right', size=30)
    ax.set_yticklabels(regs, rotation=0, size=30) 
    xleft, xright = ax.get_xlim()
    ybottom, ytop = ax.get_ylim()
    ax.set_aspect(abs((xright-xleft)/(ybottom-ytop)))

plt.tight_layout()
plt.savefig('num_common_genes_pvals_p01.png')
plt.close('all')     

