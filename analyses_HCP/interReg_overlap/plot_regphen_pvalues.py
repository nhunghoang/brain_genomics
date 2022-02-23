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
from matplotlib.ticker import FormatStrFormatter

out_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/interReg_overlap/data_matrices'
png_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/interReg_overlap/num_common_genes_pvals-single-regphen-FDR.png'

phens = ['alff', 'regional_homogeneity', 'gm_volume'] 
phens_short = ['ALFF', 'ReHo', 'GMVol']

reg_idx = np.array([5,0,6,9,2,8,7,1,4,3])
regs = os.listdir('/data1/rubinov_lab/brain_genomics/analyses_HCP/assoc/null_pvals_alff')
regs = np.array(sorted(regs)) ## abc order
regs = regs[reg_idx] 

phen_regs = [(phen,reg) for phen in phens for reg in regs] 

## load observed genes from single-gene assoc
reg_genes = {} ## k: reg, v: all regional genes 
obsv_genes = {} ## k: (phen,reg), v: significant genes 
for (phen,reg) in phen_regs: 
    pval_file = '/data1/rubinov_lab/brain_genomics/analyses_HCP/assoc/pvals_{}/{}.txt'.format(phen,reg)
    genes = np.loadtxt(pval_file, delimiter='\t', skiprows=1, usecols=[0], dtype=bytes)
    if reg not in reg_genes.keys(): reg_genes[reg] = genes
    pval_data = np.loadtxt(pval_file, delimiter='\t', skiprows=1, usecols=[3])
    pval_mask = np.zeros_like(pval_data, dtype=bool)
    pval_mask[pval_data <= 0.05] = True
    obsv_genes[(phen,reg)] = genes[pval_mask]

## sort genes by phenotype p-values 
## also record number of significant single genes per permutation 
perm_genes_sort = {} ## k: (phen, reg), v: sorted gene names (perms * genes)  
perm_diag_count = {} ## k: (phen, reg), v: (number of significant genes)
for (PHN, REG) in phen_regs:
    sorted_genes = [] ## (perms * genes) 
    perm_diag = np.zeros(100) 
    for PRM in range(100): 
        pval_file = '/data1/rubinov_lab/brain_genomics/analyses_HCP/assoc/null_pvals_{}/{}/{}.txt'.format(PHN, REG, PRM)
        pval_data = np.loadtxt(pval_file, delimiter='\t', skiprows=1, usecols=[3])
        pval_sort = np.argsort(pval_data) 
        sorted_genes.append(reg_genes[REG][pval_sort])
        perm_diag[PRM] = (pval_data <= 0.05).size
    perm_genes_sort[(PHN,REG)] = np.array(sorted_genes)
    perm_diag_count[(PHN,REG)] = perm_diag

## permutations: take top N genes based on observed regional counts 
perm_genes = {} ## k: (phen,reg), v: gene names (perms * num observed genes)
for (phen,reg) in phen_regs: 
    n_obsv = obsv_genes[(phen,reg)].size
    perm_genes[(phen,reg)] = perm_genes_sort[(phen,reg)][:,:n_obsv]

## count number of interregional gene overlaps per permutation (for top N genes)
mat_perm = np.zeros((100,30,30), dtype=int) ## intersection matrix
for p in range(100):
    for i, (phen1,reg1) in enumerate(phen_regs):
        mat_perm[p][i][i] = perm_genes[(phen1,reg1)][p].size
        for ii, (phen2,reg2) in enumerate(phen_regs[i+1:]):
            j = i + ii + 1
            genes1 = perm_genes[(phen1,reg1)][p]
            genes2 = perm_genes[(phen2,reg2)][p]
            intersect = np.intersect1d(genes1, genes2).size
            mat_perm[p][i][j] = intersect; mat_perm[p][j][i] = intersect

## count observed interregional gene overlap
mat_obsv = np.zeros((30,30), dtype=int) ## intersection matrix
for i, (phen1,reg1) in enumerate(phen_regs):
    mat_obsv[i][i] = obsv_genes[(phen1,reg1)].size
    for ii, (phen2,reg2) in enumerate(phen_regs[i+1:]):
        j = i + ii + 1
        genes1 = obsv_genes[(phen1,reg1)]
        genes2 = obsv_genes[(phen2,reg2)]
        intersect = np.intersect1d(genes1, genes2).size
        mat_obsv[i][j] = intersect; mat_obsv[j][i] = intersect

## heatmap data: effect size
effect_size = mat_obsv / mat_perm.mean(axis=0)

## heatmap data: p-value for regional counts on diagonal, zero otherwise
diag_obsv = np.diag(mat_obsv)
pvals_raw_count = np.ones((30,30)) * 1e-5
for x, (phen,reg) in enumerate(phen_regs):
    pval = (perm_diag_count[(phen,reg)] >= diag_obsv[x]).sum() / 100
    pvals_raw_count[x][x] = pval

## heatmap data: p-value for interregional counts on off-diagonal, zero otherwise
pvals_bootstrap = np.ones((30,30)) * 1e-5
for i, (phen1,reg1) in enumerate(phen_regs):
    for ii, (phen2,reg2) in enumerate(phen_regs[i+1:]):
        j = i + ii + 1
        pval = (mat_perm[:,i,j] >= mat_obsv[i,j]).sum() / 100
        pvals_bootstrap[i][j] = pval; pvals_bootstrap[j][i] = pval

## write matrices to file 
with h5py.File('{}/regphen_pvals_FDR.hdf5'.format(out_path), 'w') as f: 
    f['permutation_overlaps'] = mat_perm 
    f['observed_overlaps'] = mat_obsv

    f['effect_size'] = effect_size 
    f['pvals_raw_count'] = pvals_raw_count 
    f['pvals_bootstrap'] = pvals_bootstrap 

## load data 
with h5py.File('{}/regphen_pvals_FDR.hdf5'.format(out_path), 'r') as f: 
    effect_size = np.array(f['effect_size'])  
    pvals_raw_count = np.array(f['pvals_raw_count']) 
    pvals_bootstrap = np.array(f['pvals_bootstrap']) 

## create figure
fig, axes = plt.subplots(1, 2, figsize=(60,20))
kwargs = {'annot_kws':{'size':15}, 'linewidth':1, 'linecolor':'w'}

## plot effect size
sns.heatmap(effect_size, ax=axes[0], fmt='.2f', annot=True, cmap='Purples', vmin=1, **kwargs)
cbar = axes[0].collections[0].colorbar
cbar.ax.tick_params(labelsize=30)
cbar.ax.set_ylabel('effect size', size=40)
cbar.ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

## plot regional p-values (diagonal) and interregional p-values (off-diagonal)
sns.heatmap(pvals_raw_count, mask=1-np.eye(30), ax=axes[1], \
            fmt='.2f', **kwargs, cbar=False, cmap='Blues_r', annot=False)

#notsig_labels = np.around(pvals_bootstrap, 2).astype(str)
#notsig_labels[pvals_bootstrap < (-1*np.log10(0.05))] = ''
notsig_labels = np.around(pvals_bootstrap, 2).astype(str)
notsig_labels[pvals_bootstrap <= 0.05] = ''
sns.heatmap(pvals_bootstrap, mask=np.eye(30), ax=axes[1],  \
            fmt='', annot=notsig_labels, cmap='Blues_r', cbar=False, \
            annot_kws={'size':15, 'color':'gray'})

#sig_labels = np.around(pvals_bootstrap, 2).astype(str)
#sig_labels[pvals_bootstrap >= (-1*np.log10(0.05))] = ''
sig_labels = np.around(pvals_bootstrap, 2).astype(str)
sig_labels[pvals_bootstrap > 0.05] = ''
sns.heatmap(pvals_bootstrap, mask=np.eye(30), ax=axes[1], \
            fmt='', annot=sig_labels, cmap='Blues_r', \
            annot_kws={'size':18, 'weight':'bold'})

cbar = axes[1].collections[2].colorbar
cbar.ax.set_ylabel('p-value', size=40)
cbar.ax.tick_params(labelsize=30)
cbar.ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

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
plt.savefig('num_common_genes_pvals-single-regphen-FDR.png')
plt.close('all')

