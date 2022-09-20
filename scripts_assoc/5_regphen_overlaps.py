'''
Aggregate genes of interest (i.e. significant genes) 
for every region and every phenotype. Use this data to 
analyze similarity of selected genes for every region 
pair and every phenotype pair. 

- Nhung, updated Sept 2022
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
from matplotlib.colors import ListedColormap

## paths 
path_assoc = '/data1/rubinov_lab/brain_genomics/analyses_HCP/assoc_1M'
path_out = '/data1/rubinov_lab/brain_genomics/analyses_HCP/assoc_1M'

## names 
nperms = 10000
phens = ['gm_volume', 'myelination', 'alff', 'reho_noGS', 'connmean_noGS']

regs = ['frontal-pole', 'anterior-cingulate', 'caudate', 'putamen', 'nucleus-accumbens', \
        'hippocampus', 'amygdala', 'hypothalamus', 'substantia-nigra', 'cerebellar-hemisphere']
reg_names = ['frontal pole', 'ant. cingulate', 'caudate', 'putamen', 'n. accumbens', \
             'hippocampus', 'amygdala', 'hypothalamus', 's. nigra', 'cerebellum'] 
set_type = {'phens': phens, 'regs': regs} 

## parse association data 
regphen_genes = {k:[] for k in regs+phens} ## k: reg or phen, v: unique gene array 
gene_counts = {} ## k: (reg, phen), v: number of sig genes 
reg_genes = {} ## k: reg, v: gene array 
for phn in phens: 
    for reg in regs: 
        afile = '{}/pvals_{}/{}_1M.hdf5'.format(path_assoc, phn, reg) 
        with h5py.File(afile, 'r') as f: 
            genes = np.array(f['genes']).astype(str) 
            genes = np.array([g.split('.')[0] for g in genes]) 
            pvals = np.array(f['pearson'])[:,1] 
        reg_genes[reg] = genes 
    
        mask = (pvals <= 0.005) 
        gene_counts[(reg, phn)] = mask.sum() 
        regphen_genes[reg] = np.union1d(regphen_genes[reg], genes[mask])
        regphen_genes[phn] = np.union1d(regphen_genes[phn], genes[mask])

## parse null association data 
null_genes = {i:{k:[] for k in regs+phens} \
              for i in range(nperms)} ## k: perm, v: {reg/phen: unique gene array}
null_genes_top = {i:{k:[] for k in regs+phens} \
                  for i in range(nperms)} ## k: perm, v: {reg/phen: unique gene array}

for phn in phens: 
    for reg in regs: 
        afile = '{}/nulls/pvals_{}/{}.hdf5'.format(path_assoc, phn, reg) 
        with h5py.File(afile, 'r') as f: 
            pvals = np.array(f['null_pearson'])[:,:,1] ## (perms, genes)   
        masks = (pvals <= 0.005)
        for m, mask in enumerate(masks): 
            null_genes[m][reg] = np.union1d(null_genes[m][reg], reg_genes[reg][mask])
            null_genes[m][phn] = np.union1d(null_genes[m][phn], reg_genes[reg][mask])

            topn = gene_counts[(reg,phn)]
            num_nans = np.isnan(pvals[m]).sum()
            a = -(num_nans + topn); b = -num_nans
            idxs = np.argsort(pvals[m])[a:b]

            null_genes_top[m][reg] = np.union1d(null_genes_top[m][reg], reg_genes[reg][idxs])
            null_genes_top[m][phn] = np.union1d(null_genes_top[m][phn], reg_genes[reg][idxs])

## function: compute intersection 
def compute_int(set_keys, set_dict): 
    n_sets = len(set_keys) 
    mat = np.ones((n_sets, n_sets), dtype=int) 
    for i, key1 in enumerate(set_keys): 
        try: mat[i][i] = set_dict[key1].size 
        except: mat[i][i] = 0 
        for ii, key2 in enumerate(set_keys[i+1:]): 
            j = i + ii + 1 
            genes1 = set_dict[key1]
            genes2 = set_dict[key2]
            intersect = np.intersect1d(genes1, genes2).size 
            mat[i][j] = intersect 
            mat[j][i] = intersect
    return mat 

## function: compute normalized jaccard 
def compute_jac(set_keys, set_dict): 
    n_sets = len(set_keys) 
    mat = np.ones((n_sets, n_sets)) 
    for i, key1 in enumerate(set_keys): 
        for ii, key2 in enumerate(set_keys[i+1:]): 
            j = i + ii + 1 
            genes1 = set_dict[key1]
            genes2 = set_dict[key2]

            intersect = np.intersect1d(genes1, genes2).size 
            union = np.union1d(genes1, genes2).size
            if union == 0: union = 1e-6

            jaccard = intersect / union 
            best_jac = min(np.array(genes1).size, np.array(genes2).size) / union 
            if best_jac == 0: best_jac = 1e-6

            norm_jac = jaccard / best_jac 
            mat[i][j] = norm_jac 
            mat[j][i] = norm_jac 
    return mat 

##################################################################################

for grp in ['regs', 'phens']:

    ## gather observed mats 
    int_mat = compute_int(set_type[grp], regphen_genes) 
    jac_mat = compute_jac(set_type[grp], regphen_genes) 

    ## gather null mats 
    if grp == 'regs': n = 10 
    else: n = 5

    int_nulls = np.empty((nperms, n, n))
    jac_nulls = np.empty((nperms, n, n))

    int_nulls_top = np.empty((nperms, n, n))
    jac_nulls_top = np.empty((nperms, n, n))

    for i in range(nperms): 
        int_nulls[i] = compute_int(set_type[grp], null_genes[i])
        jac_nulls[i] = compute_jac(set_type[grp], null_genes[i])

        int_nulls_top[i] = compute_int(set_type[grp], null_genes_top[i])
        jac_nulls_top[i] = compute_jac(set_type[grp], null_genes_top[i])

        if i%100 == 0: print(grp, i)

    ## form int-jac mats 
    mask_int = np.triu(np.ones(n)).astype(bool) ## upper triangle and diagonal 
    mask_jac = np.tril(np.ones(n), k=-1).astype(bool) ## lower triangle 

    obsv_mat = np.where(mask_int, int_mat, 0)
    obsv_mat = np.where(mask_jac, jac_mat, obsv_mat)

    null_mat = np.where(mask_int, int_nulls, 0)
    null_mat = np.where(mask_jac, jac_nulls, null_mat)

    null_mat_top = np.where(mask_int, int_nulls_top, 0)
    null_mat_top = np.where(mask_jac, jac_nulls_top, null_mat_top)

    ## compute pvals 
    pvals = np.mean(null_mat >= obsv_mat, axis=0)
    pvals_top = np.mean(null_mat_top >= obsv_mat, axis=0)

    ## create asterick mats 
    pvals_ast = np.where(pvals <= 0.05, '*', '')
    pvals_ast = np.where(pvals <= 0.01, '**', pvals_ast)

    pvals_ast_top = np.where(pvals_top <= 0.05, '*', '')
    pvals_ast_top = np.where(pvals_top <= 0.01, '**', pvals_ast_top)

    ## save mats 
    with h5py.File('data_inter{}.hdf5'.format(grp[:-1]), 'w') as f: 
        f['obsv_vals'] = obsv_mat 
        f['null_vals'] = null_mat 
        f['pvals'] = pvals 
        f['pvals_topn'] = pvals_top 
        f['pvals_ast'] = pvals_ast.astype(bytes) 
        f['pvals_ast_topn'] = pvals_ast_top.astype(bytes)


