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
path_main = '/data1/rubinov_lab/brain_genomics/scripts_twas'
path_assoc = path_main + '/outputs_HCP/assoc_1M'
path_out = path_main + '/outputs_HCP/regphen_similarities'

## variables 
nperms = 10000
phens = ['gm_volume', 'myelination', 'alff', 'reho_noGS', 'connmean_noGS']
regs = ['frontal-pole', 'anterior-cingulate', 'caudate', 'putamen', 'nucleus-accumbens', \
        'hippocampus', 'amygdala', 'hypothalamus', 'substantia-nigra', 'cerebellar-hemisphere']
regphens = [(reg,phen) for reg in regs for phen in phens]
set_type = {'phens': phens, 'regs': regs} 

## parse association results 
regphen_genes = {k:[] for k in regs+phens} ## k: reg or phen, v: unique gene array 
gene_counts = {} ## k: (reg, phen), v: number of sig genes 
reg_genes = {} ## k: reg, v: gene array 

for (reg, phn) in regphens:
    afile = '{}/pvals_{}/{}.hdf5'.format(path_assoc, phn, reg) 
    with h5py.File(afile, 'r') as f: 
        genes = np.array(f['genes']).astype(str) 
        pvals = np.array(f['pearson'])[:,1] 
    reg_genes[reg] = genes 

    mask = (pvals <= 0.005) 
    gene_counts[(reg, phn)] = mask.sum() 
    regphen_genes[reg] = np.union1d(regphen_genes[reg], genes[mask])
    regphen_genes[phn] = np.union1d(regphen_genes[phn], genes[mask])

## parse null association results 
null_genes_sig = {i:{k:[] for k in regs+phens} \
                  for i in range(nperms)} ## k: perm, v: {reg/phen: unique gene array}
null_genes_top = {i:{k:[] for k in regs+phens} \
                  for i in range(nperms)} ## k: perm, v: {reg/phen: unique gene array}

for (reg, phn) in regphens:
    afile = '{}/nulls/pvals_{}/{}.hdf5'.format(path_assoc, phn, reg) 
    with h5py.File(afile, 'r') as f: 
        pvals = np.array(f['null_pearson'])[:,:,1] ## (perms, genes)   

    masks = (pvals <= 0.005)
    for m, mask in enumerate(masks): 

        ## permutation genes based on (assoc p<= 0.005)
        null_genes_sig[m][reg] = np.union1d(null_genes[m][reg], reg_genes[reg][mask])
        null_genes_sig[m][phn] = np.union1d(null_genes[m][phn], reg_genes[reg][mask])

        ## top K permutation genes (ordered by signif), where 
        ## K is the num of signif genes in observed assocs 
        topn = gene_counts[(reg,phn)]
        num_nans = np.isnan(pvals[m]).sum()
        a = -(num_nans + topn); b = -num_nans
        idxs = np.argsort(pvals[m])[a:b]

        null_genes_top[m][reg] = np.union1d(null_genes_top[m][reg], reg_genes[reg][idxs])
        null_genes_top[m][phn] = np.union1d(null_genes_top[m][phn], reg_genes[reg][idxs])

##################################################################################

## function: count intersection size between gene sets  
## return the set of common genes too, if counting observed data 
def compute_int(set_keys, set_dict, obsv=False): 
    n_sets = len(set_keys) 
    mat = np.ones((n_sets, n_sets), dtype=int) 
    iset = {} ## k: (key1, key2), v: common genes 
    for i, key1 in enumerate(set_keys): 
        try: mat[i][i] = set_dict[key1].size 
        except: mat[i][i] = 0 
        for ii, key2 in enumerate(set_keys[i+1:]): 
            j = i + ii + 1 
            genes1 = set_dict[key1]
            genes2 = set_dict[key2]
            intersect = np.intersect1d(genes1, genes2) 
            mat[i][j] = intersect.size 
            mat[j][i] = intersect.size

            if obsv: iset[(key1, key2)] = intersect 

    if obsv: return mat, iset
    else: return mat

## function: compute normalized jaccard for gene sets  
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
    int_mat, int_sets = compute_int(set_type[grp], regphen_genes, obsv=True) 
    jac_mat = compute_jac(set_type[grp], regphen_genes) 

    ## gather null mats 
    if grp == 'regs': n = 10 
    else: n = 5

    int_nulls_sig = np.empty((nperms, n, n))
    jac_nulls_sig = np.empty((nperms, n, n))

    int_nulls_top = np.empty((nperms, n, n))
    jac_nulls_top = np.empty((nperms, n, n))

    for i in range(nperms): 
        int_nulls_sig[i] = compute_int(set_type[grp], null_genes_sig[i])
        jac_nulls_sig[i] = compute_jac(set_type[grp], null_genes_sig[i])

        int_nulls_top[i] = compute_int(set_type[grp], null_genes_top[i])
        jac_nulls_top[i] = compute_jac(set_type[grp], null_genes_top[i])

        if i%100 == 0: print(grp, i)

    ## form int-jac mats 
    mask_int = np.triu(np.ones(n)).astype(bool) ## upper triangle and diagonal 
    mask_jac = np.tril(np.ones(n), k=-1).astype(bool) ## lower triangle 

    obsv_mat = np.where(mask_int, int_mat, 0)
    obsv_mat = np.where(mask_jac, jac_mat, obsv_mat)

    null_mat_sig = np.where(mask_int, int_nulls, 0)
    null_mat_sig = np.where(mask_jac, jac_nulls_sig, null_mat_sig)

    null_mat_top = np.where(mask_int, int_nulls_top, 0)
    null_mat_top = np.where(mask_jac, jac_nulls_top, null_mat_top)

    ## compute pvals 
    pvals_sig = np.mean(null_mat_sig >= obsv_mat, axis=0)
    pvals_top = np.mean(null_mat_top >= obsv_mat, axis=0)

    ## save results  
    ofile = '{}/mats_inter{}.hdf5'.format(path_out, grp[:-1]) 
    with h5py.File(ofile, 'w') as f: 
        f['obsv_vals'] = obsv_mat 
        f['null_vals_sig'] = null_mat_sig 
        f['null_vals_top'] = null_mat_tosigp 

        f['pvals_sig'] = pvals_sig 
        f['pvals_top'] = pvals_top 

        ## save gene sets  
        for (k1, k2), iset in int_sets.items():
            k = 'genes_{}_{}'.format(k1, k2)
            f[k] = iset.astype(bytes)


