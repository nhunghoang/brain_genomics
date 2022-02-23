'''
Inter-regional gene overlaps: 
- all PrediXcan models 
- PrediXcan models after quality filter 
- singularly significant genes 
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
import bct 

regs = os.listdir('/data1/rubinov_lab/brain_genomics/analyses_HCP/lofo_assoc/train_pvals_alff')
regs = np.array(sorted(regs)) ## abc order 

########################################################################################

## function: count inter-regional gene overlaps 
## and compute Jaccard similarity  
def count_genes(genes): 
    count_int = np.zeros((10,10), dtype=int) ## intersection matrix 
    count_jac = np.zeros((10,10), dtype=float) ## jaccard matrix 

    for i, reg1 in enumerate(regs):
        count_int[i][i] = genes[reg1].size 
        count_jac[i][i] = 1 

        for ii, reg2 in enumerate(regs[i+1:]):
            j = i + ii + 1
            genes1 = genes[reg1]
            genes2 = genes[reg2] 
        
            ## compute gene intersection and normalized jaccard
            intersect = np.intersect1d(genes1, genes2).size
            union = np.union1d(genes1, genes2).size
            jaccard = intersect / union 
            best_jac = min(genes1.size, genes2.size) / union 
            norm_jac = jaccard / best_jac 

            count_int[i][j] = intersect; count_int[j][i] = intersect 
            count_jac[i][j] = norm_jac; count_jac[j][i] = norm_jac
    return {'count_int':count_int, 'count_jac':count_jac}

## function: plot gene counts as heatmap 
def plot_gene_counts(counts, out_path):
    count_int = counts['count_int']
    count_jac = counts['count_jac']

    ## clustering by similarity (Louvain)
    #for g in np.arange(1,10):
    #    clusters_jac, Q_jac = bct.community_louvain(count_jac, gamma=g)
    #    print('gamma: {:.2f} ({} clusters)'.format(g, np.unique(clusters_jac).size))
    #    if np.unique(clusters_jac).size >= 3: break

    ## reorder matrix to visually group clusters 
    #_, idx_jac = bct.grid_communities(clusters_jac) 

    idx = np.array([5,0,6,9,2,8,7,1,4,3])
    reordered_int = count_int[idx][:,idx]
    reordered_jac = count_jac[idx][:,idx]

    ## plot
    off_diag = (np.ones((10,10)) - np.eye(10)).astype(bool)
    vmax_int = (int(reordered_int[off_diag].max() / 10) + 1) * 10 
    #vmax_jac = ((reordered_jac[off_diag].max() * 10).astype(int) + 1) / 10
    vmax_jac = 0.3
    vmin = 0.1

    fig, axes = plt.subplots(1, 2, figsize=(60,20))
    kwargs_int = {'vmax':vmax_int, 'fmt':'d', 'cmap':'Greens', 'annot_kws':{'size':35}, 'linewidth':1, 'linecolor':'w'}
    kwargs_jac = {'vmax':vmax_jac, 'fmt':'.2f', 'cmap':'Oranges', 'annot_kws':{'size':35}, 'linewidth':1, 'linecolor':'w'}

    ## plot full matrix that clustering is based on 
    sns.heatmap(reordered_jac, annot=True, vmin=vmin, **kwargs_jac, ax=axes[0]) 
    axes[0].collections[0].colorbar.ax.set_ylabel('normalized Jaccard similarity', size=40)
    axes[0].collections[0].colorbar.ax.tick_params(labelsize=30)
    axes[0].collections[0].colorbar.ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

    ## plot intersection and jaccard off diagonals   
    int_mask = np.triu(np.ones((10,10)), k=1).astype(bool) ## lower triangle and diagonal  
    jac_mask = np.tril(np.ones((10,10)), k=0).astype(bool) ## upper triangle 

    sns.heatmap(reordered_int, mask=int_mask, vmin=0, annot=True, **kwargs_int) 
    sns.heatmap(reordered_jac, mask=jac_mask, vmin=vmin, annot=True, **kwargs_jac, cbar=False) 
    axes[1].collections[0].colorbar.ax.set_ylabel('intersect count', size=40)
    axes[1].collections[0].colorbar.ax.tick_params(labelsize=30)

    for ax in axes:
        ax.set_xticks(np.arange(10) + 0.5)
        ax.set_yticks(np.arange(10) + 0.5)
        ax.set_xticklabels(regs[idx], rotation=20, ha='right', size=30)
        ax.set_yticklabels(regs[idx], rotation=0, size=30) 
        xleft, xright = ax.get_xlim()
        ybottom, ytop = ax.get_ylim()
        ax.set_aspect(abs((xright-xleft)/(ybottom-ytop)))
    plt.tight_layout()
    plt.savefig(out_path)
    plt.close('all')     

########################################################################################

def main(): 

    ## all PrediXcan genes 
    #expr_path = '/data1/rubinov_lab/brain_genomics/data_HCP/expression'
    #out_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/interReg_overlap/num_common_genes-predixcan.png'
    #genes = {} ## k: reg, v: gene list 
    #for reg in regs:
    #    with h5py.File('{}/{}.hdf5'.format(expr_path, reg), 'r') as f:
    #        genes[reg] = np.array(f['genes']) 
    #counts = count_genes(genes) 
    #plot_gene_counts(counts, out_path) 

    ## quality-filtered PrediXcan genes
    #expr_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/assoc/pvals_alff'
    #out_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/interReg_overlap/num_common_genes-filtered-jac-ordered.png'
    #genes = {} ## k: reg, v: gene list 
    #for reg in regs:
    #    with open('{}/{}.txt'.format(expr_path, reg), 'r') as f: 
    #        f.readline()
    #        glist = [i.split('\t')[0].split('.')[0] for i in f.readlines()]
    #        genes[reg] = np.array(glist)
    #counts = count_genes(genes) 
    #plot_gene_counts(counts, out_path) 

    ## independently selected genes 
    #expr_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/interReg_overlap/single_gene_unions'
    #out_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/interReg_overlap/num_common_genes-single-jac-ordered.png'
    #genes = {} ## k: reg, v: gene list 
    #for reg in regs:
    #   with h5py.File('{}/{}.hdf5'.format(expr_path, reg), 'r') as f:
    #       genes[reg] = np.array(f['genes']) 
    #counts = count_genes(genes) 
    #plot_gene_counts(counts, out_path) 

    ## selection in gene sets 
    #expr_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/DATA_OUTPUT/CCC_expr_regress'
    #out_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/interReg_overlap/num_common_genes-multi-jac-ordered.png'
    #genes = {} ## k: reg, v: gene list 
    #for reg in regs:
    #   with h5py.File('{}/{}.hdf5'.format(expr_path, reg), 'r') as f:
    #       genes[reg] = np.array(f['genes']) 
    #counts = count_genes(genes) 
    #plot_gene_counts(counts, out_path) 

    ## all PrediXcan genes - single assoc
    expr_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/interReg_overlap/single_gene_unions_FDR'
    out_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/interReg_overlap/num_common_genes-single-FDR.png'
    genes = {} ## k: reg, v: gene list 
    for reg in regs:
       with h5py.File('{}/{}.hdf5'.format(expr_path, reg), 'r') as f:
           genes[reg] = np.array(f['genes']) 
    counts = count_genes(genes) 
    plot_gene_counts(counts, out_path) 

    ## all PrediXcan genes - multi assoc
    #expr_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/DATA_OUTPUT/CCC_expr_regress-allgenes'
    #out_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/interReg_overlap/num_common_genes-multi-allgenes.png'
    #genes = {} ## k: reg, v: gene list 
    #for reg in regs:
    #   with h5py.File('{}/{}.hdf5'.format(expr_path, reg), 'r') as f:
    #       genes[reg] = np.array(f['genes']) 
    #counts = count_genes(genes) 
    #plot_gene_counts(counts, out_path) 

main()
