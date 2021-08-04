'''
For a given regional phenotype, compute corr(phen, expr) for every gene. 
Then, FDR-correct the p-values of these gene-specific correlations, and 
plot some of them (expr on x, phen on y, subjs are points). 

- Nhung, June 2021  
'''

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import os
import numpy as np
import h5py
from scipy.stats import pearsonr, spearmanr
from statsmodels.stats.multitest import fdrcorrection
import read_data as rd 

region = sys.argv[1] 
phenotype = sys.argv[2] 

## output path  
cpath = '/data1/rubinov_lab/brain_genomics/analyses_HCP_v7/plots_single_genes'
if not os.path.exists(cpath): os.mkdir(cpath) 

## align sample & subject IDs 
sfile = '/data1/rubinov_lab/brain_genomics_accre/scripts/coexpr_coactivity/sim_annealing/sa_variables.hdf5'
with h5py.File(sfile, 'r') as f:
    nid_order = np.array(f['nid_order'])
    gid_order = np.array(f['gid_order'])

## fetch regional phenotype 
regions = ['amygdala', 'hippocampus', 'putamenBG', 'aCingulateCortex', \
            'frontalCortex', 'nAccumbensBG', 'caudateBG', 'cerebellum']
r_idx = regions.index(region) 
pfile = '/data1/rubinov_lab/brain_genomics_accre/data_HCP/phenotypes2.hdf5'
with h5py.File(pfile, 'r') as f:
    phen_mat = np.array(f[phenotype])
phen = phen_mat[1:-1][:,r_idx] 
if phenotype == 'avgconn': phen /= 228 ## (115 regions * 2 - self) 

## fetch regional expression and gene names   
expr_gids, expr_dict = rd.HCP_v7_expr([region])
expr = expr_dict[region]['exprs'][gid_order][1:-1].T ## (genes * subjs)  
genes = expr_dict[region]['genes']

## compute gene-to-phen correlations  
data = np.zeros((genes.shape[0], 3))
for g,gene in enumerate(genes):
    gene_expr = expr[g]

    ## return nan for genes with no expression variance
    if np.var(gene_expr) == 0:
        data[g,:] = np.nan
        continue

    ## otherwise, compute correlation
    rho, pval = pearsonr(phen, gene_expr)
    data[g,0] = rho; data[g,1] = pval
genes_psig = genes[data[:,1] <= 0.05] 
print('{} p sig'.format(genes_psig.shape[0]))

## compute FDR-correct pvals
pvals = data[:,1]
pvals = pvals[~np.isnan(pvals)]
rejected, corrected = fdrcorrection(pvals)
c = 0
for g in range(genes.shape[0]):
    if np.isnan(data[g,2]): continue
    else: data[g,2] = corrected[c]; c += 1
genes_fsig = genes[data[:,2] <= 0.05] 
print('{} f sig'.format(genes_fsig.shape[0]))

## scatter-plot some genes 
n_plot = 0  
for gf in genes_fsig:
    if n_plot > 5: break
    n_plot += 1
    gf_idx = np.argwhere(genes == gf)[0][0]
    gf_expr = expr[gf_idx]
    rho, pval = pearsonr(gf_expr, phen)
    b, m = np.polynomial.polynomial.polyfit(gf_expr, phen, 1)

    fig, ax = plt.subplots(1,1,figsize=(15,15))
    ax.scatter(gf_expr, phen, c='k')
    ax.plot(gf_expr, b+(m*gf_expr), '-', c='y')
    ax.set_xlabel('\nexpression', fontsize=30)
    ax.set_ylabel(phenotype + '\n', fontsize=30)
    ax.tick_params(labelsize=40)

    xleft, xright = ax.get_xlim()
    ybottom, ytop = ax.get_ylim()
    ax.set_aspect(abs((xright-xleft)/(ybottom-ytop)))

    title = '{}, {}\nr={:.3f} (p ≤ {:.3f})\n'.format(gf, region, rho, pval)
    ax.set_title(title, size=30)
    fname = '{}/{}_{}_{}_{:.3f}.png'.format(cpath, phenotype, region, gf, rho)
    plt.savefig(fname)
    plt.close('all')

