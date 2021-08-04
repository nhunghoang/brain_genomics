'''
Generate two scatter plots for a given regional phenotype: 
- phenotype on y-axis; points are subjects 
- on x-axis: sum(all_genes); sum(selected_genes) based on sim_annealing runs 

Nhung, June 2021 
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
cpath = '/data1/rubinov_lab/brain_genomics/analyses_HCP_v7/plots_multi_genes'
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

## scatter-plot all genes against phenotype 
expr_all = np.sum(expr, axis=0)
rho, pval = pearsonr(expr_all, phen)
b, m = np.polynomial.polynomial.polyfit(expr_all, phen, 1)

fig, ax = plt.subplots(1,1,figsize=(15,15))
ax.scatter(expr_all, phen, c='k')
ax.plot(expr_all, b+(m*expr_all), '-', c='y')
ax.set_xlabel('\nexpression (all genes)', fontsize=30)
ax.set_ylabel(phenotype + '\n', fontsize=30)
ax.tick_params(labelsize=40)

xleft, xright = ax.get_xlim()
ybottom, ytop = ax.get_ylim()
ax.set_aspect(abs((xright-xleft)/(ybottom-ytop)))

title = 'All Genes, {}\nr={:.3f} (p ≤ {:.3f})\n'.format(region, rho, pval)
ax.set_title(title, size=30)
fname = '{}/{}_{}_allgenes_{:.3f}.png'.format(cpath, phenotype, region, rho)
plt.savefig(fname)
plt.close('all')

## scatter-plot selected genes against phenotype 
gfile = '/data1/rubinov_lab/brain_genomics_accre/scripts/coexpr_coactivity/sim_annealing/2021_{}/bestEmprCorrs.hdf5'.format(phenotype) 
with h5py.File(gfile, 'r') as f: 
    expr_sel = np.array(f['{}-128'.format(region)])[0] 

rho, pval = pearsonr(expr_sel, phen)
b, m = np.polynomial.polynomial.polyfit(expr_sel, phen, 1)

fig, ax = plt.subplots(1,1,figsize=(15,15))
ax.scatter(expr_sel, phen, c='k')
ax.plot(expr_sel, b+(m*expr_sel), '-', c='y')
ax.set_xlabel('\nexpression (128 genes)', fontsize=30)
ax.set_ylabel(phenotype + '\n', fontsize=30)
ax.tick_params(labelsize=40)

xleft, xright = ax.get_xlim()
ybottom, ytop = ax.get_ylim()
ax.set_aspect(abs((xright-xleft)/(ybottom-ytop)))

title = 'Selected Genes, {}\nr={:.3f} (p ≤ {:.3f})\n'.format(region, rho, pval)
ax.set_title(title, size=30)
fname = '{}/{}_{}_128genes_{:.3f}.png'.format(cpath, phenotype, region, rho)
plt.savefig(fname)
plt.close('all')
