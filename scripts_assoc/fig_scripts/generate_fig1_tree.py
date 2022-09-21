'''
Figure 1E: Dendrogram of gr-expression for HCP cohort  

Nhung, updated Aug 2022
'''

import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.cluster.hierarchy import linkage
import sys
import h5py
import seaborn as sns 

expr_dir = '/Users/nhunghoang/Desktop/local_platypus/expr_noregr' 
dem_file = '/Users/nhunghoang/Desktop/local_platypus/cohort890_demographics.txt'
cov_file = '/Users/nhunghoang/Desktop/local_platypus/covariates.txt'

## parse demographics  
with open(dem_file, 'r') as f: 

    ## subjID sampID famID gender age race ethnicity 
    info = np.loadtxt(dem_file, delimiter='\t', skiprows=1, usecols=[2,5], dtype=str)
    fam = info[:,0]
    pop = info[:,1]

pop_color = {'White':'r', 
             'Black or African Am.':'g',
             'Asian/Nat. Hawaiian/Othr Pacific Is.':'b',
             'Unknown or Not Reported': 'k',   
             'Am. Indian/Alaskan Nat.': 'm', 
             'More than one': 'y'}

## gather expression
all_expr = [] 
for regf in os.listdir(expr_dir): 
    with h5py.File('{}/{}'.format(expr_dir, regf), 'r') as f: 
        reg_expr = np.array(f['pred_expr']).T
        all_expr.append(reg_expr)
cat_expr = np.concatenate(all_expr, axis=1) ## (subjs, genes)

## compute linkage 
link_mat = linkage(cat_expr, method='ward')

## generate clustermap 
plt.ion()

pop_row_colors = np.array([pop_color[p] for p in pop])
cm = sns.clustermap(cat_expr, row_linkage=link_mat, \
                    row_cluster=True, col_cluster=False, \
                    row_colors=pop_row_colors, yticklabels=1, \
                    tree_kws={'colors':'k'})

## remove heatmap and colorbar 
plt.delaxes(cm.ax_heatmap)
plt.delaxes(cm.cax)

## save dendrogram (pdf takes too long to generate) 
plt.savefig('fig_plots/fig1_tree.png', dpi=400) 

## generate dendrogram inset 
dax = cm.ax_row_dendrogram
rax = cm.ax_row_colors

## arrange family IDs based on tree indices 
fdict = {f:i for i,f in enumerate(np.unique(fam))}
idxs = cm.dendrogram_row.reordered_ind
fnum = [fdict[f] for f in fam[idxs]]

## set row and label colors based on ordering 
pop_clr = pop_row_colors[idxs]
fam_lab = ['fam {}'.format(f) for f in fnum] 

## set family labels 
rax.set_yticks(np.arange(890)+0.5)
rax.set_yticklabels(fam_lab, fontsize=14)
for tlabel, tcolor in zip(rax.get_yticklabels(), pop_clr):
    tlabel.set_color(tcolor)

## zoom into location of interest
dax.set_xlim([70, 0])

d1 = 1550; d2 = 1350
ratio = rax.get_ylim()[0] / dax.get_ylim()[0]
dax.set_xlim([70, 0])
dax.set_ylim([d1, d2])
rax.set_ylim([d1*ratio, d2*ratio])
plt.tight_layout()

## save inset 
plt.savefig('fig_plots/fig1_tree-inset.png', dpi=400)

