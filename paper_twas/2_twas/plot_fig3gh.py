'''
Figure 3G-H: Heatmaps for cross-region TWAS results 

- Nhung, Aug 2023 
'''

import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt
import seaborn as sns

import numpy as np 
import pandas as pd 
import sys 

from mpl_toolkits.axes_grid1 import make_axes_locatable
from statsmodels.stats.multitest import multipletests as sm
from sklearn.preprocessing import minmax_scale as mm_scale
from sctriangulate.colors import build_custom_continuous_cmap as BCCC

## args 
group = 'UKB' 
phens = 'vol_mean'  
ptype = 'bonf' ## 'bonf'

## paths
main_path = f'/Users/nhunghoang/Desktop/remote_platypus/paper_twas/outputs_{group}'

twas_path = f'{main_path}/twas_JTI/{phens}'
itwa_path = f'{main_path}/twas_JTI/cross_regs/{phens}'

plot_path = f'{main_path}/plots'
plot_file = f'{plot_path}/fig3gh.pdf'
if ptype == 'bonf': 
    plot_file = f'{plot_path}/fig_s2cd.pdf'

## regions
regs = ['dlpfc', 'anterior-cingulate', 'amygdala', 'hippocampus', \
        'caudate', 'putamen', 'nucleus-accumbens', 'cerebellar-hemisphere']

name = ['DLPFC', 'ant. cingulate', 'amygdala', 'hippocampus', \
        'caudate', 'putamen', 'nuc. accumbens', 'cerebellum']

################################################################################

## load twas results 
nr = len(regs)
num_mod = np.zeros((nr, nr), dtype=int)
num_sig = np.zeros((nr, nr), dtype=int)
sig_gen = {} ## k: (grex reg, phen reg), v: array of sig genes 

for ig, grex_reg in enumerate(regs): 
    for ip, phen_reg in enumerate(regs): 

        if grex_reg == phen_reg: 
            tfile = f'{twas_path}/{grex_reg}.txt' 
        else: 
            tfile = f'{itwa_path}/grex_{grex_reg}_phen_{phen_reg}.txt'

        df = pd.read_table(tfile, sep='\t', usecols=['gene', 'pvalue'])

        ## FDR/BON correction 
        if ptype == 'bonf':
            keep, fval, _, _ = sm(df['pvalue'], method='bonferroni', alpha=0.05)
        else:
            keep, fval, _, _ = sm(df['pvalue'], method='fdr_bh', alpha=0.05)

        ## store sig results  
        num_mod[ip][ig] = keep.size
        num_sig[ip][ig] = keep.sum()
        sig_gen[(grex_reg, phen_reg)] = df.loc[keep]['gene'].values 

################################################################################

## func: compute normalized Jaccard
def jac_norm(g1, g2):
    ic = np.intersect1d(g1, g2).size
    un = np.union1d(g1, g2).size
    jc = ic / un

    jbest = np.min([g1.size, g2.size]) / un
    jnorm = jc / jbest

    return jnorm

## compute Jaccard (per row, wrt to diagonal) 
jac_mat = np.zeros((nr, nr), dtype=float)
for ig, reg in enumerate(regs): 

    ## sig genes for this reg
    diag_arr = sig_gen[(reg, reg)] 

    for ip, phen_reg in enumerate(regs): 

        ## sig genes for this phenotype (wrt to grex) 
        phen_arr = sig_gen[(reg, phen_reg)]

        ## jac 
        jac_mat[ip][ig] = jac_norm(diag_arr, phen_arr)

################################################################################

## init plot 
plt.ion()
fig, axes = plt.subplots(1,2, figsize=(8,3))
kwargs = {'linecolor': '#708090', 'linewidths': 0.01, 'square': True, 'clip_on': True, \
          'cbar_kws': {'shrink': 0.4, 'orientation': 'horizontal'}, 'annot_kws': {'size': 8}} 

pink_cmap = BCCC([255, 247, 247], [254, 180, 198], [255, 52, 104], [149, 53, 78])

## make space for colorbars
divide0 = make_axes_locatable(axes[0])
div_num = divide0.append_axes('top', size='3%', pad=0.05)

divide1 = make_axes_locatable(axes[1])
div_jac = divide1.append_axes('top', size='3%', pad=0.05)

################################################################################

## heatmap: number of significant genes 
_ = sns.heatmap(num_sig, ax=axes[0], annot=num_sig, fmt='g', \
                cbar_ax=div_num, cmap=pink_cmap, **kwargs)

## heatmap: Jaccard (per row, wrt to diagonal) 
def format_float(val): return ('{:.2f}'.format(val).lstrip('0') or '0')
jm = np.vectorize(format_float)(jac_mat)
for i in range(len(regs)): 
    jm[i][i] = '1.0'

_ = sns.heatmap(jac_mat, ax=axes[1], annot=jm, fmt='s', vmin=0, \
                cbar_ax=div_jac, cmap='Oranges', **kwargs)

## formatting 
for ax in axes: 
    ax.set_xlabel('GREx', size=7) 
    ax.set_ylabel('mean volume', size=7)
    ax.set_xticklabels(name, size=7, rotation=90, ha='right')
    ax.set_yticklabels(name, size=7, rotation=0)
axes[1].set_yticks([])
axes[1].set_yticklabels('')

ca0 = axes[0].collections[0].colorbar.ax
ca1 = axes[1].collections[0].colorbar.ax

if ptype == 'bonf': rng = np.arange(15, 51, 5)
else: rng = np.arange(80, 201, 40) 

ca0.set_xticks(rng)
ca1.set_xticks(np.arange(0, 1.1, 0.2))

tps = {'labelsize': 7, 'top': True, 'bottom': False, \
       'labeltop': True, 'labelbottom': False}
ca0.tick_params(**tps)
ca1.tick_params(**tps)

## save 
plt.tight_layout()
plt.savefig(plot_file, format='pdf')



