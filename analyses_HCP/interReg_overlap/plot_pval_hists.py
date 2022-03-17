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
import matplotlib.patches as mpatches
import pandas as pd 

phens = ['alff', 'regional_homogeneity', 'gm_volume']
phens_short = ['ALFF', 'ReHo', 'GMVol']
phens_color = ['blue', 'orange', 'green'] 

reg_idx = np.array([5,0,6,9,2,8,7,1,4,3])
regs = os.listdir('/data1/rubinov_lab/brain_genomics/analyses_HCP/assoc/null_pvals_alff')
regs = np.array(sorted(regs)) ## abc order
regs = regs[reg_idx]

pval_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/assoc/pvals'  
pval_col = 3 ## 2 for uncorrected, 3 for FDR-corrected 

## create figure 
fig, axes = plt.subplots(2, 5, figsize=(60,20))

## gather p-values and plot histograms 
for r,reg in enumerate(regs): 
    y = int(r/5)
    x = int(r%5) 
    ax = axes[y][x]

    ## gather p-values 
    phen_data = [] 
    for p,phen in enumerate(phens):
        pval_file = '{}_{}/{}.txt'.format(pval_path, phen, reg)
        pval_data = np.loadtxt(pval_file, delimiter='\t', skiprows=1, usecols=[pval_col]) 
        pval_data = pd.DataFrame(pval_data, columns=['p-value']).assign(phenotype=phens_short[p])
        phen_data.append(pval_data) 

    ## plot p-values on histogram 
    hist_data = pd.concat(phen_data)
    sns.histplot(hist_data, x='p-value', hue='phenotype', \
        element='step', ax=ax, zorder=3, legend=True, \
        fill=False, linewidth=3, bins=20)#, binrange=[0,0.05])

    ax.tick_params(labelsize=18)
    ax.set_xlabel('p-value', size=18, weight='bold')
    ax.set_ylabel('count', size=18, weight='bold')
    ax.set_title(reg, size=24, weight='bold')

    ax.grid(color='gray', linestyle='-', alpha=0.4, zorder=1)
    
#phen_patches = [] 
#for p in range(3): 
#    patch = mpatches.Patch(facecolor=phens_color[p], alpha=0.5, \
#            edgecolor=phens_color[p], linewidth=2, \
#            label=phens_short[p])
#fig.legend(handles=phen_patches, loc='lower center')

plt.savefig('phen_pval_hist_FDR.png') 
plt.close('all') 
sys.exit()

