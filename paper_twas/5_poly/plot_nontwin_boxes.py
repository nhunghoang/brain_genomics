'''
In line with Fig 5C, plot the 'top k gene' 
boxes for the nonTwin distribution. 

- Nhung, Feb 2024
'''

import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

import numpy as np
import pandas as pd
import h5py

## paths 
main_path = '/Users/nhunghoang/Desktop/remote_platypus/paper_twas'
colr_path = f'{main_path}/aux_files/color_dict.txt'

euro_path = f'{main_path}/outputs_HCP/allEuro/polygenic_models/split_r2s.hdf5'
ntwn_path = f'{main_path}/outputs_HCP/nonTwin/polygenic_models/split_r2s.hdf5'

plot_path = f'{main_path}/outputs_HCP/plots/nonTwin_boxes.pdf'

## params 
regs = ['dlpfc', 'anterior-cingulate', 'amygdala', 'hippocampus', \
        'caudate', 'putamen', 'nucleus-accumbens', 'cerebellar-hemisphere']
phens = ['vol_mean', 'alff_mean', 'reho_noGS_mean', 'connmean_noGS_mean']

## load data 
eu_r2s = {} ## k: (reg, phen), v: (n_split,) 
nt_r2s = {} ## k: (reg, phen), v: (n_split,) 

ef = h5py.File(euro_path, 'r') 
nf = h5py.File(ntwn_path, 'r') 
for key in ef.keys(): 
    [reg, phen] = key.split('X')
    eu_r2s[(reg, phen)] = ef[key][()]
    nt_r2s[(reg, phen)] = nf[key][()]
ef.close(); nf.close()

## load colors 
colors = pd.read_table(colr_path, index_col='label').to_dict()['hex']
colors = [colors[r] for r in regs]

## func: customize boxplot
def customize(bp, which):
    if which == 'face':
        _ = [box.set_facecolor(rc) for box,rc in zip(bp['boxes'], colors)]
        _ = [box.set_alpha(0.5) for box in bp['boxes']]
        _ = [box.set_edgecolor('none') for box in bp['boxes']]
        _ = [whk.set_color('none') for whk in bp['whiskers']]
        _ = [med.set_color('none') for med in bp['medians']]
    if which == 'edge':
        _ = [box.set_facecolor('none') for box in bp['boxes']]
        _ = [med.set_color('k') for med in bp['medians']]

## plot params 
plt.ion() 
fig, axes_ = plt.subplots(4, 1, figsize=(4,8))
axes = {p: ax for p,ax in zip(phens, axes_.flatten())}

xpos = np.arange(1, 17, 2)
fly_props = {'marker': '.', 'markersize': 3, 'markeredgecolor': 'none', 'alpha': 1}

## plot 
for phen, ax in axes.items():
    eu_data = [eu_r2s[(reg, phen)] for reg in regs]
    nt_data = [nt_r2s[(reg, phen)] for reg in regs]

    ## euro
    ebp = ax.boxplot(eu_data, positions=xpos, patch_artist=True, \
                    showcaps=False, showmeans=False, flierprops=fly_props)

    _ = [box.set_facecolor(rc) for box,rc in zip(ebp['boxes'], colors)]
    _ = [box.set_alpha(0.5) for box in ebp['boxes']]
    _ = [box.set_edgecolor('none') for box in ebp['boxes']]
    _ = [whk.set_color('none') for whk in ebp['whiskers']]
    _ = [med.set_color('none') for med in ebp['medians']]

    ebe = ax.boxplot(eu_data, positions=xpos, patch_artist=True, flierprops=fly_props)
    _ = [box.set_facecolor('none') for box in ebe['boxes']]
    _ = [med.set_color('k') for med in ebe['medians']]

    ## nontwin
    nbp = ax.boxplot(nt_data, positions=xpos+0.75, patch_artist=True, flierprops=fly_props)
    _ = [box.set_facecolor('none') for box in nbp['boxes']]
    _ = [box.set_edgecolor(rc) for box,rc in zip(nbp['boxes'], colors)]
    _ = [med.set_color(rc) for med,rc in zip(nbp['medians'], colors)]

    cc = [c for c in colors for _ in range(2)]
    _ = [whk.set_color(c) for whk,c in zip(nbp['whiskers'], cc)]
    _ = [cap.set_color(c) for cap,c in zip(nbp['caps'], cc)]

    ## format 
    ax.set_xticks([i + 0.25 for i in xpos])
    ax.set_xticklabels([reg[:5] for reg in regs], size=8)
    ax.set_ylabel('$r^2$', size=10)
    ax.set_title(phen, size=8)

## legend 
ax = axes[phens[0]]  
eu = mpatches.Patch(facecolor='#bbbbbb', edgecolor='k', label='allEuro')  
nt = mpatches.Patch(facecolor='none', edgecolor='k', label='nonTwin')  

patches = [eu, nt]
ax.legend(handles=patches, bbox_to_anchor=(0.001, 1), loc='upper left', ncol=2, fontsize=7)

## save 
plt.tight_layout()
plt.savefig(plot_path, format='pdf')
    

