'''
Plot the replication distribution of HCP nonTwins and UKB subsampling. 

- Nhung, Feb 2024
'''

import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

import numpy as np
import pandas as pd

## paths 
main_path = '/Users/nhunghoang/Desktop/remote_platypus/paper_twas'
colr_path = f'{main_path}/aux_files/color_dict.txt'

repl_path = f'{main_path}/outputs_UKB/HCPnonTwin_summary.csv'

plot_path = f'{main_path}/outputs_UKB/plots/nonTwinrepl_boxes.pdf'

## params 
regs = ['dlpfc', 'anterior-cingulate', 'amygdala', 'hippocampus', \
        'caudate', 'putamen', 'nucleus-accumbens', 'cerebellar-hemisphere']
phens = ['vol_mean', 'alff_mean', 'reho_noGS_mean', 'connmean_noGS_mean']

## load data 
df = pd.read_csv(repl_path).groupby('reg')
hcp_data = [] 
ukb_data = [] 
for reg in regs:    
    rdf = df.get_group(reg)
    hcp_data.append(rdf['hcp_nontwin'].values)
    ukb_data.append(rdf['ukb_repl'].values)

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
fig, ax = plt.subplots(1, 1, figsize=(4,2))

xpos = np.arange(1, 17, 2)
fly_props = {'marker': '.', 'markersize': 3, 'markeredgecolor': 'none', 'alpha': 1}

## plot 
bp0 = ax.boxplot(hcp_data, positions=xpos, patch_artist=True, \
                showcaps=False, showmeans=False, flierprops=fly_props)

_ = [box.set_facecolor(rc) for box,rc in zip(bp0['boxes'], colors)]
_ = [box.set_alpha(0.5) for box in bp0['boxes']]
_ = [box.set_edgecolor('none') for box in bp0['boxes']]
_ = [whk.set_color('none') for whk in bp0['whiskers']]
_ = [med.set_color('none') for med in bp0['medians']]

be0 = ax.boxplot(hcp_data, positions=xpos, patch_artist=True, flierprops=fly_props)
_ = [box.set_facecolor('none') for box in be0['boxes']]
_ = [med.set_color('k') for med in be0['medians']]

bp1 = ax.boxplot(ukb_data, positions=xpos+0.75, patch_artist=True, flierprops=fly_props)
_ = [box.set_facecolor('none') for box in bp1['boxes']]
_ = [box.set_edgecolor(rc) for box,rc in zip(bp1['boxes'], colors)]
_ = [med.set_color(rc) for med,rc in zip(bp1['medians'], colors)]

cc = [c for c in colors for _ in range(2)]
_ = [whk.set_color(c) for whk,c in zip(bp1['whiskers'], cc)]
_ = [cap.set_color(c) for cap,c in zip(bp1['caps'], cc)]

## format 
ax.set_xticks([i + 0.25 for i in xpos])
ax.set_xticklabels([reg[:5] for reg in regs], size=8)
ax.set_ylabel('$r^2$', size=10)
ax.set_title('fraction of replication', size=8)

## legend 
hcp = mpatches.Patch(facecolor='#bbbbbb', edgecolor='k', label='HCP nonTwin')  
ukb = mpatches.Patch(facecolor='none', edgecolor='k', label='UKB sample')  

patches = [hcp, ukb]
ax.legend(handles=patches, bbox_to_anchor=(0.001, 1), loc='upper left', ncol=2, fontsize=7)

## save 
plt.tight_layout()
plt.savefig(plot_path, format='pdf')
    

