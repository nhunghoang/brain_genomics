'''
Plot box plots of the replications. 

- Nhung, Feb 2024
'''

import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt

import pandas as pd 
import numpy as np 

regs = ['dlpfc', 'anterior-cingulate', 'amygdala', 'hippocampus', \
        'caudate', 'putamen', 'nucleus-accumbens', 'cerebellar-hemisphere']

## paths 
main_path = '/Users/nhunghoang/Desktop/remote_platypus/paper_twas'
twas_path = f'{main_path}/outputs_UKB/replication_twas.csv' 

colr_path = f'{main_path}/aux_files/color_dict.txt' 

## load colors 
df = pd.read_table(colr_path, index_col='label')
colors = df.to_dict()['hex']

## load twas data 
repl_sizes = ['HCP', '772', '2000', '5000', '15000'] 
data = pd.read_csv(twas_path)

## init plots 
plt.ion()
fig, axes_ = plt.subplots(4, 2, figsize=(7,8))
axes = {r:a for r,a in zip(regs, axes_.flatten())}

## loop
xpos = np.arange(len(repl_sizes)) + 1
fly_props = {'marker': '.', 'markersize': 3, 'markeredgecolor': 'none', 'alpha': 1}

for reg, rdf in data.groupby('vol'): 
    ax = axes[reg] 
    mc = colors[reg] 

    df = rdf.pivot(index='itr', columns='repl_set', values='perc_repl')
    df = df[repl_sizes]

    bp = ax.boxplot(df.values, positions=xpos, patch_artist=True, \
                    showcaps=False, showmeans=False, flierprops=fly_props)
    _ = [box.set_facecolor(mc) for box in bp['boxes']]
    _ = [box.set_alpha(0.5) for box in bp['boxes']]
    _ = [box.set_edgecolor('none') for box in bp['boxes']]
    _ = [whk.set_color('none') for whk in bp['whiskers']]
    _ = [med.set_color('none') for med in bp['medians']]

    be = ax.boxplot(df.values, positions=xpos, patch_artist=True, flierprops=fly_props)
    _ = [box.set_facecolor('none') for box in be['boxes']]
    _ = [med.set_color('k') for med in be['medians']]

for ax in axes_.flatten(): 
    ax.set_xticks(xpos ,size=8)
    ax.set_xticklabels(repl_sizes, size=8)
    ax.set_xlabel('size of replication set', size=9)

    ax.set_ylim([-0.02, 0.9])

plt.suptitle('y-axis: fraction of discovery genes that were replicated', size=10)
plt.tight_layout()
plt.savefig('repl_plot.pdf', format='pdf')

