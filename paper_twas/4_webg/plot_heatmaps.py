'''

Plot phecat * vol heatmap 
and gene * tissue heatmap

'''

import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd
import seaborn as sns
import sys

from mpl_toolkits.axes_grid1 import make_axes_locatable
from sctriangulate.colors import build_custom_continuous_cmap as BCCC

## params
ptype = sys.argv[1] ## FDR, BON
pink_cmap = BCCC([255, 247, 247], [254, 180, 198], [255, 52, 104], [149, 53, 78])

## paths 
main_path = '/Users/nhunghoang/Desktop/remote_platypus/paper_twas'
gpos_path = f'{main_path}/aux_files/gene_positions.bed'
colr_path = f'{main_path}/aux_files/color_dict.txt'

gxxt_path = f'{main_path}/outputs_UKB/pdx_overlap/{ptype}_table_genxtis.csv'
pxxv_path = f'{main_path}/outputs_UKB/pdx_overlap/{ptype}_table_phexvol.csv'

## out paths 
plot_path = f'{main_path}/outputs_UKB/plots/fig_s3ab.pdf'

## regions
regs = ['dlpfc', 'anterior-cingulate', 'amygdala', 'hippocampus', \
        'caudate', 'putamen', 'nucleus-accumbens', 'cerebellar-hemisphere']
name = ['DLPFC', 'ant. cingulate', 'amygdala', 'hippocampus', \
        'caudate', 'putamen', 'nuc. accumbens', 'cerebellar hemi.']
rmap = {r: n for r, n in zip(regs, name)}

##############################################################

## load gene positions
pos = pd.read_table(gpos_path, usecols=['chrom', 'gene_id', 'rank'])
rankings = dict(zip(pos['gene_id'], pos['rank']))

## load regional colors 
cdf = pd.read_table(colr_path, index_col='label')
reg_color = {r: cdf.loc[r]['hex'] for r in regs}

##############################################################

## load heatmap data 
genxtis = pd.read_csv(gxxt_path, index_col='sym')[regs]
phexvol = pd.read_csv(pxxv_path, index_col='phecate')[regs]

##############################################################

## init plot 
plt.ion()
fig, axes_ = plt.subplots(1, 2, figsize=(6, 24))
axes = {'phe': axes_[0], 'gen': axes_[1]}

daxes = {}
for key, ax in axes.items(): 
    div = make_axes_locatable(ax)
    dax = div.append_axes('top', size='1%', pad=0.05)
    daxes[key] = dax

kws = {'square': True, 'fmt': '', 'linewidth': 0.5, 'linecolor': '#E5E4E2', \
       'cbar_kws': {'orientation': 'horizontal'}}
tsize = 4
if ptype == 'BON': tsize = 8

##############################################################

## sort phecats by sum count
tot = phexvol[regs].sum(axis=1)
idx = np.argsort(tot)[::-1]
phesort = phexvol.index[idx]
phexvol = phexvol.reindex(phesort)

## init plot for phecat * volume
ax = axes['phe']
data = phexvol.values.astype(int)
vlog = np.where(data > 0, np.log10(data), -0.5)
vtxt = np.where(data > 0, data, '')

## band-aid: add faux squares to make 
## same size as genxtis mat
diff = genxtis.shape[0] - phexvol.shape[0]
vlog = np.concatenate([vlog, np.zeros((diff, 8))])
vtxt = np.concatenate([vtxt, np.zeros((diff, 8), dtype=str)])
mask = np.zeros_like(vlog, dtype=bool)
mask[vlog.shape[0]:] = True 

## plot phecat * volume
sns.heatmap(vlog, ax=ax, annot=vtxt, mask=mask, \
            cmap='Oranges', cbar_ax=daxes['phe'], \
            annot_kws={'size': tsize}, **kws)

## ax format
yformat = {'axis': 'y', 'labelright': True, 'right': True, 'labelleft': False, 'left': False}
ax.tick_params(size=0, pad=2, **yformat)

ax.set_yticks(np.arange(phesort.size) + 0.5)
ax.set_yticklabels(phesort, rotation=0, size=tsize)

ax.set_xticks(np.arange(8) + 1)
ax.set_xticklabels(regs, size=tsize, rotation=90, ha='right')

## cbar format 
cax = ax.collections[0].colorbar.ax
cax.tick_params(axis='x', top=True, labeltop=True, \
                bottom=False, labelbottom=False, size=2)

cax.set_xlim([0, cax.get_xlim()[1]])

if ptype == 'BON': xticks = [1, 4, 16] 
else: xticks = [1, 4, 16, 64]

cax.set_xticks(np.log10(xticks))
cax.set_xticklabels(xticks, size=tsize)

##############################################################

## sort genes by sum count
tot = genxtis[regs].sum(axis=1)
idx = np.argsort(tot)[::-1]
gensort = genxtis.index[idx]
genxtis = genxtis.reindex(gensort)

## bandaid: remove last two rows (zero counts) 
#if ptype == 'FDR': genxtis = genxtis.drop(['PPP2CA', 'LTBP1'])
#elif ptype == 'bonf': genxtis = genxtis.drop(['ANKS1A'])
#gensort = genxtis.index

## plot gene * tissue
ax = axes['gen']
data = genxtis.values.astype(int)
vlog = np.where(data > 0, np.log10(data), -0.5)
vtxt = np.where(data > 0, data, '')

sns.heatmap(vlog, ax=ax, annot=vtxt, \
            cmap=pink_cmap, cbar_ax=daxes['gen'], \
            annot_kws={'size': tsize}, **kws)

## ax format
ax.tick_params(size=0, pad=2)
ax.set_yticks(np.arange(gensort.size) + 0.5)
ax.set_yticklabels(gensort, rotation=0, size=tsize)

ax.set_xticks(np.arange(8) + 1)
ax.set_xticklabels(regs, size=tsize, rotation=90, ha='right')

## cbar format 
cax = ax.collections[0].colorbar.ax
cax.tick_params(axis='x', top=True, labeltop=True, \
                bottom=False, labelbottom=False, size=2)

cax.set_xlim([0, cax.get_xlim()[1]])

if ptype == 'BON': xticks = [1, 2, 4, 8] 
else: xticks = [1, 4, 16]


cax.set_xticks(np.log10(xticks))
cax.set_xticklabels(xticks, size=tsize)

##############################################################

## save 
plt.savefig(plot_path, format='pdf')




