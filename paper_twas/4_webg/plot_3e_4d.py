'''
Scatter plot of enrichment analyses (GWAS Catalog and PrediXVU). 

- Nhung, updated Feb 2024
''' 

import matplotlib 
matplotlib.use('tkagg') 
import matplotlib.pyplot as plt  

import numpy as np 
import pandas as pd 

import sys 

enric = sys.argv[1] ## gwas_catalog 
ptype = sys.argv[2] ## FDR, BON
try: 
    group = sys.argv[3] ## KB
    phens = sys.argv[4] ## vol_mean 
except: 
    group = 'UKB'
    phens = 'vol_mean'

regs = ['dlpfc', 'anterior_cingulate', 'amygdala', 'hippocampus', \
        'caudate', 'putamen', 'nucleus_accumbens', 'cerebellar_hemisphere']

name = ['DLPFC', 'ant. cingulate', 'amygdala', 'hippocampus', \
        'caudate', 'putamen', 'nuc. accumbens', 'cerebellar hemi.']

## paths 
main_path = '/Users/nhunghoang/Desktop/remote_platypus/paper_twas'
colr_path = f'{main_path}/aux_files/color_dict.txt'

if ptype == 'BON': tag = 'fig_s3c' #'fig_s2a'
else: tag = 'fig4'

webg_path = f'{main_path}/outputs_{group}/enrich_{enric}/{ptype}_{phens}/enrichment_results' ## _{reg}.txt
plot_path = f'{main_path}/outputs_{group}/plots/{tag}_{enric}_{phens}.pdf'

####################################################################################################

## load enrichment results 
edata = {} ## k: reg, v: df w/ annots and pvals 
cdict = {'geneSet': 'phecode', 'description': 'phename', 'FDR': 'FDR'}
for reg in regs: 

    rfile = f'{webg_path}_{reg}.txt'
    df = pd.read_table(rfile, sep='\t', usecols=cdict.keys())
    df = df.rename(columns=cdict)

    df = df.loc[~df['phename'].isna()]

    logf = np.where(df['FDR'] > 0, -np.log10(df['FDR']), -1)
    logf = np.where(logf == -1, logf.max(), logf)
    df['logf'] = logf

    df = df.sort_values('phecode')
    edata[reg] = df

####################################################################################################

## load plot colors 
colors = pd.read_table(colr_path, sep='\t', index_col='label')
colors = colors.to_dict()['hex']

## init plot 
plt.ion()
fig, ax = plt.subplots(1, 1, figsize=(3,3))

## regional loop 
y0 = 0 
yticks = [] 
for r, reg in enumerate(regs[::-1]): 
    
    mcolr = colors[reg.replace('_', '-')]
    xdata = edata[reg]['logf']
    yrnge = np.arange(y0, y0 + xdata.size, 1) + (r*1.5)

    ## scatter plot 
    sc = ax.scatter(xdata, yrnge, c=mcolr, s=20, alpha=0.25, edgecolor='none')
    ec = ax.scatter(xdata, yrnge, c='none', s=20, edgecolor=mcolr, linewidth=0.5) 

    ## update yticks
    yt = y0 + (xdata.size / 2)
    y0 += xdata.size + 5 #100
    yticks.append(yt)

    ## print sig annots (for labeling) 
    mask = (xdata > 2.3)
    desc = edata[reg]['phename'][mask] 
    logf = xdata[mask] 
    print('\n' +  reg) 
    for f, d in zip(logf, desc): 
        print('{:.2f} {}'.format(f,d))

## set ax limits 
if enric == 'gwas_catalog': 
    ax.set_xlim([0, 14])
    ax.set_ylim([-100, 337])

    if ptype == 'BON':
        ax.set_xlim([0, 14])
        ax.set_ylim([-45, 1620])

if enric == 'pdx_nom0.001':
    ax.set_xlim([-0.1, 11])
    ax.set_ylim([-20, 675])

    if ptype == 'BON':
        ax.set_xlim([-0.1, 11])
        ax.set_ylim([-8, 220])

xlim = ax.get_xlim()
ylim = ax.get_ylim()

## set pval lines 
h0 = ax.vlines(-np.log10(0.01), ylim[0], ylim[1], linewidth=0.8, linestyle='solid', edgecolor='k')
h1 = ax.vlines(-np.log10(0.05), ylim[0], ylim[1], linewidth=0.5, linestyle='solid', edgecolor='k')

## set ax ticks
ax.set_yticks(yticks)
ax.set_yticklabels(name[::-1], size=7, ha='right', rotation=0)
## yticks = yticks[:3] + [96, 120, 144, 170, 197]; ax.set_yticks(yticks)

ax.set_xlabel('-log10(FDR)', size=7)
ax.set_title('{} {}'.format(group, phens), size=7)

## square
xleft, xright = ax.get_xlim()
ybottom, ytop = ax.get_ylim()
ax.set_aspect(abs((xright-xleft)/(ybottom-ytop)))

plt.tight_layout() 
plt.savefig(plot_path, format='pdf')

