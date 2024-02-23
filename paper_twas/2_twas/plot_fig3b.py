'''
Manhattan plots (Fig 3B and S1). 

- Nhung, Feb 2024
'''

import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd
import sys

from scipy.stats import spearmanr

## regions 
regs = ['dlpfc', 'anterior-cingulate', 'amygdala', 'hippocampus', \
        'caudate', 'putamen', 'nucleus-accumbens', 'cerebellar-hemisphere']
regs = ['dlpfc_psychencode']

## args 
group = sys.argv[1] ## UKB 
phens = sys.argv[2] ## vol_mean

main_path = '/Users/nhunghoang/Desktop/remote_platypus/paper_twas'

gpos_path = f'{main_path}/aux_files/gene_positions.bed'
colr_path = f'{main_path}/aux_files/color_dict.txt'

#tgwa_path = f'{main_path}/outputs_UKB/TWAS_GWAS_summary.csv'
#sigs_path = f'{main_path}/outputs_UKB/TWAS_GWAS_sigbars.csv'
tgwa_path = f'{main_path}/outputs_UKB/TWAS_GWAS_dlpfcPE_summary.csv'
sigs_path = f'{main_path}/outputs_UKB/TWAS_GWAS_dlpfcPE_sigbars.csv'

plot_path = f'{main_path}/outputs_UKB/plots/fig3b_s1'

##############################################################

## manh: gene positions 
pos = pd.read_table(gpos_path, usecols=['chrom', 'gene_id', 'rank'])
rankings = dict(zip(pos['gene_id'], pos['rank']))
chr_dict = dict(zip(pos['gene_id'], pos['chrom']))

## scat: UKB JTI summary table
table = pd.read_csv(tgwa_path).dropna()

## load sig data
sig_data = pd.read_csv(sigs_path, index_col='phenotype').to_dict()

##############################################################

## fig size 
#fsize = (7, 9.6) #(7, 1.2) 
fsize = (7, 1.2) 

## tick params
tsize = 7 
tstep = {r: 4 for r in regs} 
tstep['hippocampus'] = 8 

## marker params
markr = '.'
msize = 10 
alpha = {'sig': 0.3, 'not': 0.05, 'enot': 0.5}
lw = 0.3

## base colors 
ch_colors = ['#cccccc', '#a7a7a7']

## regional colors 
mdf = pd.read_table(colr_path, index_col='label')
reg_color = {r: mdf.loc[r]['hex'] for r in regs}

## ymax 
ymax = np.ceil(-np.log10(table[['phenotype', 'pval_GWAS', 'pval_TWAS']] \
               .groupby('phenotype').min()) + 1).max(axis=1).to_dict()

### init plot 
#plt.ion()
#fig = plt.figure(figsize=fsize)
#grid = plt.GridSpec(8, 3, figure=fig, width_ratios=[2/3, 1/6, 1/6])
#
### add axes 
#axes = {} 
#for r, reg in enumerate(regs):
#    axes[f'manh_{reg}'] = plt.subplot(grid[r,0])
#    axes[f'many_{reg}'] = plt.subplot(grid[r,1])
#    axes[f'ones_{reg}'] = plt.subplot(grid[r,2])

###################################################################

## func: plot manh 
def plot_manh(data, reg):
    ax = axes[f'manh_{reg}']
    ch_ticks = []
    start_idx = 0
    pmax = 0

    ## sig thresholds 
    benh = -np.log10(sig_data['TWAS_FDR'][reg])
    bonf = -np.log10(sig_data['TWAS_BON'][reg])

    ## edge marker kws 
    em_kws = {'color': 'none', 's': msize, 'marker': markr, \
              'alpha': alpha['enot'], 'linewidth': lw}

    ## face marker kws
    fm_kws = {'s': msize, 'marker': markr, 'edgecolor': 'none'}

    ## plot by chrom 
    for ch, cd in data:

        logp = cd['logp'].values 
        sig = (logp >= benh)
        xrange = np.arange(start_idx, start_idx + logp.size)

        ## TWAS not sig fill 
        sc0 = ax.scatter(xrange[~sig], logp[~sig], c=ch_colors[int(ch)%2], alpha=alpha['not'], **fm_kws)

        ## TWAS not sig edge 
        ec0 = ax.scatter(xrange[~sig], logp[~sig], edgecolor=ch_colors[int(ch)%2], **em_kws)

        ## TWAS sig fill
        sc1 = ax.scatter(xrange[sig], logp[sig], c=reg_color[reg], alpha=alpha['sig'], **fm_kws)

        ## TWAS sig edge 
        sc2 = ax.scatter(xrange[sig], logp[sig], edgecolor=reg_color[reg], **em_kws)

        ## chr ticks 
        ch_ticks.append(start_idx + (logp.size/2) - 1)
        start_idx += logp.size
        pmax = np.max([pmax, logp.max()])

    ## x-axis 
    ax.set_xlim(0, start_idx)
    #ax.set_xticks(ch_ticks)
    #ax.set_xticklabels([str(c) for c in range(1,23)], size=7)

    ## y-axis
    ax.set_ylim([0, ymax[reg]])
    ax.set_yticks(np.arange(0, ymax[reg], tstep[reg]))

    ## both axes 
    ax.tick_params(labelsize=tsize)
    ax.margins(x=0.005, y=0)

    ## plot sig lines
    ax.hlines(benh, 0, start_idx, color='k', linewidth=lw)
    ax.hlines(bonf, 0, start_idx, color='k', linewidth=lw, linestyle='dashed')

    return

###################################################################

## plot pval scat
def plot_pvals(data, reg, ax):
    colr = reg_color[reg]

    ## lines
    kws = {'color': 'k', 'linewidth': lw}

    gbon = -np.log10(sig_data['GWAS_BON'][reg]) 
    tbon = -np.log10(sig_data['TWAS_BON'][reg]) 
    tben = -np.log10(sig_data['TWAS_FDR'][reg]) 

    ax.vlines(gbon, 0, ymax[reg], linestyles='dashed', **kws)
    ax.hlines(tbon, 0, ymax[reg], linestyles='dashed', **kws)
    ax.hlines(tben, 0, ymax[reg], linestyles='solid', **kws)

    ## TWAS FDR bar
    tsig = (data['twas'] >= tben)

    ## edge marker kws
    em_kws = {'color': 'none', 's': msize, 'alpha': alpha['enot'], 'linewidth': lw, 'marker': markr}
    
    ## face marker kws
    fm_kws = {'s': msize, 'marker': markr, 'edgecolor':'none'}

    ## not TWAS FDR-sig 
    sc0 = ax.scatter(data['gwas'][~tsig], data['twas'][~tsig], c=ch_colors[1], alpha=alpha['not'], **fm_kws)
    ec0 = ax.scatter(data['gwas'][~tsig], data['twas'][~tsig], edgecolor=ch_colors[0], **em_kws) 

    ## TWAS FDR-sig
    sc1 = ax.scatter(data['gwas'][tsig], data['twas'][tsig], c=reg_color[reg], alpha=alpha['sig'], **fm_kws)
    ec2 = ax.scatter(data['gwas'][tsig], data['twas'][tsig], edgecolor=reg_color[reg], **em_kws) 

    ## limits 
    ax.set_xlim([0, ymax[reg]])
    ax.set_ylim([0, ymax[reg]])
    ax.margins(0)

    ## ticks 
    ax.set_xticks(np.arange(0, ymax[reg], tstep[reg]))
    ax.set_yticks(np.arange(0, ymax[reg], tstep[reg]))
    ax.tick_params(labelsize=tsize)


    xleft, xright = ax.get_xlim()
    ybottom, ytop = ax.get_ylim()
    ax.set_aspect(abs((xright-xleft)/(ybottom-ytop)))
    return

###################################################################

## main loop 
plt.ion()

reg_tables = table.groupby('phenotype')
for reg, rtable in table.groupby('phenotype'):

    ## init plot 
    fig = plt.figure(figsize=fsize)
    grid = plt.GridSpec(1, 3, figure=fig, width_ratios=[2/3, 1/6, 1/6])
    
    ## add axes 
    axes = {} 
    axes[f'manh_{reg}'] = plt.subplot(grid[0,0])
    axes[f'many_{reg}'] = plt.subplot(grid[0,1])
    axes[f'ones_{reg}'] = plt.subplot(grid[0,2])

    ## isolate twas 
    df = rtable[['gene', 'pval_TWAS']].drop_duplicates()
    df['logp'] = -np.log10(df['pval_TWAS'])

    ###### MANH ######

    ## add position info 
    df['rank'] = df['gene'].map(rankings)
    df['chrom'] = df['gene'].map(chr_dict)

    ## sort by gene position 
    df.sort_values('rank', inplace=True)
    df = df.loc[~df['rank'].isnull()]

    ## organize by chrom 
    cdata = df.groupby('chrom')

    ## plot manh 
    plot_manh(cdata, reg)

    ###### SCAT ######

    ## get the Spearman correlation of the entire table 
    all_rho, all_pvl = spearmanr(rtable['pval_GWAS'], rtable['pval_TWAS'])

    ## plot all pval points 
    all_data = {'gwas': -np.log10(rtable['pval_GWAS']), \
                'twas': -np.log10(rtable['pval_TWAS'])}
    plot_pvals(all_data, reg, axes[f'many_{reg}'])

    ## keep the best GWAS for every (TWAS) gene 
    pdf = rtable.sort_values('pval_GWAS', ascending=True).drop_duplicates('gene')

    ## case: multiple genes share the same best GWAS
    ## keep the best TWAS gene for every SNP
    pdf = pdf.sort_values('pval_TWAS', ascending=True).drop_duplicates('rsid')

    ## get the Spearman correlation of the remaining table
    top_rho, top_pvl = spearmanr(pdf['pval_GWAS'], pdf['pval_TWAS'])

    ## plot this scatter 
    top_data = {'gwas': -np.log10(pdf['pval_GWAS']), \
                'twas': -np.log10(pdf['pval_TWAS'])}
    plot_pvals(top_data, reg, axes[f'ones_{reg}'])

    ## standardize height 
    plt.tight_layout()
    sx, sy, sw, sh = axes[f'many_{reg}'].get_position().bounds
    mx, my, mw, mh = axes[f'manh_{reg}'].get_position().bounds
    mc_pos = [mx, sy, mw, sh]
    axes[f'manh_{reg}'].set_position(mc_pos)

    ## save 
    plt.savefig(f'{plot_path}_{reg}.pdf', format='pdf')

