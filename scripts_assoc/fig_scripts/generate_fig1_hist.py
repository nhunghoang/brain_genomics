'''
Figure 1C: histogram of PrediXcan model p-values 

- Nhung, updated Sept 2022
'''

import numpy as np 
import h5py 
import matplotlib.pyplot as plt 

regs = ['frontal-pole', 'anterior-cingulate', 'caudate', 'putamen', 'nucleus-accumbens', \
        'hippocampus', 'amygdala', 'hypothalamus', 'substantia-nigra', 'cerebellar-hemisphere'] 
rcolors = {'frontal-pole': '#EC7F82', \
           'anterior-cingulate': '#F1A043', \
           'caudate': '#9FCF63', \
           'putamen': '#6FAE61', \
           'nucleus-accumbens': '#F2C164', \
           'hippocampus': '#B176F7', \
           'amygdala': '#5F2F8D', \
           'hypothalamus': '#77B3D6', \
           'substantia-nigra': '#2F6EBA', \
           'cerebellar-hemisphere': '#FE3A3E'}

## get plot data 
file_r2s = '/data1/rubinov_lab/brain_genomics/models_PrediXcan_v8/all_r2_distribution.hdf5'
file_pvs = '/data1/rubinov_lab/brain_genomics/models_PrediXcan_v8/all_pv_distribution.hdf5'

with h5py.File(file_pvs, 'r') as f: 
    reg_pvs = {} ## k: reg, v: -log10(pvals) 
    for reg in regs: 
        pvs = f[reg][()] ## (all PrediXcan models,) 
        reg_pvs[reg] = -np.log10(pvs) 

## plot histogram and inset 
plt.ion() 
fig, ax = plt.subplots(1, 1, figsize=(5,3))

bounds = [18, 260, 35, 500]
ins = ax.inset_axes(bounds, transform=ax.transData)

## define histogram bins 
_, bin_edges = np.histogram(np.arange(63), bins=60)
xdata = bin_edges[:-1] + (bin_edges[1]-bin_edges[0])/2
for x in xdata:
    ax.vlines(x, 0, 1000, linewidth=0.3, color='k', alpha=0.1, zorder=1)
    ins.vlines(x, 0, 1000, linewidth=0.3, color='k', alpha=0.1, zorder=1)

## vertical line at p=0.01 threshold 
vl = ax.vlines(-np.log10(0.01), 0, 1000, linewidth=1, color='k', zorder=2)
ax.text(3, 860, '$p=0.01$', fontsize=12, weight='bold', zorder=3)

## regional histograms 
for reg, color in rcolors.items(): 
    hdata, _ = np.histogram(reg_pvs[reg], bins=bin_edges)
    a0 = ax.scatter(xdata, hdata, s=10, color=color, zorder=4)  
    i0 = ins.scatter(xdata, hdata, s=10, color=color, zorder=4) 

## set plot limits and margins 
ax.margins(0) 
ax.set_ylim([0,950])

ins.set_xlim([25,60])
ins.set_ylim([0,40])
    
## set ticks 
tick_size = 10

ax.set_xticks([10, 20, 30, 40, 50, 60]) 
ax.set_yticks([0, 200, 400, 600, 800])

ax.set_xticklabels(['10', '', '30', '', '50', ''], fontsize=tick_size)  
ax.set_yticklabels(['0', '', '400', '', '800'], fontsize=tick_size)

ins.set_xticks([30, 35, 40, 45, 50, 55, 60])
ins.set_yticks([0, 10, 20, 30, 40])

ins.set_xticklabels(['30', '', '', '45', '', '', '60'], fontsize=tick_size)
ins.set_yticklabels(['0', '', '20', '', '40'], fontsize=tick_size)

## set labels 
label_size = 12 

ax.set_xlabel("gene model prediction score ($-\log_{10}{p}$)", fontsize=label_size)
ax.set_ylabel("frequency (# gene models)", fontsize=label_size)

## indicate zoom 
patch, clines = ax.indicate_inset_zoom(ins, edgecolor='k')
clines[0].set_visible(True); clines[1].set_visible(False)
clines[2].set_visible(False); clines[3].set_visible(True)

## tighten and save  
plt.tight_layout()
plt.savefig('fig_plots/fig1_hist.pdf', format='pdf') 

