'''
Compare expression levels to number of 
sampling sites per region. 

- Nhung, March 2022 
'''

import os 
import numpy as np 
import h5py 
from scipy.stats import pearsonr, spearmanr
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

## paths 
abi_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/predixcan_validation/compare_allen_hcp/allen_expr_standardize.hdf5'

## gather region names (in abc order)  
rand_regs_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/assoc/null_pvals_alff'
regions = sorted(os.listdir(rand_regs_path))

## gather Allen regional expression and genes (and subject num of sites per reg)  
with h5py.File(abi_path, 'r') as f: 
    abi_subj = np.array([s.decode() for s in np.array(f['subjects'])])
    abi_gene = np.array([g.decode() for g in np.array(f['genes'])])
    abi_expr = {reg:np.array(f[reg]) for reg in regions} ## k: reg, v: expr (subjs * genes) 
    abi_site = {(s,r): np.array(f['{}-{}-sites'.format(s,r)]) \
                for s in abi_subj for r in regions} ## k: (subj,reg), v: num sites 

## gather data to plot 
coexpr = [] 
nsites = [] 
for reg1 in regions: 
    if reg1[:3] == 'cer': continue
    expr1 = abi_expr[reg1] ## (subjs * genes)
    reg_coexpr = [] 
    for reg2 in regions: 
        if reg2[:3] == 'cer': continue
        expr2 = abi_expr[reg2]
        subj_coexpr = [spearmanr(expr1[s], expr2[s])[0] for s in range(abi_subj.size)]
        reg_coexpr.append(np.array(subj_coexpr).mean()) 
    coexpr.append(np.array(reg_coexpr).sum())
    nsites.append(np.mean([abi_site[(s,reg1)] for s in abi_subj]))

## scatter plot  
fig, ax = plt.subplots(1,1)
ax.scatter(nsites, coexpr, c='k') 

ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

ax.tick_params(labelsize=10, length=5)
ax.set_xlabel('average number of sites', size=10)
ax.set_ylabel('average co-expression sum', size=10)

xleft, xright = ax.get_xlim()
ybottom, ytop = ax.get_ylim()
ax.set_aspect(abs((xright-xleft)/(ybottom-ytop)))

plt.savefig('allen_site_bias')
plt.close('all')
