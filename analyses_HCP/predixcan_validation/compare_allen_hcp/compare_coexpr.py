'''
Compare co-expression of same genes across 
average Allen (ABI) and HCP individuals. 

- Nhung, updated Feb 2022
'''

import os 
import numpy as np 
import h5py 
from scipy.stats import pearsonr, spearmanr
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

tag = '_standardize' 

## paths 
hcp_path = '/data1/rubinov_lab/brain_genomics/data_HCP/expression'
abi_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/predixcan_validation/compare_allen_hcp/allen_expr{}.hdf5'.format(tag)

## gather region names (in abc order)  
rand_regs_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/assoc/null_pvals_alff'
regions = sorted(os.listdir(rand_regs_path))
#regions.remove('hypothalamus')
#regions.remove('substantia-nigra')

## gather HCP regional expression and genes
hcp_expr = {} ## k: reg, v: expr (subjs * genes)   
hcp_gene = {} ## k: reg, v: genes 
for reg in regions: 
    with h5py.File('{}/{}.hdf5'.format(hcp_path, reg), 'r') as f: 
        hcp_expr[reg] = np.array(f['pred_expr']).T
        genes = [g.decode().split('.')[0] for g in np.array(f['genes'])]
        hcp_gene[reg] = np.array(genes) 

## gather Allen regional expression and genes (and subject num of sites per reg)  
with h5py.File(abi_path, 'r') as f: 
    abi_subj = np.array([s.decode() for s in np.array(f['subjects'])])
    abi_gene = np.array([g.decode() for g in np.array(f['genes'])])
    abi_expr = {reg:np.array(f[reg]) for reg in regions} ## k: reg, v: expr (subjs * genes) 
    abi_site = {(s,r): np.array(f['{}-{}-sites'.format(s,r)]) \
                for s in abi_subj for r in regions} ## k: (subj,reg), v: num sites 

## function: compute subject-wise co-expression 
## source: https://tinyurl.com/wd7vvt7t 
def subj_coexpr(reg_expr1, reg_expr2):

    ## spearman version 
    return np.array([spearmanr(r1,r2)[0] for r1,r2 in zip(reg_expr1,reg_expr2)]) 
    
    ## mean center 
    c1 = reg_expr1 - reg_expr1.mean(axis=1)[:,None]
    c2 = reg_expr2 - reg_expr2.mean(axis=1)[:,None]
    ## sum of squares 
    s1 = (c1**2).sum(axis=1)
    s2 = (c2**2).sum(axis=1)
    ## Pearson's coefficients 
    top = np.dot(c1, c2.T)
    bot = np.sqrt( np.dot(s1[:,None], s2[None]) ) 
    return np.diag(top/bot)

## initialize pairwise subplots for HCP/Allen co-expression 
fig_hcp, axes_hcp = plt.subplots(5,9, figsize=(40,25)) ## 45 region pairs 
fig_abi, axes_abi = plt.subplots(5,9, figsize=(40,25)) ## 45 region pairs 
acc = 0 ## use to compute subplot coordinates 

## compute co-expression  
hcp_coexpr = [] 
abi_coexpr = [] 
for i,reg1 in enumerate(regions): 
    for reg2 in regions[i+1:]:

        #if 'cerebellar-hemisphere' in [reg1,reg]: continue  

        ## find pairwise common genes available in HCP and Allen 
        g1 = hcp_gene[reg1]; g2 = hcp_gene[reg2] 
        hcp_common, idx1, idx2 = np.intersect1d(g1, g2, return_indices=True)  

        common_genes = np.intersect1d(hcp_common, abi_gene) 
        hcp_idx = np.in1d(hcp_common, common_genes) 
        abi_idx = np.in1d(abi_gene, common_genes)  

        ## compute HCP co-expression 
        hcp_expr1 = hcp_expr[reg1][:,idx1][:,hcp_idx] ## (subjs * genes) 
        hcp_expr2 = hcp_expr[reg2][:,idx2][:,hcp_idx]  
        hcp_subj_coexpr = subj_coexpr(hcp_expr1, hcp_expr2)
        hcp_coexpr.append(hcp_subj_coexpr.mean())
        
        ## compute Allen co-expression 
        abi_expr1 = abi_expr[reg1][:,abi_idx] ## (subjs * genes) 
        abi_expr2 = abi_expr[reg2][:,abi_idx]  
        abi_subj_coexpr = subj_coexpr(abi_expr1, abi_expr2)
        abi_coexpr.append(abi_subj_coexpr.mean())

        ## normalize by region "size"?
        #sites1 = [abi_site[(s,reg1)] for s in abi_subj]
        #sites2 = [abi_site[(s,reg2)] for s in abi_subj]
        #normed = abi_subj_coexpr / (np.sqrt(sites1) * np.sqrt(sites2))
        #abi_coexpr.append(normed.mean())

        ## co-expression subplots
        ## axes: regs, points: average subject gene expression 
        row = int(acc/9); col = acc % 9
        axes_hcp[row][col].scatter(hcp_expr1.mean(axis=0), hcp_expr2.mean(axis=0), c='b', s=0.6, alpha=0.35)
        axes_abi[row][col].scatter(abi_expr1.mean(axis=0), abi_expr2.mean(axis=0), c='g', s=0.6, alpha=0.35)

        axes_hcp[row][col].set_title('{:.3f}'.format(hcp_coexpr[-1]), size=30) 
        axes_abi[row][col].set_title('{:.3f}'.format(abi_coexpr[-1]), size=30)

        for ax in [axes_hcp[row][col], axes_abi[row][col]]: 
            ax.set_xlabel(reg1[:3], size=20) 
            ax.set_ylabel(reg2[:3], size=20) 
            xleft, xright = ax.get_xlim()
            ybottom, ytop = ax.get_ylim()
            ax.set_aspect(abs((xright-xleft)/(ybottom-ytop)))
        acc += 1 

fig_hcp.tight_layout(); fig_abi.tight_layout() 
fig_hcp.savefig('plot_HCP_coexpr_spearman{}.png'.format(tag))
fig_abi.savefig('plot_ABI_coexpr_spearman{}.png'.format(tag))
plt.close('all')
        
rho, pval = spearmanr(hcp_coexpr, abi_coexpr)
print('{:.3f} ({:.5f})'.format(rho, pval))

## scatter plot HCP and Allen co-expression 
fig, ax = plt.subplots(1,1)
ax.scatter(hcp_coexpr, abi_coexpr, c='k') 

ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

ax.tick_params(labelsize=10, length=5)
ax.set_xlabel('HCP co-expression', size=10)
ax.set_ylabel('Allen co-expression', size=10)

ax.set_title('r = {:.3f} (p ≤ {:.3f})'.format(rho, pval), size=15)

xleft, xright = ax.get_xlim()
ybottom, ytop = ax.get_ylim()
ax.set_aspect(abs((xright-xleft)/(ybottom-ytop)))

plt.savefig('plot_coexpr_spearman{}.png'.format(tag))
plt.close('all')
