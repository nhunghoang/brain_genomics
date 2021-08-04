'''
Save observed-CCC simulated annealing results (from 100 runs) into a 
cohesive file. Then analyze gene selection results from across the runs. 

- Nhung, June 2021 
'''

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt 
import os 
import numpy as np 
import h5py 

results_dir = '/data1/rubinov_lab/brain_genomics/analyses_HCP/multi_gene_assoc/results_ccc_obs32'
plots_dir_out = '/data1/rubinov_lab/brain_genomics/analyses_HCP/multi_gene_assoc/plots_ccc32'

## gather common genes per region pair 
cgenes_file = '/data1/rubinov_lab/brain_genomics/analyses_HCP/multi_gene_assoc/pairwise_common_genes.txt' 
common_genes = {} ## k: reg1_reg2, v: array of common genes 
with open(cgenes_file, 'r') as f: lines = f.readlines()
for line in lines: 
    [reg_pair, genes] = line.split(':')
    gene_arr = np.array(genes.strip().split(','))
    common_genes[reg_pair] = gene_arr

## parse log files 
CCC_rp = {} ## k: region pair, v: CCC array of len 100 
nsel_rp = {} ## k: region pair, v: nsel array of len 100 
for res_subdir in os.listdir(results_dir): 
    if res_subdir[:4] != 'logs': continue 

    region_pair = res_subdir[5:] 
    CCC = []; nsel = [] 
    for r in range(100): 
        log_file = '{}/{}/{}.log'.format(results_dir, res_subdir, r) 
        with open(log_file, 'r') as f: 
            result = f.readlines()[1].strip().split(' ')
            nsel.append(int(result[0]))
            CCC.append(float(result[4]))

    CCC_rp[region_pair] = np.array(CCC)
    nsel_rp[region_pair] = np.array(nsel)
             
## parse weight files, gather genes selected 
selected_genes_flat = {} ## k: reg1_reg2, v: list of selected genes across all runs  
selected_genes = {} ## k: reg1_reg2, v: list of selected genes per run 
for res_subdir in os.listdir(results_dir): 
    if res_subdir[:7] != 'weights': continue 

    region_pair = res_subdir[8:]
    selc_genes_flat = [] 
    selc_genes = [] 
    for r in range(100): 
        wgt_file = '{}/{}/{}.hdf5'.format(results_dir, res_subdir, r) 
        with h5py.File(wgt_file, 'r') as f: 
            W = np.array(f['W']).astype(bool)
        genes = common_genes[region_pair][W]
        selc_genes.append(genes) 
        selc_genes_flat.extend(genes) 
    selected_genes[region_pair] = selc_genes    
    selected_genes_flat[region_pair] = selc_genes_flat

## distribution statistics  
best_select = {} ## key: region pair, val: genes selected in best model 
for reg_pair in CCC_rp.keys(): 
    
    [reg1, reg2] = reg_pair.split('_')

    CCC_distrib = CCC_rp[reg_pair] 
    nsel_distrib = nsel_rp[reg_pair] 
    
    best_idx = np.argmax(CCC_distrib) 
    best_CCC = CCC_distrib[best_idx] 
    best_nsel = nsel_distrib[best_idx] 
    best_select[reg_pair] = selected_genes[reg_pair][best_idx] 

    mean_CCC = np.mean(CCC_distrib) 
    mean_nsel = np.mean(nsel_distrib) 

    std_CCC = np.std(CCC_distrib) 
    std_nsel = np.std(nsel_distrib) 

    n_total = len(common_genes[reg_pair])
    print('{} & {}:'.format(reg1, reg2))
    print('best CCC: {:.3f} ({} / {} genes selected)'.format(best_CCC, best_nsel, n_total))
    print('mean CCC: {:.3f} ± {:.3f}'.format(mean_CCC, std_CCC))
    #print('mean nsel: {:.3f} ± {:.3f}\n'.format(mean_nsel, std_nsel))

## scatter plot: (num genes total, num genes selected); points are region pairs 
'''
reg_pairs = list(common_genes.keys())
n_total = [common_genes[rp].size for rp in reg_pairs]
n_selec = [best_select[rp].size for rp in reg_pairs] 

fig, ax = plt.subplots(1,1,figsize=(15,15))
ax.scatter(n_total, n_selec, c='k')
ax.set_xlabel('number of total genes', fontsize=30)
ax.set_ylabel('number of selected genes', fontsize=30)
ax.tick_params(labelsize=30)

xleft, xright = ax.get_xlim()
ybottom, ytop = ax.get_ylim()
ax.set_aspect(abs((xright-xleft)/(ybottom-ytop)))

title = 'Proportion of Genes Selected in Best CCC Models' 
ax.set_title(title, size=30)
fname = '{}/proportion_genes_selected.png'.format(plots_dir_out)
plt.savefig(fname)
plt.close('all')
'''

## histograms of selected genes 
# x: different genes, y: freq across 100 runs, diff color: in best model 
for reg_pair in selected_genes_flat.keys(): 
    
    [reg1, reg2] = reg_pair.split('_')

    data = selected_genes_flat[reg_pair] 
    chosen = best_select[reg_pair] 

    genes, counts = np.unique(data, return_counts=True) 
    colors = ['b' if (g in chosen) else 'y' for g in genes] 

    fig, ax = plt.subplots(1,1,figsize=(30,5))
    #n, bins, patches = ax.hist(data, color='y')  
    ax.scatter(genes, counts, c=colors) 
    
    
    '''
    values =  np.random.randint(51, 140, 1000)
    n, bins, patches = plt.hist(values, bins=np.arange(50, 140, 2), align='left', color='g')
    patches[40].set_fc('r')
    plt.show()
    '''

    ax.set_xlabel('genes', fontsize=30) 
    ax.set_ylabel('frequency', fontsize=30) 
    ax.set_xticklabels(genes, rotation=15, ha="right")

    title = '{} & {}: Gene Selection Across 100 Runs'.format(reg1, reg2)
    ax.set_title(title, size=30) 
    fname = '{}/gene-freq_{}.png'.format(plots_dir_out, reg_pair)
    plt.savefig(fname) 
    plt.close('all') 
