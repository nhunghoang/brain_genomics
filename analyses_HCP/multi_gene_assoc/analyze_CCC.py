'''
'''

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt 
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import numpy as np 
import h5py
import sys 
import os 

best_models_file = '/data1/rubinov_lab/brain_genomics/analyses_HCP/multi_gene_assoc/ccc_best_models.hdf5'

## Parse the results of the CCC observed and null model analyses 
## Save the CCC and corresponding weight sets of the best observed and null shuffle runs  
## hdf5 key format: reg1-reg2-[ccc/weights]-[observed/null]
def write_best_models():

    ## organize observed correlations 
    obsv_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/multi_gene_assoc/results_ccc_obs32' 
    obsv_cccs = {} ## k: region pair, v: best correlation    
    obsv_wgts = {} ## k: region pair, v: corresponding weight set 

    ## loop through region pairs (100 observed runs per pair)  
    for log_dir in os.listdir(obsv_path): 
        if log_dir[:4] != 'logs': continue  
        cccs = np.zeros(100, dtype=float)  
        for i in range(100): 
            log_path = '{}/{}/{}.log'.format(obsv_path, log_dir, i)
            with open(log_path, 'r') as f: log = f.readlines()  
            log_ccc = float(log[1].split(' ')[5])
            cccs[i] = log_ccc 

        ## select the best observed run 
        reg_pair = log_dir[5:] ## format: reg1_reg2
        best_idx = np.argmax(cccs) 
        wgt_path = '{}/weights_{}/{}.hdf5'.format(obsv_path, reg_pair, i)
        with h5py.File(wgt_path, 'r') as f: best_weights = np.array(f['W'])

        obsv_cccs[reg_pair] = cccs[best_idx] 
        obsv_wgts[reg_pair] = best_weights  

    ## organize null correlations 
    null_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/multi_gene_assoc/results_ccc_null32' 
    null_cccs = {} ## k: region pair, v: best correlation per shuffle  
    null_wgts = {} ## k: region pair, v: corresponding weight sets  

    sh_rn = [] 

    ## loop through region pairs (100 shuffles x 100 runs)  
    for reg_pair in obsv_cccs.keys():
        log_dir = 'logs_{}'.format(reg_pair)
        cccs = np.zeros((100,100), dtype=float) ## rows are shuffles, cols are runs  
        for sh in range(100): 
            for rn in range(100):
                log_path = '{}/{}/{}_{}.log'.format(null_path, log_dir, sh, rn)
                try:
                    with open(log_path, 'r') as f: log = f.readlines()  
                except: 
                    print(log_path)
                    sh_rn.append((sh*100) + rn)
                    continue 
                log_ccc = float(log[1].split(' ')[5])
                cccs[sh][rn] = log_ccc 

        ## select the best observed run per shuffle  
        best_idxs = np.argmax(cccs, axis=1) 
        best_weights = [] 
        for sh, rn in enumerate(best_idxs):
            wgt_path = '{}/weights_{}/{}_{}.hdf5'.format(null_path, reg_pair, sh, rn)
            with h5py.File(wgt_path, 'r') as f: best_weights.append(np.array(f['W']))

        null_cccs[reg_pair] = np.max(cccs, axis=1) 
        null_wgts[reg_pair] = best_weights  

    for x in np.unique(sh_rn):
        print(x)

    ## save models 
    with h5py.File(best_models_file, 'w') as f: 
        for reg_pair in obsv_cccs.keys():
            f[reg_pair + '.ccc-observed'] = obsv_cccs[reg_pair]
            f[reg_pair + '.weights-observed'] = obsv_wgts[reg_pair]
            f[reg_pair + '.ccc-null'] = null_cccs[reg_pair]
            f[reg_pair + '.weights-null'] = null_wgts[reg_pair]

#############################################################################################################

def main(): 

    #write_best_models()
    
    ## read best models 
    obsv_cccs = {}; obsv_wgts = {} 
    null_cccs = {}; null_wgts = {} 
    with h5py.File(best_models_file, 'r') as f: 
        for key in f.keys():
            data = np.array(f[key])
            reg_pair, model = key.split('.')
            if model == 'ccc-observed': obsv_cccs[reg_pair] = data
            elif model == 'weights-observed': obsv_wgts[reg_pair] = data
            elif model == 'ccc-null': null_cccs[reg_pair] = data
            elif model == 'weights-null': null_wgts[reg_pair] = data

    reg_pairs = list(obsv_cccs.keys())

    ## create box plots (two, dividing num of reg pairs in half) 
    def plot_box(RPs, count):

        fig_size = (50,10)
        title_size = 30
        label_size = 30
        tick_size = 20
        legend_size = 25

        null = np.array([null_cccs[rp] for rp in RPs]).T 
        obsv = [obsv_cccs[rp] for rp in RPs]
        fig, ax = plt.subplots(1,1,figsize=fig_size)

        bp = ax.boxplot(null)
        for median in bp['medians']:
            median.set(color='black')
        for flier in bp['fliers']:
            flier.set(marker='o', markersize=2, markerfacecolor='k')

        ax.scatter(np.arange(count)+1, obsv, marker='x', color='r')

        title = '[{}] CCC Best Simulated Annealing Models'.format(count)
        #xlabels = [r1+' - '+r2 for i,r1 in enumerate(REGIONS_SHORT) for r2 in REGIONS_SHORT[i+1:]]
        xlabels = RPs
        ax.set_title(title, size=title_size)
        ax.set_xticklabels(xlabels, size=tick_size, rotation=10, ha='right')
        ax.set_ylabel('Best CCC', size=label_size)
        ax.tick_params(axis='y', labelsize=tick_size)

        colors = ['red', 'black']
        lines = [Line2D([0], [0], color=c, linewidth=2, linestyle='-') for c in colors]
        labels = ['observed', 'null']
        plt.legend(lines, labels, loc=2, prop={'size':legend_size})

        fname = '/data1/rubinov_lab/brain_genomics/analyses_HCP/multi_gene_assoc/plots_ccc32/null-boxplot{}.png'.format(count)
        plt.savefig(fname)
        plt.close('all')

    plot_box(reg_pairs[:23], 23)
    plot_box(reg_pairs[23:], 22)

main()
