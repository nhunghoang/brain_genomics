''' 
Analyze the prediction performances of the PrediXcan 
gene models, in terms of r^2 and p-values. 

Use analysis to decide quality thresholds for the initial 
set of genes per region (prior to phenotype correlation). 

In addition, look at the distribution of variance across HCP 
for all genes per tissue, differentiating between models with 
all SNPs vs missing (HCP) SNPs.  

- Nhung, June 2021 
'''

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os 
import h5py 
from statsmodels.stats.multitest import fdrcorrection 

## input and output paths 
data_dir = '/data1/rubinov_lab/brain_genomics/data_HCP/expression'
outs_dir = '/data1/rubinov_lab/brain_genomics/models_PrediXcan_v8/plots_model_quality'

## read genes per model that HCP doesn't have all SNPs for 
unkept_file = '/data1/rubinov_lab/brain_genomics/data_HCP/expression/all_genes_not_kept.txt'
with open(unkept_file, 'r') as f: lines = f.readlines() 
unkept_models = {} ## k: region, v: list of genes
for line in lines: 
    data = line.split('\t') 
    reg = data[0]
    mod = data[1] 
    try: unkept_models[reg].append(mod) 
    except KeyError: unkept_models[reg] = [mod] 

## r^2 and p-value thresholds of interest 
#r2_incr = np.arange(0.1, 1, 0.1) ## hm-big 
r2_incr = np.arange(0, 0.16, 0.025) ## hm-small
#pv_incr = np.arange(0.01, 0.051, 0.01)
pv_incr = np.array([10e-6, 10e-5, 10e-4, 10e-3, 10e-2])

## heatmap details 
r2_labels = ['≥' + str(round(r2,3)) for r2 in r2_incr]
pv_labels = ['≤' + str(round(pv,5)) for pv in pv_incr]
hm_kwargs = {'cmap':'plasma', 'xticklabels':r2_labels, 'yticklabels':pv_labels, \
            'fmt':'4d', 'square': True, 'annot':True, 'annot_kws':{'size':16}}
fig_size = (20,15)
title_size = 20
label_size = 22

## loop through regions 
for data_file in os.listdir(data_dir): 
    if data_file[-4:] != '.log': continue 
    with open('{}/{}'.format(data_dir, data_file), 'r') as f: 
        f.readline() ## header 
        lines = f.readlines() 
        n_models = len(lines) 
    reg = data_file.split('.')[0]
    print('- {}'.format(reg))
    
    ## gather name, r^2, and p-value per gene model 
    models = []
    r2s = np.zeros(n_models, dtype=float)
    pvs = np.zeros(n_models, dtype=float)
    for i,line in enumerate(lines): 
        data = line.strip().split('\t')
        models.append(data[0])
        r2s[i] = data[-2]
        pvs[i] = data[-1]

    models = np.array(models) 
    r2_shape = r2_incr.shape[0]
    pv_shape = pv_incr.shape[0]
    counts = np.zeros((pv_shape, r2_shape), dtype=int) 

    ## apply FDR correction to p-values 
    pvs_rejected, pvs_corrected = fdrcorrection(pvs) 

    ## scatter plot of r^2 vs p-value (uncorrected & corrected), with points representing models 
    '''
    pvs_log10 = np.log10(pvs) 
    pcs_log10 = np.log10(pvs_corrected)
    plt.scatter(r2s, pvs_log10, s=1, c='r', label='uncorrected')
    plt.scatter(r2s, pcs_log10, s=1, c='b', label='FDR-corrected')
    plt.ylim(top=0)
    '''

    ## add line for p ≤ 0.01 and p ≤ 10e-6
    '''
    plt.plot(r2s, [np.log10(0.01) for r in r2s], linewidth=0.5, c='k')
    plt.plot(r2s, [np.log10(10e-6) for r in r2s], linewidth=0.5, c='k')
    yticks_new = list(plt.yticks()[0]) + [np.log10(0.01), np.log10(10e-6)]
    ylabels_new = [str(a) for a in yticks_new[:-2]] + ['p ≤ 0.01', 'p ≤ 10e-6']
    plt.yticks(yticks_new, ylabels_new) 
    '''

    ## scatter plot text
    '''
    plt.title('{}\n({} total genes)'.format(reg.upper(), n_models), size=title_size-5, fontweight='bold') 
    plt.xlabel('r^2 values', fontsize=label_size-5) 
    plt.ylabel('log10(p-values)', fontsize=label_size-5) 
    plt.legend()
    plt.savefig('{}/sc-{}.png'.format(outs_dir, reg))
    plt.close('all')
    '''

    ## for this specific threshold pair, report yield and the number of models that are missing HCP SNPs
    '''
    r2_models = models[r2s >= 0.05]
    pv_models = models[pvs_corrected <= 0.01]
    common_models = np.intersect1d(r2_models, pv_models)
    incomplete_models = np.intersect1d(common_models, np.array(unkept_models[reg]))

    count = common_models.shape[0]
    inc_count = incomplete_models.shape[0] 
    pyield = (float(count) / n_models) * 100
    print('{:>25s}: {} / {} = {:.2f}% yield, includes {} incomplete models'.format(reg, count, n_models, pyield, inc_count))
    ''' 

    ## count the number of models that meet every pair of (r^2,p) thresholds 
    '''
    for i,pv in enumerate(pv_incr): 
        for j,r2 in enumerate(r2_incr): 
            r2_models = models[r2s >= r2]
            pv_models = models[pvs_corrected <= pv]
            count = np.intersect1d(r2_models, pv_models).shape[0]
            counts[i,j] = count  
    '''

    ## heatmap showing the number of models passing every threshold pair
    '''
    fig, ax = plt.subplots(1,1,figsize=fig_size) 
    hm = sns.heatmap(counts, **hm_kwargs, ax=ax)
    hm.set_xticklabels(r2_labels, rotation=20)
    hm.set_yticklabels(pv_labels, rotation=0)

    ax.tick_params(labelsize=label_size) 
    ax.set_title('{}\n({} total genes)'.format(reg.upper(), n_models), size=title_size, fontweight='bold') 
    ax.set_xlabel('\nr^2 values', fontsize=label_size) 
    ax.set_ylabel('p-values (FDR-corrected)\n', fontsize=label_size) 

    plt.savefig('{}/hm-small-{}.png'.format(outs_dir, reg))
    plt.close('all')
    '''

    ## compute variance across HCP subjects per gene
    with h5py.File('{}/{}.hdf5'.format(data_dir, reg), 'r') as df: 
        expr = np.array(df['pred_expr'])
        genes = np.array([g.decode("utf-8") for g in np.array(df['genes'])])
    all_var = np.var(expr, axis=1) 
    print('({:.3f}, {:.3f})'.format(all_var.min(), all_var.max()))

    ## separate complete vs incomplete models  
    ## note: this intersect call only works in numpy version ≥ 1.15.0
    incomp = np.array(unkept_models[reg])
    common, all_idx, incomp_idx = np.intersect1d(genes, \
                                                incomp, \
                                                assume_unique=True, \
                                                return_indices=True)  
    incomp_var = all_var[all_idx] 
    mask = np.ones_like(all_var, dtype=bool)
    mask[all_idx] = False 
    comp_var = all_var[mask] 

    ## histogram of HCP gene variance 
    fig, ax = plt.subplots(1,1,figsize=(15,6))
    counts, bins, _ = ax.hist([comp_var, incomp_var], \
                            bins=30, \
                            range=(0,0.65), \
                            stacked=True, \
                            linewidth=1, \
                            edgecolor='w', \
                            alpha=0.75, \
                            color=['b','r'], \
                            label=['complete models', 'missing HCP SNPs'])

    ## adjust plot aspect ratio  
    xleft, xright = ax.get_xlim()
    ybottom, ytop = ax.get_ylim()
    ax.set_aspect(abs((xright-xleft)/(ybottom-ytop)*0.25))

    ## histogram text 
    title = '{}\n({}/{} incomplete)'.format(reg.upper(), incomp_var.shape[0], n_models)
    ax.set_title(title, size=title_size, fontweight='bold')
    ax.tick_params(labelsize=label_size) 
    ax.set_xlabel('\nexpression variance across HCP', fontsize=label_size) 
    ax.set_ylabel('frequency (# genes)\n', fontsize=label_size) 
    ax.legend() 

    ## show counts above bars 
    for compc, allc, b in zip(counts[0], counts[1], bins):
        incompc = allc - compc 
        if allc > 0:
            txt = '{}\n{}'.format(int(incompc), int(compc))
            plt.gca().text(b, allc, txt) 

    plt.savefig('{}/hist-{}.png'.format(outs_dir, reg))
    plt.close('all')

