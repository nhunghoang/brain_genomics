'''
Analyze the prediction performances of the PrediXcan 
gene models, in terms of r^2 and p-values. 

Use analysis to decide quality thresholds for the initial 
set of genes per region (prior to phenotype correlation). 

- Nhung, June 2021 
'''

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os 

## input and output paths 
data_dir = '/data1/rubinov_lab/brain_genomics/data_HCP/expression'
outs_dir = '/data1/rubinov_lab/brain_genomics/models_PrediXcan_v8/plots_r2_pval'

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
#r2_incr = np.arange(0.1, 0.71, 0.1) 
r2_incr = np.arange(0, 0.16, 0.025) 
pv_incr = np.arange(0.01, 0.051, 0.01)

## heatmap details 
r2_labels = ['≥' + str(r2) for r2 in r2_incr]
pv_labels = ['≤' + str(pv) for pv in pv_incr]
hm_kwargs = {'cmap':'plasma', 'xticklabels':r2_labels, 'yticklabels':pv_labels, \
            'fmt':'4d', 'square': True, 'annot':True, 'annot_kws':{'size':16}}
fig_size = (20,12)
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

    ## scatter plot of r^2 vs p-value, with points representing models 
    #plt.scatter(r2s, pvs, s=1)
    #plt.title('{}\n({} total genes)'.format(reg.upper(), n_models), size=title_size-5, fontweight='bold') 
    #plt.xlabel('r^2 values', fontsize=label_size-5) 
    #plt.ylabel('p-values', fontsize=label_size-5) 
    #plt.savefig('{}/sc-{}.png'.format(outs_dir, reg))
    #plt.close('all')

    ## for this specific threshold pair, report yield and the number of models that are missing HCP SNPs
    r2_models = models[r2s >= 0.05]
    pv_models = models[pvs <= 0.01]
    common_models = np.intersect1d(r2_models, pv_models)
    incomplete_models = np.intersect1d(common_models, np.array(unkept_models[reg]))

    count = common_models.shape[0]
    inc_count = incomplete_models.shape[0] 
    pyield = (float(count) / n_models) * 100
    print('{:>25s}: {} / {} = {:.2f}% yield, includes {} incomplete models'.format(reg, count, n_models, pyield, inc_count))

    ## count the number of models that meet every pair of (r^2,p) thresholds 
    for i,pv in enumerate(pv_incr): 
        for j,r2 in enumerate(r2_incr): 
            r2_models = models[r2s >= r2]
            pv_models = models[pvs <= pv]
            count = np.intersect1d(r2_models, pv_models).shape[0]
            counts[i,j] = count  

    ## heatmap showing the number of models passing every threshold pair
    fig, ax = plt.subplots(1,1,figsize=fig_size) 
    hm = sns.heatmap(counts, **hm_kwargs, ax=ax)
    hm.set_yticklabels(pv_labels, rotation=0)

    ax.tick_params(labelsize=label_size) 
    ax.set_title('{}\n({} total genes)'.format(reg.upper(), n_models), size=title_size, fontweight='bold') 
    ax.set_xlabel('r^2 values', fontsize=label_size) 
    ax.set_ylabel('p-values', fontsize=label_size) 

    plt.savefig('{}/hm-{}.png'.format(outs_dir, reg))
    plt.close('all')
