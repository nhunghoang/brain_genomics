'''
Parse the PrediXcan expression imputation logs, and filter gene models 
that do not pass the given r^2 and p-value thresholds. For the models that 
do pass, make a note of their SNP-completeness (proportion of weights missing).  

- Nhung, June 2021 

TODO: 
- decide if it's worth to keep file of unkept genes 
- record missing weight proportion next to kept genes files 
'''

import os 
import sys
import numpy as np 
import h5py 

r2_thr = float(sys.argv[1])
pv_thr = float(sys.argv[2])

predix_dir = '/data1/rubinov_lab/brain_genomics/data_HCP/expression'
filter_dir = '/data1/rubinov_lab/brain_genomics/data_HCP/expression/filtered_quality_r{}_p{}'.format(r2_thr, pv_thr)
gmodel_dir = '/data1/rubinov_lab/brain_genomics/data_HCP/predixcan_v8/gene_models_quality_r{}_p{}'.format(r2_thr, pv_thr)

if not os.path.exists(filter_dir): 
    os.mkdir(filter_dir)
    os.mkdir(gmodel_dir)

#no_keep_file = '/data1/rubinov_lab/brain_genomics/data_HCP/expression/not_kept_model_quality.txt' 
#no_keep = open(no_keep_file, 'w') 

## loop through brain regions 
for log in os.listdir(predix_dir): 
    if log[-4:] != '.log': continue 
    region = log.split('.')[0] 

    ## parse PrediXcan imputation log 
    with open('{}/{}'.format(predix_dir, log), 'r') as f: 
        f.readline() ## header 
        lines = f.readlines() 
    data = np.array([line.strip().split('\t') for line in lines])
    genes = data[:,0] 
    gidxs = np.arange(genes.shape[0])
    r2s = data[:,4].astype(float)
    pvs = data[:,5].astype(float)

    ## filter gene models based on r^2 and p-value thresholds 
    r2_filter = gidxs[r2s >= r2_thr]  
    pv_filter = gidxs[pvs <= pv_thr] 
    both = np.intersect1d(r2_filter, pv_filter, assume_unique=True) 
    genes_keep = genes[both]  
    
    ## save region-specific file of gene models kept 
    ## TODO: record missing weight proportion next to each gene 
    with open('{}/{}.txt'.format(gmodel_dir, region), 'w') as f1: 
        for gene in genes_keep: 
            f1.write(gene + '\n') 

    ## write new expression hdf5 files with just the valid genes 
    expr_data = {} 
    with h5py.File('{}/{}.hdf5'.format(predix_dir, region), 'r') as f2: 
        expr_data = {key:np.array(f2[key]) for key in f2.keys()}
    with h5py.File('{}/{}.hdf5'.format(filter_dir, region), 'w') as f3: 
        f3['genes'] = expr_data['genes'][both]
        f3['pred_expr'] = expr_data['pred_expr'][both] 
        f3['samples'] = expr_data['samples']   

    print('{:22}: {} / {} genes kept'.format(region, genes_keep.shape[0], genes.shape[0]))
