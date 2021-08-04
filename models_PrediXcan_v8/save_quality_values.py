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
r2_path = '/data1/rubinov_lab/brain_genomics/models_PrediXcan_v8/filtered_r2_distribution.hdf5'
pv_path = '/data1/rubinov_lab/brain_genomics/models_PrediXcan_v8/filtered_pv_distribution.hdf5'

if os.path.exists(r2_path): os.remove(r2_path) 
if os.path.exists(pv_path): os.remove(pv_path) 

all_r2_path = '/data1/rubinov_lab/brain_genomics/models_PrediXcan_v8/all_r2_distribution.hdf5'
all_pv_path = '/data1/rubinov_lab/brain_genomics/models_PrediXcan_v8/all_pv_distribution.hdf5'

## loop through brain regions 
for log in os.listdir(predix_dir): 
    if log[-4:] != '.log': continue 
    region = log.split('.')[0] 
    if region == 'cerebellum': continue
    if region == 'cortex': continue 
    if region == 'spinal-cord': continue 

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

    r2_pass = r2s[both] 
    pv_pass = pvs[both] 
    
    ## save region-specific file of gene models kept 
    with h5py.File(r2_path, 'a') as fr2:
        fr2[region] = r2_pass 
    with h5py.File(pv_path, 'a') as fpv:
        fpv[region] = pv_pass 

    with h5py.File(all_r2_path, 'a') as f0: 
        f0[region] = r2s
    with h5py.File(all_pv_path, 'a') as f1: 
        f1[region] = pvs
        
    print('{:22}: {} / {} genes kept ({} {})'\
        .format(region, genes_keep.shape[0], genes.shape[0], r2_pass.size, pv_pass.size))
