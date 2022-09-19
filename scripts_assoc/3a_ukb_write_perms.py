'''
For UKB - 
.Generate N perms (for reproducible gene set analyses). 
.Permutations are entirely random (b/c no family structures). 

- Nhung, updated June 2022  
'''

import numpy as np 
import h5py 

nperms = int(1e5)  

## paths 
file_sub = '/data1/rubinov_lab/brain_genomics/data_UKB/UKB_subjs_941.txt'
file_out = '/data1/rubinov_lab/brain_genomics/analyses_UKB/DATA_OUTPUT/null_permutations_latest_100k.hdf5'

## get subject indices 
with open(file_sub, 'r') as f: 
    nsubjs = len(f.readlines())

## set sample order and permutations of subject order 
samp_order = np.arange(nsubjs, dtype=int)  
subj_order = np.zeros((nperms, nsubjs), dtype=int) 
for i in range(nperms): 
    subj_order[i] = np.random.permutation(samp_order) 

## save permutation indices   
with h5py.File(file_out, 'w') as f: 
    f['samp_idx'] = samp_order 
    f['subj_idx'] = subj_order 

