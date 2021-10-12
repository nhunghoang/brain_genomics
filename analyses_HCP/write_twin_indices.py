'''
Script for saving indices of twins with shape (# twin pairs, 2) 
and indices of non-twins with shape (# non-twins). 
- This is for the HCP cohort (n = 890)
- Very useful for the 10k permutations for the single-gene associations. 

[10.06.21] For now, doing this on whole cohort (not just training cohort).  
[10.07.21] Training cohort option has been implemented. 

- Nhung, Oct 2021 
'''

import numpy as np 
import h5py 

train = True 

out_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/twin_indices.hdf5'
if train:
    out_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/train_twin_indices.hdf5'

## get sample IDs in cohort order 
pairs_file = '/data1/rubinov_lab/brain_genomics/analyses_HCP/subj_samp_assoc_order.hdf5'
split_file = '/data1/rubinov_lab/brain_genomics/analyses_HCP/train_test_assoc_split.hdf5'
if train:
    with h5py.File(split_file, 'r') as f: 
        all_samples = np.array(f['train_samples'])
else:
    with h5py.File(pairs_file, 'r') as f:
        all_samples = np.array(f['samples'])

## parse demographics for twin details
demo_file = '/data1/rubinov_lab/brain_genomics/data_HCP/subject_demographics/sampleID_race_familyID.txt'
twin_status = {} ## k: twin sample id (w/o 'MH0'), v: family id
with open(demo_file, 'r') as f:
    f.readline()
    for line in f.readlines():
        [sampID, race, famID] = line.strip().split('\t')
        if famID[-1] == 'T':
            sID = int(sampID[3:])
            twin_status[sID] = famID

## split sample IDs into twin and non-twin groups
twin_pairs = {} ## k: family id, v: [twin1_idx, twin2_idx]
non_twin_idx = []
for idx,samp in enumerate(all_samples): 

    try: ## is a twin
        fid = twin_status[samp]
        try: twin_pairs[fid].append(idx) ## is second twin
        except KeyError: twin_pairs[fid] = [idx] ## is first twin

    except KeyError: ## is not a twin
        non_twin_idx.append(idx)
        
## move singular twins to non-twin group 
## while creating the (# twin pairs, 2) array 
twin_idx = [] 
fid_rm = []

for fid, twins in twin_pairs.items():
    if len(twins) != 2: fid_rm.append(fid) ## singular twin 
    else: twin_idx.append(np.array(twins)) ## twin pair 
         
for fid in fid_rm:
    tidx = twin_pairs.pop(fid)
    non_twin_idx.extend(tidx)

## save the twin and non-twin indices 
with h5py.File(out_path, 'w') as f: 
    f['twins_idx'] = np.array(twin_idx)
    f['non_twins_idx'] = np.array(non_twin_idx) 

## print counts 
print('# twin pairs: {}'.format(len(twin_idx)))
print(' # non-twins: {}'.format(len(non_twin_idx)))
