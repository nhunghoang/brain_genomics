'''
For HCP - 
.Sort subjects into twin and non-twin groups (for associations).
.Generate N perms (for reproducible gene set analyses). 
.Maintain twin relationships (i.e. shuffle twin pairs 
 with other pairs, non-twins with other non-twins). 

- Nhung, updated June 2022  
'''

import numpy as np 
import h5py 

nperms = int(1e5)  

## paths 
file_dem = '/data1/rubinov_lab/brain_genomics/data_HCP/subject_demographics/sampleID_race_familyID.txt'
file_ids = '/data1/rubinov_lab/brain_genomics/analyses_HCP/DATA_OUTPUT/subj_samp_assoc_order.hdf5'
file_out = '/data1/rubinov_lab/brain_genomics/analyses_HCP/DATA_OUTPUT/null_permutations_latest_100k.hdf5'

## parse demographics for twin details 
twin_status = {} ## k: twin sample id (w/o 'MH0'), v: family id 
with open(file_dem, 'r') as f: 
    f.readline() 
    for line in f.readlines(): 
        [sampID, race, famID] = line.strip().split('\t') 
        if famID[-1] == 'T': 
            sID = int(sampID[3:])
            twin_status[sID] = famID 

## get corresponding sample & subject IDs 
with h5py.File(file_ids, 'r') as f: 
    all_samples = np.array(f['samples'])
    all_subjects = np.array(f['subjects'])

## split subject IDs into twin and non-twin groups  
twin_idx = {} ## k: family id, v: [twin1_idx, twin2_idx] 
non_twin_idx = [] 
for idx, sample in enumerate(all_samples):

    try: ## is a twin 
        fid = twin_status[sample] 
        try: twin_idx[fid].append(idx) ## is second twin 
        except KeyError: twin_idx[fid] = [idx] ## is first twin 

    except KeyError: ## is not a twin 
        non_twin_idx.append(idx) 

## move singular twins to non-twin group 
fid_rm = [] 
for fid, twins in twin_idx.items(): 
    if len(twins) != 2: 
        fid_rm.append(fid)
for fid in fid_rm: 
    tidx = twin_idx.pop(fid) 
    non_twin_idx.extend(tidx) 
            
## rearrange sample and subject orders so that twin pairs are placed first 
twin_indices = np.array([tidx for fid, tidx in twin_idx.items()], dtype=int) ## (twin pairs, 2) 
nontwin_indices = np.array(non_twin_idx, dtype=int) 
new_subj_order = np.concatenate((twin_indices.flatten(), nontwin_indices)) 

## function: generate a permutation 
def get_permutation(): 
    twin_perm = np.random.permutation(twin_indices) 
    nontwin_perm = np.random.permutation(nontwin_indices)
    new_perm = np.concatenate((twin_perm.flatten(), nontwin_perm))
    return new_perm

## save n shuffles of sample order, as well as the new subject order   
num_subjects = len(new_subj_order) 
nulls = np.zeros((nperms, num_subjects), dtype=int) 
#nulls_of_nulls = np.zeros((nperms, nperms, num_subjects), dtype=int) 

for i in range(nperms): 
    nulls[i] = get_permutation() 
    #for ii in range(nperms): 
    #    nulls_of_nulls[i][ii] = get_permutation() 

    if (i+1)%100 == 0: 
        perc = int(((i+1)/nperms)*100)
        print('{:d} %'.format(perc))

## NOTE: apply permutations to phenotypes now, 
## not expression matrices 
## so here, I reassign the indices
with h5py.File(file_out, 'w') as f: 
    f['samp_idx'] = new_subj_order 
    f['subj_idx'] = nulls 
    #f['null_idx'] = nulls_of_nulls 
    f['twin_idx'] = twin_indices 
    f['nontwin_idx'] = nontwin_indices 


