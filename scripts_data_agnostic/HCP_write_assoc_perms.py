'''
Generate 1000 shuffles for the single-gene associations. 
Generate another 10k shuffles, to be used as the permutations 
for these 1000 shuffles.  

Maintain twin relationships (i.e. shuffle twin pairs with other pairs, 
non-twins with other non-twins). 

Yield a subject-ID order where twin pairs are grouped at the front and
non-twins are grouped at the end. Rearrange the sample IDs in the same 
order before separately shuffling the twin and non-twin groups.   

- Nhung, updated April 2022  
'''

import numpy as np 
import h5py 

shuffle_file = '/data1/rubinov_lab/brain_genomics/analyses_HCP/DATA_OUTPUT/null_permutations.hdf5'

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

## get corresponding sample & subject IDs 
pair_file = '/data1/rubinov_lab/brain_genomics/analyses_HCP/DATA_OUTPUT/subj_samp_assoc_order.hdf5'
with h5py.File(pair_file, 'r') as f: 
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
fids = list(twin_idx.keys())
new_subj_order = [] 
for fid in fids: 
    new_subj_order.extend(twin_idx[fid])
new_subj_order.extend(non_twin_idx) 

## function: generate a permutation 
def get_permutation(): 
    perm = [] 
    
    ## shuffle twin pairs 
    np.random.shuffle(fids)
    for fid in fids: 
        perm.extend(twin_idx[fid])
    
    ## shuffle non-twins 
    np.random.shuffle(non_twin_idx) 
    perm.extend(non_twin_idx) 

    return perm 

## save n shuffles of sample order, as well as the new subject order   
num_subjects = len(new_subj_order) 
nulls = np.zeros((1000, num_subjects), dtype=int) 
nulls_of_nulls = np.zeros((1000, 1000, num_subjects), dtype=int) 

for i in range(1000): 
    nulls[i] = get_permutation() 
    for ii in range(1000): 
        nulls_of_nulls[i][ii] = get_permutation() 

    if (i+1)%100 == 0: 
        perc = int(((i+1)/1000)*100)
        print('{:d} %'.format(perc))

with h5py.File(shuffle_file, 'w') as f: 
    f['subj_idx'] = new_subj_order 
    f['samp_idx'] = nulls 
    f['null_idx'] = nulls_of_nulls 


