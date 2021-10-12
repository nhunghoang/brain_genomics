'''
Generate 100 sample-ID shuffles (i.e. shuffle on the genetics side) 
for the CCC simulated annealing null distribution. 

Maintain twin relationships (i.e. shuffle twin pairs with other pairs, 
non-twins with other non-twins). 

Yield a subject-ID order where twin pairs are grouped at the front and
non-twins are grouped at the end. Rearrange the sample IDs in the same 
order before separately shuffling the twin and non-twin groups.   

- Nhung, updated June 2021 

Yield a shuffling for the training group only. Index permutations will 
be with respect to the entire group (n = 890). 

TODO: Yield a shuffling for the whole group as well. 
'''

import numpy as np 
import h5py 

num_shuffles = 100 
shuffle_file = '/data1/rubinov_lab/brain_genomics/analyses_HCP/null_permutations.hdf5'

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
pair_file = '/data1/rubinov_lab/brain_genomics/analyses_HCP/subj_samp_assoc_order.hdf5'
with h5py.File(pair_file, 'r') as f: 
    all_samples = np.array(f['samples'])
    all_subjects = np.array(f['subjects'])

## consider the train/test group indices separately 
split_file = '/data1/rubinov_lab/brain_genomics/analyses_HCP/train_test_assoc_split.hdf5'
with h5py.File(split_file, 'r') as f: 
    train_idx = np.array(f['train_idx_890'])
    train_samples = np.array(f['train_samples'])

    test_idx = np.array(f['test_idx_890'])
    test_samples = np.array(f['test_samples'])

assert(np.array_equal(train_samples, all_samples[train_idx])) ## assert same order
assert(np.array_equal(test_samples, all_samples[test_idx])) ## assert same order

## split subject IDs into twin and non-twin groups  
train_twin_idx = {} ## k: family id, v: [twin1_idx, twin2_idx] 
train_non_twin_idx = [] 
for idx in train_idx: 
    samp = all_samples[idx]
    subj = all_subjects[idx] 

    try: ## is a twin 
        fid = twin_status[samp] 
        try: train_twin_idx[fid].append(idx) ## is second twin 
        except KeyError: train_twin_idx[fid] = [idx] ## is first twin 

    except KeyError: ## is not a twin 
        train_non_twin_idx.append(idx) 

test_twin_idx = {} ## k: family id, v: [twin1_idx, twin2_idx] 
test_non_twin_idx = [] 
for idx in test_idx: 
    samp = all_samples[idx]
    subj = all_subjects[idx] 

    try: ## is a twin 
        fid = twin_status[samp] 
        try: test_twin_idx[fid].append(idx) ## is second twin 
        except KeyError: test_twin_idx[fid] = [idx] ## is first twin 

    except KeyError: ## is not a twin 
        test_non_twin_idx.append(idx) 

## move singular twins to non-twin group 
train_fid_rm = [] 
for fid, twins in train_twin_idx.items(): 
    if len(twins) != 2: 
        train_fid_rm.append(fid)
for fid in train_fid_rm: 
    tidx = train_twin_idx.pop(fid) 
    train_non_twin_idx.extend(tidx) 
            
test_fid_rm = [] 
for fid, twins in test_twin_idx.items(): 
    if len(twins) != 2: 
        test_fid_rm.append(fid)
for fid in test_fid_rm: 
    tidx = test_twin_idx.pop(fid) 
    test_non_twin_idx.extend(tidx) 
            
## rearrange sample and subject orders so that twin pairs are placed first 
train_fids = list(train_twin_idx.keys())
train_new_subj_order = [] 
for fid in train_fids: 
    train_new_subj_order.extend(train_twin_idx[fid])
train_new_subj_order.extend(train_non_twin_idx) 

test_fids = list(test_twin_idx.keys())
test_new_subj_order = [] 
for fid in test_fids: 
    test_new_subj_order.extend(test_twin_idx[fid])
test_new_subj_order.extend(test_non_twin_idx) 

## save n shuffles of sample order, as well as the new subject order   
with h5py.File(shuffle_file, 'w') as f: 
    f['subj_train_idx'] = train_new_subj_order 
    f['subj_test_idx'] = test_new_subj_order 

    for n in range(num_shuffles): 
        train_new_samp_order = [] 
        
        ## shuffle twin pairs 
        np.random.shuffle(train_fids)
        for fid in train_fids: 
            train_new_samp_order.extend(train_twin_idx[fid])
        
        ## shuffle non-twins 
        np.random.shuffle(train_non_twin_idx) 
        train_new_samp_order.extend(train_non_twin_idx) 

        f['samp_train_idx_{}'.format(n)] = train_new_samp_order 

    for n in range(num_shuffles): 
        test_new_samp_order = [] 
        
        ## shuffle twin pairs 
        np.random.shuffle(test_fids)
        for fid in test_fids: 
            test_new_samp_order.extend(test_twin_idx[fid])
        
        ## shuffle non-twins 
        np.random.shuffle(test_non_twin_idx) 
        test_new_samp_order.extend(test_non_twin_idx) 

        f['samp_test_idx_{}'.format(n)] = test_new_samp_order 


