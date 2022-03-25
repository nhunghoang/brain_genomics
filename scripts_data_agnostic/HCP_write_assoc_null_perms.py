'''
Generate 100 sample-ID shuffles (i.e. shuffle on the genetics side) 
for the CCC simulated annealing null distribution. 

Maintain twin relationships (i.e. shuffle twin pairs with other pairs, 
non-twins with other non-twins). 

Yield a subject-ID order where twin pairs are grouped at the front and
non-twins are grouped at the end. Rearrange the sample IDs in the same 
order before separately shuffling the twin and non-twin groups.   

- Nhung, updated June 2021 

[11.15.21] Update: Generate 10k shuffles (needed for single-gene assocs). 
For the multi-gene assocs, randomly sample 100 perms from this 10k set. 
This is happening over the entire cohort (no split). 
THIS ONE FOR THE NULLS HERE. 
'''

import numpy as np 
import h5py 

#####

## DATA_OUTPUT/null_permutations_10k.hdf5: 10k perms used to compute p-values for observed single-gene assocs 
## assoc/nulls_for_nulls.hdf5: 10k perms used to compute p-values for null single-gene assocs 
## assoc/nulls_for_multi.hdf5: 100 perms used to compute p-values for observed multi-gene assocs 

#####

num_shuffles = 900 
shuffle_file = '/data1/rubinov_lab/brain_genomics/analyses_HCP/assoc/nulls_of_multi2.hdf5'

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

## save n shuffles of sample order, as well as the new subject order   
with h5py.File(shuffle_file, 'w') as f: 
    f['subj_idx'] = new_subj_order 

    for n in range(num_shuffles): 
        new_samp_order = [] 
        
        ## shuffle twin pairs 
        np.random.shuffle(fids)
        for fid in fids: 
            new_samp_order.extend(twin_idx[fid])
        
        ## shuffle non-twins 
        np.random.shuffle(non_twin_idx) 
        new_samp_order.extend(non_twin_idx) 

        #f['samp_idx_{}'.format(n)] = new_samp_order 
        f['samp_idx_{}'.format(n+100)] = new_samp_order ## only for multi2 

        if (n+1)%1000 == 0: 
            perc = int(((n+1)/num_shuffles)*100)
            print('{:d} %'.format(perc))
