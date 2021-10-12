'''
For the multi-gene association analyses, 
we apply an 80/20 train/test split on the 
HCP cohort (n = 890) in order to assess the 
generalizability of the selected gene sets. 

This random split is applied on family groups, 
such that no family is separated into different 
groups. 

- Nhung, Aug 2021
'''

import numpy as np 
import h5py 
import sys 

## train ratio 
TR = 0.75

## output path 
SPLIT_FILE = '/data1/rubinov_lab/brain_genomics/analyses_HCP/DATA_OUTPUT/train_test_assoc_split75.hdf5' 

## input paths 
DEMO_FILE = '/data1/rubinov_lab/brain_genomics/data_HCP/subject_demographics/sampleID_race_familyID.txt' 
SAMP_FILE = '/data1/rubinov_lab/brain_genomics/analyses_HCP/DATA_OUTPUT/subj_samp_assoc_order.hdf5' 

## parse sample cohort (890 samps) 
with h5py.File(SAMP_FILE, 'r') as f: 
    cohort = np.array(f['samples'])
    cohort_size = cohort.size 
train_size = int(cohort_size*TR) 

## parse demographic info (1142 samps * [sampID, race, famID]) 
demo_all = np.loadtxt(DEMO_FILE, dtype=str, skiprows=1, delimiter='\t')
for d in demo_all: 
    if d[2][-1] == 'T': d[2] = d[2][:-1]
    d[0] = d[0][3:]
demo_dict = {int(d[0]):d for d in demo_all} ## k: sampID, v: [sampID, race, famID]

## only keep cohort info 
demo_data = np.array([demo_dict[samp] for samp in cohort])

## shuffle family groups 
fam_ids, counts = np.unique(demo_data[:,2], return_counts=True) 
shuffler = np.random.permutation(fam_ids.size)  
fam_ids = fam_ids[shuffler]
counts = counts[shuffler]

## gather training group 
train_fams = []
train_count = 0  
c = 0 
while train_count != train_size: 
    if (counts[c] + train_count) <= train_size: 
        train_fams.append(fam_ids[c])
        train_count += counts[c]  
    c += 1 
    #print('[{}]: {} people'.format(c-1, train_count))
print('c = {}\n'.format(c))

## save  
train_bool = np.isin(demo_data[:,2], train_fams) 
train_idx = np.arange(cohort_size)[train_bool]
test_idx = np.arange(cohort_size)[~train_bool]
with h5py.File(SPLIT_FILE, 'w') as f: 
    f['train_samples'] = cohort[train_idx]
    f['test_samples'] = cohort[test_idx] 
    f['train_idx_{}'.format(cohort_size)] = train_idx
    f['test_idx_{}'.format(cohort_size)] = test_idx 

## check demo percentages 
samp_race = demo_data[:,1]
train_race = samp_race[train_idx]
test_race = samp_race[test_idx] 

print('TRAIN')
us, cs = np.unique(train_race, return_counts=True)
for u,c in zip(us,cs): 
    print('{}: {:>.3f}'.format(u, c/train_idx.size))
print('\nTEST')
us, cs = np.unique(test_race, return_counts=True)
for u,c in zip(us,cs): 
    print('{}: {:>.3f}'.format(u, c/test_idx.size))
