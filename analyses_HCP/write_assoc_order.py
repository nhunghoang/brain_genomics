'''
For the multi- and single-gene association analyses, 
we only consider individuals with both subject and 
sample IDs (n = 891). This script saves two ordered 
sets of indices: 

- subject_idx: ordered indices relative to phenotype data (n = 939)
- sample_idx: corresponding indices, but relative to expression data (n = 1142)  

These index sets should be used to reorder data in the 
association scripts. 

Nhung, July 2021
'''

import numpy as np 
import h5py

## map sample IDs to subject IDs
id_file = '/data1/rubinov_lab/brain_genomics/data_HCP/neuro_genetic_ids.txt'
samp_to_subj = {}
with open(id_file, 'r') as f:
    f.readline()
    for line in f.readlines():
        [subj, samp] = line.strip().split('\t')
        samp_to_subj[samp] = subj

## read subject order in which phenotypes are saved in (n = 939)   
subj_file = '/data1/rubinov_lab/brain_genomics/data_HCP/timeseries_order.hdf5'
with h5py.File(subj_file, 'r') as f:
    subj_order_all = np.array([str(s) for s in np.array(f['subjects'])])

## read sample order in which expression are saved in (n = 1142)  
samp_file = '/data1/rubinov_lab/brain_genomics/data_HCP/predixcan_v8/sample_ids.txt'
with open(samp_file, 'r') as f:
    samp_order_all = np.array([s.split('\t')[0] for s in f.readlines()])

## only keep individuals with valid subject & sample IDs (n = 891) 
subjects = []; samples = [] 
subj_idx = []; samp_idx = []
for i,samp in enumerate(samp_order_all): 
    try:
        subj_of_samp = samp_to_subj[samp_order_all[i]]
        idx_of_subj = np.argwhere(subj_order_all == subj_of_samp)[0][0]

        subjects.append(subj_of_samp)
        samples.append(samp[3:]) ## ignore 'MH0' for hdf5 write  

        samp_idx.append(i)
        subj_idx.append(idx_of_subj)
    except:
        continue

## save the index sets 
idx_file = '/data1/rubinov_lab/brain_genomics/analyses_HCP/subj_samp_assoc_order.hdf5'
with h5py.File(idx_file, 'w') as f: 
    f['subjects'] = np.array(subjects, dtype=int) 
    f['samples'] = np.array(samples, dtype=int)
    f['subject_idx_939'] = subj_idx 
    f['sample_idx_1142'] = samp_idx  

