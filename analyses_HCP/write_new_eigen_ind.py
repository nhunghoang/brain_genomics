'''
Script for writing a new .ind file (needed for EIGENSTRAT) 
such that the label column is formatted as [train/test][race]. 

This will allow us to run EIGENSTRAT on the given train/test 
groups separately (using the -w flag) without having to regenerate 
the .snp and .geno files. 

doc: https://github.com/DReichLab/EIG/tree/master/EIGENSTRAT 

- Nhung, Oct 2021 
'''

import numpy as np 
import h5py 

ind_file_0 = '/data1/rubinov_lab/brain_genomics/analyses_HCP/DATA_OUTPUT/eigendata/hcp_cohort.ind' 

split_file = '/data1/rubinov_lab/brain_genomics/analyses_HCP/DATA_OUTPUT/train_test_assoc_split75.hdf5'
split_outf = '/data1/rubinov_lab/brain_genomics/analyses_HCP/DATA_OUTPUT/eigendata/hcp_split75.ind'

with h5py.File(split_file, 'r') as f: 
    test_idx = np.array(f['test_idx_890'])

ind_all = np.loadtxt(ind_file_0, delimiter='\t', dtype=str) 
for i in range(ind_all.shape[0]): 
    if np.isin(i, test_idx): 
        ind_all[i][2] = 'test-' + ind_all[i][2] 
    else: 
        ind_all[i][2] = 'train-' + ind_all[i][2] 

np.savetxt(split_outf, ind_all, delimiter='\t', fmt='%s')

## poplist.txt for -w flag of EIGENSTRAT 
train_pop = '/data1/rubinov_lab/brain_genomics/analyses_HCP/DATA_OUTPUT/eigendata/train_pop75.txt'
testt_pop = '/data1/rubinov_lab/brain_genomics/analyses_HCP/DATA_OUTPUT/eigendata/test_pop75.txt'

pops = ['Asian/Nat. Hawaiian/Othr Pacifi', 'Black or African Am.', 'More than one', 'Unknown or Not Reported', 'White', 'Am. Indian/Alaskan Nat.']

with open(train_pop, 'w') as f: 
    for pop in pops: 
        f.write('train-{}\n'.format(pop))

with open(testt_pop, 'w') as f: 
    for pop in pops[:-1]: ## last pop not rep'd in test group 
        f.write('test-{}\n'.format(pop))
