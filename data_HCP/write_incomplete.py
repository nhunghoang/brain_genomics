'''
'''

import sys
import h5py 
import numpy as np 
from scipy.io import savemat

data_dir = sys.argv[1]

with h5py.File('/data1/rubinov_lab/brain_genomics/analyses_HCP/multi_gene_assoc/ccc_constants.hdf5', 'r') as f: 
    old_subjs = np.array(f['IDs'], dtype=str)[1]

with h5py.File('/data1/rubinov_lab/brain_genomics/data_HCP/preprocessed_schmel/schmel_nofilt/timeseries_order.hdf5', 'r') as f: 
    old_scans = {s:np.array(f[s]) for s in old_subjs}

with h5py.File('/data1/rubinov_lab/brain_genomics/data_HCP/{}/timeseries_order.hdf5'.format(data_dir), 'r') as f: 
    new_subjs = np.array(f['subjects'], dtype=str) 
    new_scans = {s:np.array(f[s]) for s in new_subjs}
    

inc_scans = {} 
non_scans = [] 
for old_subj in old_subjs: 
    try:
        read_scans = new_scans[old_subj]
        unread_scans = np.setdiff1d(old_scans[old_subj], read_scans)
        if unread_scans.size:
            inc_scans['s' + old_subj] = unread_scans + 1  
    except KeyError: 
        #non_scans['s' + old_subj] = old_scans[old_subj] + 1
        non_scans.append(old_subj)

for a, b in inc_scans.items():
    print('{}: {}'.format(a,b))

inc_file = '/data1/rubinov_lab/brain_genomics/data_HCP/{}/incomplete_scans.mat'.format(data_dir)
savemat(inc_file, inc_scans)
print('{} subjects have incomplete scans'.format(len(inc_scans)))

#non_file = '/data1/rubinov_lab/brain_genomics/data_HCP_hoacer/nonexist_scans.mat'
#savemat(non_file, non_scans)
#print('{} subjects have non-existent scans'.format(len(non_scans)))
non_file = '/data1/rubinov_lab/brain_genomics/data_HCP/{}/nonexist_scans.txt'.format(data_dir)
with open(non_file, 'w') as f: 
    for s in non_scans:
        f.write('{}\n'.format(s))
print('{} subjects have non-existent scans'.format(len(non_scans)))
