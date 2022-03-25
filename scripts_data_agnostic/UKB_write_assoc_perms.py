'''
Generate 10k shuffles (for the single-gene assocs). 
For the multi-gene assocs, randomly sample 100 perms from this 10k set. 

Note: unlike HCP, there are no twin relationships to constrain in the UKB data. 

- Tim, Mar 15 2022  
'''


import numpy as np 
import h5py 
import sys 

dset = 'UKB'

num_shuffles = 10000 
shuffle_file = '/data1/rubinov_lab/brain_genomics/analyses_{}/DATA_OUTPUT/null_permutations_10k.hdf5'.format(dset)

sID = np.loadtxt("/data1/rubinov_lab/brain_genomics/data_UKB/ordered_subject_list.txt", dtype=int)
new_samp_order = np.array(list(range(sID.shape[0])))

## save n shuffles of sample order, as well as the new subject order   
with h5py.File(shuffle_file, 'w') as f: 
    f['subj_idx'] = new_samp_order

    for n in range(num_shuffles): 
        np.random.shuffle(new_samp_order)

        f['samp_idx_{}'.format(n)] = new_samp_order.tolist()

        if (n+1)%1000 == 0: 
            perc = int(((n+1)/num_shuffles)*100)
            print('{:d} %'.format(perc))
