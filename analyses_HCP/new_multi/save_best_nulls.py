'''
Each phenotype has a file. 
Keys: regions 
Data: (100 perms, best train & corresponding test) --> (100,2) 
'''

import numpy as np 
import h5py 
import os 

mpath = '/data1/rubinov_lab/brain_genomics/analyses_HCP/new_multi'
out_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/new_multi/best_nulls'

phens = ['alff', 'regional_homogeneity', 'gm_volume']

for phen in phens: 
    pout = '{}/{}.hdf5'.format(out_path, phen) 
    pf = h5py.File(pout, 'w') 

    pdir = '{}/null_results_{}'.format(mpath, phen)
    for reg in os.listdir(pdir): 

        rtrain = np.zeros((100,10)) ## 100 shuffles, 10 runs
        rtestt = np.zeros((100,10)) ## 100 shuffles, 10 runs

        rdir = '{}/{}'.format(pdir, reg) 
        for log in os.listdir(rdir): 
            with open('{}/{}'.format(rdir, log), 'r') as f: 
                lines = f.readlines()
            [sh, rn] = log.split('.')[0].split('_')[1:]
            sh = int(sh); rn = int(rn)
            ltrain = float(lines[1].split('|')[1].split(' ')[1]) 
            ltestt = float(lines[2].split('|')[1].split(' ')[1]) 
            rtrain[sh][rn] = ltrain 
            rtestt[sh][rn] = ltestt

        best_nulls = np.zeros((100,2)) ## best train rho per shuffle & corresponding test rho 
        idx = np.argmax(rtrain, axis=1)
        for j,i in enumerate(idx): 
            best_nulls[j][0] = rtrain[j][i]
            best_nulls[j][1] = rtestt[j][i]

        pf[reg] = best_nulls 
    pf.close()
