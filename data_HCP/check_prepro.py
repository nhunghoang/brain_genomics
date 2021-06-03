'''
Script for checking the rs-fMRI preprocessing, including 
the number of processed timeseries compared to the number 
of available scans. 

- Nhung, May 2021 
'''

import numpy as np 
import h5py 
import os 

## count the number of available scans (voxel form) per subject 
rs = ['rfMRI_REST1_LR', 'rfMRI_REST1_RL', 'rfMRI_REST2_LR', 'rfMRI_REST2_RL'] 
subjs = os.listdir('/data1/datasets/hcp/')
n_exists = {s:0 for s in subjs} # k: subj, v: num exists 
for s in subjs: 
    for state in rs: 
        dpath = '/data1/datasets/hcp/{}/MNINonLinear/Results/{}/{}_hp2000_clean.nii.gz'.format(s, state, state)
        if os.path.exists(dpath): n_exists[s] += 1 

## count the number of timeseries (regions x timepoints) per subject 
mat_dir = '/data1/rubinov_lab/brain_genomics/data_HCP/timeseries_seq' 
#mat_dir = '/data1/rubinov_lab/brain_genomics/data_HCP/timeseries' 
#mat_dir = '/data1/rubinov_lab/brain_genomics/data_HCP/ts_missed2' 
mat_files = os.listdir(mat_dir) 
counts = [] # best ts shape and num of available scans 

has_at_least_one_ts = [] 
has_scans_but_no_ts = [] 

for i,mat_file in enumerate(mat_files): 
    subj = mat_file.split('.')[0]
    best_ts = np.zeros((0,0))
    with h5py.File('{}/{}'.format(mat_dir,mat_file), 'r') as f: 

        ## timeseries 
        objs = f.get('timeseries')
        for t in range(4):
            try: 
                obj = objs[0][t]
                name = h5py.h5r.get_name(obj, f.id)
                curr_ts = np.array(f[name])
                if curr_ts.size > best_ts.size: 
                    best_ts = curr_ts
            except: 
                pass  
        info = 'best ts shape: {:12} | num scans: {:1}'.format(str(best_ts.shape), n_exists[subj])
        counts.append(info) 

        if (n_exists[subj] > 0) and (best_ts.size > 2): 
            has_at_least_one_ts.append(subj) 
        
        if (n_exists[subj] > 0) and (best_ts.size == 2): 
            has_scans_but_no_ts.append(subj) 

        #if best_ts.size <= 2: 
        #    print(mat_file.split('.')[0])
    
        ## regressors 
        '''
        rg = f.get('regressors')
        for j in range(4):
            try:
                mvt = np.array(f[rg[0][j]]['Movt_RRMS'])
            except:
                pass 
        '''

## summary 
counts = np.array(counts) 
n_mats = len(mat_files) 
shapes, shape_counts = np.unique(counts, return_counts=True) 
for s, c in zip(shapes, shape_counts): 
    print('[ {:4} ] {}'.format(c, s))

print('{} subjects have scans but no timeseries'.format(len(has_scans_but_no_ts))) 
print('{}/{} subjects have at least one timeseries'.format(len(has_at_least_one_ts), n_mats))

## write subjects with at least one scan but no timeseries 
with open('subjs_missing_ts.txt', 'w') as f: 
    for s in has_scans_but_no_ts: 
        f.write(s + '\n') 
