'''
Write neural phenotypes to hdf5 files. 

- Nhung, May 2021 
'''

import numpy as np 
from scipy.stats import pearsonr
import h5py 
import os 

## (RUN ONCE) Save order of available subject timeseries  
## hdf5 keys: subject-array, subject-to-timeseries-dictionary  
def save_timeseries_order(mat_dir, mat_files, ts_file):
    subjects = [] 
    timeseries = {} 
    for mat_file in mat_files: 
        with h5py.File('{}/{}'.format(mat_dir, mat_file), 'r') as f: 
            objs = f.get('timeseries') 
            ts = [] 
            for t in range(4): 
                obj = objs[0][t] 
                name = h5py.h5r.get_name(obj, f.id) 
                ts_shape = np.array(f[name]).shape
                if ts_shape == (1200,119): 
                    ts.append(t) 
            if len(ts) != 0:
                subj = mat_file.split('.')[0]
                subjects.append(subj)
                timeseries[subj] = ts 
    with h5py.File(ts_file, 'w') as f1:
        f1['subjects'] = np.array(subjects, dtype=int)
        for subj in subjects: 
            f1[subj] = timeseries[subj]

## Load order of available subject timeseries 
## hdf5 keys: subject-array, subject-to-timeseries-dictionary  
def load_timeseries_order(ts_file): 
    with h5py.File(ts_file, 'r') as f: 
        subjects = [str(s) for s in np.array(f['subjects'])] 
        mapping = {s:np.array(f[s]) for s in subjects}
    return subjects, mapping 

## Compute and write coactivity matrices 
## Note: normalize scan and zero out self-connections 
## hdf5 key: [subject-id] // val: (nscans,119,119) matrix 
def write_coactivity_matrices(mat_dir, subjects, ts_map, coact_file): 
    cf = h5py.File(coact_file, 'w') 
    for subj in subjects: 
        coact_all = [] 
        scan_nums = ts_map[subj]
        with h5py.File('{}/{}.mat'.format(mat_dir, subj), 'r') as f: 
            objs = f.get('timeseries') 
            for sn in scan_nums: 
                obj = objs[0][sn] 
                name = h5py.h5r.get_name(obj, f.id) 
                ts = np.array(f[name]).T ## (119,1200) 
                ts_mean = np.mean(ts, axis=1)[:,None] 
                ts_norm = (ts-ts_mean)/ts_mean
                coact = np.corrcoef(ts_norm) 
                np.fill_diagonal(coact, 0) 
                coact_all.append(coact) 
        cf[subj] = np.array(coact_all)
    cf.close()

## Compute and write connectivity mean and variance for regions of interest 
## hdf5 keys: regions // val: value in subject-scan order 
def write_conn_mean_var(subjs, coact_file, reg_idx): 

    with h5py.File(coact_file, 'r') as f: 
        all_coacts = {s:np.array(f[s]) for s in subjs}

    nsubjs = len(subjs) 
    all_means = {r:np.zeros(nsubjs, dtype=float) for r in reg_idx} 
    all_varis = {r:np.zeros(nsubjs, dtype=float) for r in reg_idx} 

    for s,subj in enumerate(subjs): 
        coacts = all_coacts[subj]
        ## compute mean/var for each scan across regions THEN average scans together 
        means = np.mean(np.mean(coacts, axis=1), axis=0)
        varis = np.mean(np.var(coacts, axis=1), axis=0)
        for reg,idx in reg_idx.items(): 
            all_means[reg][s] = means[idx]
            all_varis[reg][s] = varis[idx]     

    with h5py.File('/data1/rubinov_lab/brain_genomics/data_HCP/phenotypes/connectivity_mean.hdf5', 'w') as m: 
        for reg,arr in all_means.items(): 
            m[reg] = arr  
    with h5py.File('/data1/rubinov_lab/brain_genomics/data_HCP/phenotypes/connectivity_variance.hdf5', 'w') as v: 
        for reg,arr in all_varis.items(): 
            v[reg] = arr  


## Main 
def main(): 

    region_indices = {'AMYG-lh':109, 'AMYG-rh':101, 'HIP-lh':108, 'HIP-rh':100, 'PUT-lh':114, 'PUT-rh':106, \
                    'NAc-lh':112, 'NAc-rh':104, 'CAU-lh':115, 'CAU-rh':107, 'CER-lh':116, 'CER-rh':117, 'CER-vermis':118}

    mat_dir = '/data1/rubinov_lab/brain_genomics/data_HCP_neuro' 
    mat_files = os.listdir(mat_dir) 

    ts_file = '/data1/rubinov_lab/brain_genomics/data_HCP/timeseries_order.hdf5'
    #save_timeseries_order(mat_dir, mat_files, ts_file) 
    subjs, subj_ts_map = load_timeseries_order(ts_file) 

    coact_file = '/data1/rubinov_lab/brain_genomics/data_HCP/phenotypes/coactivity_matrices.hdf5'
    #write_coactivity_matrices(mat_dir, subjs, subj_ts_map, coact_file) 

    write_conn_mean_var(subjs, coact_file, region_indices) 

main() 


