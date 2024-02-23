'''
Save random subsample indices for the UKB/HCP data. 

- Nhung, Feb 2024
'''

import numpy as np 
import h5py  

main_path = '/data1/rubinov_lab/brain_genomics/paper_twas/inputs'  

## UKB
def sample_UKB(n_shuffles=100):
    uk_samples = 39565 
    full_idxs = np.arange(uk_samples).astype(int)

    samp_idxs = {} 
    disc_sizes = [772, 2000, 5000, 15000]
    for dsc in disc_sizes: 

        ds = uk_samples - dsc
        
        idxs = np.zeros((n_shuffles, ds), dtype=int)
        for i in range(n_shuffles): 
            idx = np.random.choice(full_idxs, ds, replace=False) 
            idxs[i] = np.sort(idx) 
        
        samp_idxs[dsc] = idxs
            
    opath = f'{main_path}_UKB/twas_repl_samples.hdf5' 
    with h5py.File(opath, 'w') as f: 
        for ds, idxs in samp_idxs.items(): 
            f[str(ds)] = idxs

## UKB: sample N = HCP/nonTwin cohort 
def sample_UKB_HCP_nonTwin(n_shuffles=100):
    
    ## load existing sampling 
    ipath = f'{main_path}_UKB/twas_repl_samples_HCPnonTwin.hdf5'
    with h5py.File(ipath, 'r') as f: 
        data = f['repl_idx'][()]

    ## sample
    uk_samples = 39565 
    full_idxs = np.arange(uk_samples).astype(int)

    n_nontwin = 657
    repl_idxs = np.zeros((n_shuffles, n_nontwin), dtype=int)
    for i in range(n_shuffles): 
        idx = np.random.choice(full_idxs, n_nontwin, replace=False) 
        repl_idxs[i] = np.sort(idx)

    ## concat old and new 
    repl_idxs = np.concatenate((repl_idxs, data), axis=0)
    opath = f'{main_path}_UKB/twas_repl_samples_HCPnonTwin.hdf5'
    with h5py.File(opath, 'w') as f: 
        f['repl_idx'] = repl_idxs

## UKB: sample N = HCP/allEuro cohort 
def sample_UKB_HCP_allEuro(n_shuffles=100):
    
    ## sample
    uk_samples = 39565 
    full_idxs = np.arange(uk_samples).astype(int)

    n_nontwin = 772
    repl_idxs = np.zeros((n_shuffles, n_nontwin), dtype=int)
    for i in range(n_shuffles): 
        idx = np.random.choice(full_idxs, n_nontwin, replace=False) 
        repl_idxs[i] = np.sort(idx)

    opath = f'{main_path}_UKB/twas_repl_samples_HCPallEuro.hdf5'
    with h5py.File(opath, 'w') as f: 
        f['repl_idx'] = repl_idxs
    
## HCP 
def sample_HCP(n_shuffles=100): 

    all_euro = 772
    non_twin = 657

    full_idxs = np.arange(all_euro).astype(int)
    samp_idxs = np.zeros((n_shuffles, non_twin), dtype=int)

    for i in range(n_shuffles): 
        idx = np.random.choice(full_idxs, non_twin, replace=False)
        samp_idxs[i] = np.sort(idx)

    opath = f'{main_path}_HCP/allEuro/twas_disc_samples.hdf5'
    with h5py.File(opath, 'w') as f: 
        f['772c'] = samp_idxs

#######################################################

#sample_UKB_HCP_nonTwin(400)
sample_UKB_HCP_allEuro(500)
    
