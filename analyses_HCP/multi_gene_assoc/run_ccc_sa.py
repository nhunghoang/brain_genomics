'''
CCC Simulated Annealing (binary gene selection)  

This version does not constrain the number of genes to be 
selected. The initial weights are randomly decided. At every 
step, the weight of one random gene is considered.  

- Nhung, updated June 2021 
'''

import numpy as np
from scipy.stats import pearsonr
import h5py
from time import time
import sys
from random import shuffle
import os

## output path 
out_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/multi_gene_assoc/results_ccc_obs'

## job array params 
FID = sys.argv[1]

## simulated annealing params 
TEMP = 0.01 
DECAY = 0.99 
DECAY_RATE = 1000 ## num of iterations before temp decay (-1 if it should be n_genes)
SAME_CCC = 1e5 ## converges after CCC has stayed the same for this number of steps
TOLERANCE = 1e-7 ## CCC is considered the same if difference is within this range
STEP_LIMIT = 1e9 ## converge or max out on this number of steps

print('PARAMS: temp = {}, decay = {}'.format(TEMP, DECAY))

## gather data
ccc_file = '/data1/rubinov_lab/brain_genomics/analyses_HCP/multi_gene_assoc/ccc_variables.hdf5' 
data = {c:{} for c in ['Gi', 'Gj', 'Fs', 'GiGj', 'Gi2', 'Gj2']}
with h5py.File(ccc_file, 'r') as f: 
    for key in f.keys(): 
        [reg_pair, c] = key.split('.') 
        data[c][reg_pair] = np.array(f[key]) 
region_pairs = list(data['Fs'].keys())

## create directories on first iteration 
if FID == '0': 
    os.mkdir(out_path)  
    for reg_pair in region_pairs: 
        [reg1, reg2] = reg_pair.split('_')
        os.mkdir('{}/weights_{}_{}'.format(out_path, reg1, reg2))
        os.mkdir('{}/logs_{}_{}'.format(out_path, reg1, reg2))

## pretty-print time 
def pretty_time(tpass):
    hr = int(tpass//3600)
    mn = int(tpass%3600//60)
    sc = int(tpass%3600%60)
    return '{:d} hr, {:d} mn, {:d} sc'.format(hr, mn, sc)

## SA function 
def simulated_annealing(reg_pair):
    
    ## constants 
    GiGj = data['GiGj'][reg_pair]
    Gi2 = data['Gi2'][reg_pair]
    Gj2 = data['Gj2'][reg_pair]
    Fs = data['Fs'][reg_pair] 

    ## decay rate 
    n_genes = GiGj.shape[1] 
    decay_rate = DECAY_RATE
    if DECAY_RATE == -1: decay_rate = n_genes 
    
    ## initialize weights 
    W = np.random.choice([0,1], n_genes)

    ## initial correlation 
    wGiGj = np.sum(GiGj*W, axis=1)
    wGi2 = np.sum(Gi2*W, axis=1)
    wGj2 = np.sum(Gj2*W, axis=1)
    Gs = wGiGj / (np.sqrt(wGi2) * np.sqrt(wGj2))
    CCC, pval = pearsonr(Fs,Gs)

    ## start
    tstart = time()
    step = 0
    same_CCC = 0
    temp = TEMP

    while (same_CCC < SAME_CCC): ## i.e. stop if CCC hasn't changed in k steps

        ## gradual temp decay
        if step%(decay_rate) == 0:
            temp *= DECAY

        if step%(10*decay_rate) == 0:
            ccc_ = round(CCC,5)
            tmp_ = round(temp,5)
            wgt_ = np.sum(W)
            print('[STEP {:6d}] CCC: {:.5f} / TEMP: {:.5f} / NSEL: {:3d}'.format(step, ccc_, tmp_, wgt_))

        ## randomly select one weight to flip    
        w_idx = np.random.choice(n_genes) 

        ## delta calcs
        def delta_op(w_val, val, idx): 
            if W[idx] == 0: return w_val + np.squeeze(val[:,idx])
            if W[idx] == 1: return w_val - np.squeeze(val[:,idx])

        d_wGiGj = delta_op(wGiGj, GiGj, w_idx)
        d_wGi2 = delta_op(wGi2, Gi2, w_idx)
        d_wGj2 = delta_op(wGj2, Gj2, w_idx)
        d_Gs = d_wGiGj / (np.sqrt(d_wGi2) * np.sqrt(d_wGj2))
        new_ccc, new_pval = pearsonr(Fs, d_Gs)

        ## update by chance or if good step
        old_ccc = CCC
        ccc_diff = CCC - new_ccc
        if (CCC < new_ccc) or (np.exp(-ccc_diff/temp) > np.random.rand()):
            W[w_idx] = 1 - W[w_idx] 
            CCC = new_ccc
            pval = new_pval
            wGiGj = d_wGiGj
            wGi2 = d_wGi2
            wGj2 = d_wGj2

        ## check if CCC has changed much
        if abs(old_ccc - CCC) < TOLERANCE: same_CCC += 1
        else: same_CCC = 0

        ## avoid running forever
        step += 1
        if step > STEP_LIMIT:
            ccc_ = round(CCC,5)
            tmp_ = round(temp,5)
            wgt_ = np.sum(W)
            tpass = pretty_time(time()-tstart)
            print('[STEP {:6d}] CCC: {:.5f} / TEMP: {:.5f} / NSEL: {:3d}'.format(step, ccc_, tmp_, wgt_))
            print('NONC - {}'.format(tpass))
            return CCC, pval, W, tpass, 'NONC'

    ccc_ = round(CCC,5)
    tmp_ = round(temp,5)
    wgt_ = np.sum(W)
    tpass = pretty_time(time()-tstart)
    print('[STEP {:6d}] CCC: {:.5f} / TEMP: {:.5f} / NSEL: {:3d}'.format(step, ccc_, tmp_, wgt_))
    print('CONV ({} / {} genes) - {}'.format(wgt_, n_genes, tpass))
    return CCC, pval, W, tpass, 'CONV - step {}'.format(step)

## run  
shuffle(region_pairs) 
for reg_pair in region_pairs: 

    [reg1, reg2] = reg_pair.split('_')
    n_genes = data['Gi'][reg_pair].shape[1]
    print('> {}, {} ({} genes)'.format(reg1, reg2, n_genes))

    ## original correlation, non-weighted
    GiGj = data['GiGj'][reg_pair]
    Gi2 = data['Gi2'][reg_pair]
    Gj2 = data['Gj2'][reg_pair]
    Fs = data['Fs'][reg_pair]

    nGiGj = np.sum(GiGj, axis=1)
    nGi2 = np.sum(Gi2, axis=1)
    nGj2 = np.sum(Gj2, axis=1)
    Gs = nGiGj / (np.sqrt(nGi2) * np.sqrt(nGj2))
    nw_rho, nw_pval = pearsonr(Fs,Gs)

    ## simulated annealing
    w_rho, w_pval, W, tpass, conv_stat = simulated_annealing(reg_pair)
    
    ## save results  
    wgt = '{}/weights_{}_{}/{}.hdf5'.format(out_path, reg1, reg2, FID)
    log = '{}/logs_{}_{}/{}.log'.format(out_path, reg1, reg2, FID)

    with h5py.File(wgt, 'w') as f: 
        f['W'] = W

    with open(log, 'w') as f: 
        uw_line = '{:3d} genes: r = {:.10f} (p ≤ {:.3f})'.format(n_genes, nw_rho, nw_pval)
        wt_line = '{:3d} genes: r = {:.10f} (p ≤ {:.3f})'.format(np.sum(W), w_rho, w_pval)
        ct_line = '{} ({})\n'.format(conv_stat, tpass)  
        f.write('\n'.join([uw_line, wt_line, ct_line]))
