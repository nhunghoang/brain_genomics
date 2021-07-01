'''
Null model version of CCC Simulated Annealing (binary gene selection)  

This script (in isolation) runs SA for all 45 region pairs, 
for a given shuffle run. There are 100 shuffles, and each shuffle 
should be ran 100 times. 

This version does not constrain the number of genes to be 
selected. The initial weights are randomly decided. At every 
step, the weight of one random gene is considered.  

[06.29.21] CHANGE: Constrain the number of selected genes to be 32.

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
out_path = '/data/rubinov_lab/brain_genomics_project/platypus/results_ccc_null32'
#out_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/multi_gene_assoc/results_ccc_null'

## job array params 
shuffle_run = int(sys.argv[1])
SID = int(shuffle_run/100) ## shuffle version 
RID = shuffle_run%100 ## run of this shuffle 

## simulated annealing params 
TEMP = 0.01 
DECAY = 0.99 
DECAY_RATE = 1000 ## num of iterations before temp decay (-1 if it should be n_genes)
SAME_CCC = 1e5 ## converges after CCC has stayed the same for this number of steps
TOLERANCE = 1e-7 ## CCC is considered the same if difference is within this range
STEP_LIMIT = 1e9 ## converge or max out on this number of steps
#print('PARAMS: temp = {}, decay = {}'.format(TEMP, DECAY))

## get null shuffles
shuffle_file = '/data/rubinov_lab/brain_genomics_project/platypus/null_shuffles.hdf5'
#shuffle_file = '/data1/rubinov_lab/brain_genomics/analyses_HCP/multi_gene_assoc/null_shuffles.hdf5'
with h5py.File(shuffle_file, 'r') as f: 
    subj_order = np.array(f['subj_idx'])
    samp_order = np.array(f['samp_idx_{}'.format(SID)])

## gather data and rearrange based on null shuffle 
ccc_file = '/data/rubinov_lab/brain_genomics_project/platypus/ccc_constants.hdf5' 
#ccc_file = '/data1/rubinov_lab/brain_genomics/analyses_HCP/multi_gene_assoc/ccc_constants.hdf5' 
data = {c:{} for c in ['Gi', 'Gj', 'Fs', 'GiGj', 'Gi2', 'Gj2']}
with h5py.File(ccc_file, 'r') as f: 
    for key in f.keys(): 
        if key == 'IDs': continue 
        [reg_pair, c] = key.split('.') 
        if c == 'Fs': data[c][reg_pair] = np.array(f[key])[subj_order]
        else: data[c][reg_pair] = np.array(f[key])[samp_order] 

region_pairs = list(data['Fs'].keys())

## create directories on first iteration 
if (shuffle_run == 0) and (not os.path.exists(out_path)): 
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
    selected_idx = np.random.choice(n_genes, 32, replace=False)
    mask = np.ones(n_genes, dtype=bool)
    mask[selected_idx] = False
    unselected_idx = np.arange(n_genes)[mask]
    W = np.zeros(n_genes, dtype=int)
    W[selected_idx] = 1

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

        '''
        if step%(10*decay_rate) == 0:
            ccc_ = round(CCC,5)
            tmp_ = round(temp,5)
            wgt_ = np.sum(W)
            print('[STEP {:6d}] CCC: {:.5f} / TEMP: {:.5f} / NSEL: {:3d}'.format(step, ccc_, tmp_, wgt_))
        '''

        ## randomly select two weights to switch     
        sel2unsel = np.random.choice(selected_idx)
        unsel2sel = np.random.choice(unselected_idx)
        #w_idx = np.random.choice(n_genes) 

        ## delta calcs
        def delta_op(w_val, val, sel2unsel_idx, unsel2sel_idx): 
            return w_val - np.squeeze(val[:,sel2unsel_idx]) + np.squeeze(val[:,unsel2sel_idx])

        d_wGiGj = delta_op(wGiGj, GiGj, sel2unsel, unsel2sel)
        d_wGi2 = delta_op(wGi2, Gi2, sel2unsel, unsel2sel)
        d_wGj2 = delta_op(wGj2, Gj2, sel2unsel, unsel2sel)
        d_Gs = d_wGiGj / (np.sqrt(d_wGi2) * np.sqrt(d_wGj2))
        new_ccc, new_pval = pearsonr(Fs, d_Gs)

        ## update by chance or if good step
        old_ccc = CCC
        ccc_diff = CCC - new_ccc
        if (CCC < new_ccc) or (np.exp(-ccc_diff/temp) > np.random.rand()):
            W[sel2unsel] = 0 
            W[unsel2sel] = 1 
            tmp = sel2unsel 
            selected_idx[selected_idx==sel2unsel] = unsel2sel
            unselected_idx[unselected_idx==unsel2sel] = tmp 
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
            print('NONC - {}\n'.format(tpass))
            return CCC, pval, W, tpass, 'NONC'

    ## print details of final step before convergence 
    ccc_ = round(CCC,5)
    tmp_ = round(temp,5)
    wgt_ = np.sum(W)
    tpass = pretty_time(time()-tstart)
    print('[STEP {:6d}] CCC: {:.5f} / TEMP: {:.5f} / NSEL: {:3d}'.format(step, ccc_, tmp_, wgt_))
    print('CONV ({} / {} genes) - {}\n'.format(wgt_, n_genes, tpass))
    return CCC, pval, W, tpass, 'CONV - step {}'.format(step)

## run  
shuffle(region_pairs) 
for rpi, reg_pair in enumerate(region_pairs): 

    [reg1, reg2] = reg_pair.split('_')

    wgt = '{}/weights_{}_{}/{}_{}.hdf5'.format(out_path, reg1, reg2, SID, RID)
    log = '{}/logs_{}_{}/{}_{}.log'.format(out_path, reg1, reg2, SID, RID)
    if (os.path.exists(wgt)) and (os.path.exists(log)): 
        print('<{:>2d}> {} & {}'.format(rpi+1, reg1, reg2))
        continue 

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

    n_genes = GiGj.shape[1]
    info = '{:3d} genes: r = {:.5f} (p ≤ {:.3f})'.format(n_genes, nw_rho, nw_pval)
    print('<{:>2d}> {} & {} | {}'.format(rpi+1, reg1, reg2, info))

    ## simulated annealing
    w_rho, w_pval, W, tpass, conv_stat = simulated_annealing(reg_pair)
    
    ## save results  
    wgt = '{}/weights_{}_{}/{}_{}.hdf5'.format(out_path, reg1, reg2, SID, RID)
    log = '{}/logs_{}_{}/{}_{}.log'.format(out_path, reg1, reg2, SID, RID)

    with h5py.File(wgt, 'w') as f: 
        f['W'] = W

    with open(log, 'w') as f: 
        uw_line = '{:3d} genes: r = {:.10f} (p ≤ {:.3f})'.format(n_genes, nw_rho, nw_pval)
        wt_line = '{:3d} genes: r = {:.10f} (p ≤ {:.3f})'.format(np.sum(W), w_rho, w_pval)
        ct_line = '{} ({})\n'.format(conv_stat, tpass)  
        f.write('\n'.join([uw_line, wt_line, ct_line]))
