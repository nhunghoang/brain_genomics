'''
Regional Phenotype Simulated Annealing (binary gene selection) 

This version constrains the number of selected genes to be 32. 
The initial weights are randomly decided. At every step, the 
proposed solution considers swapping one selected gene for one 
unselected gene. 

- Nhung Hoang, updated July 2021
'''

import numpy as np
from scipy.stats import pearsonr
import h5py
from time import time
import sys
from random import shuffle 
import os

## job array params
PHN = sys.argv[1] ## phenotype 
FID = sys.argv[2] ## independent run 

## output path 
out_path = '/data/rubinov_lab/brain_genomics_project/platypus/results_{}_obs32'.format(PHN)
#out_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/multi_gene_assoc/results_{}_obs32'.format(PHN)

## simulated annealing params
TEMP = 0.01
DECAY = 0.99
DECAY_RATE = 1000 ## num of iterations before temp decay (-1 if it should be n_genes)
SAME_RHO = 1e5 ## converges after RHO has stayed the same for this number of steps
TOLERANCE = 1e-7 ## RHO is considered the same if difference is within this range
STEP_LIMIT = 1e9 ## converge or max out on this number of steps
#print('PARAMS: temp = {}, decay = {}'.format(TEMP, DECAY))

## gather data 
order_file = '/data1/rubinov_lab/brain_genomics/analyses_HCP/subj_samp_assoc_order.hdf5' 
with h5py.File(order_file, 'r') as f: 
    subj_idx = np.array(f['subject_idx'])
    samp_idx = np.array(f['sample_idx'])

phen_file = '/data1/rubinov_lab/brain_genomics/data_HCP/phenotypes/{}.hdf5'.format(PHN)
with h5py.File(phen_file, 'r') as f:
    phens = {reg: np.array(f[reg])[subj_idx] for reg in f.keys()}

expr_dir = '/data1/rubinov_lab/brain_genomics/data_HCP/expression/filtered_quality_r0.3_p0.01' 
genes = {}; exprs = {}
for expr_file in os.listdir(expr_dir):
    if expr_file[-5:] != '.hdf5': continue
    with h5py.File('{}/{}'.format(expr_dir, expr_file), 'r') as f:
        reg = expr_file.split('.')[0]
        genes[reg] = np.array([g.decode("utf-8") for g in np.array(f['genes'])])
        exprs[reg] = np.array(f['pred_expr'])[:,samp_idx].T ## (samps * genes)

regions = np.array(list(phens.keys()))[:1]

## create result directories on first iteration 
if FID == '0':
    try: 
        os.mkdir(out_path)
        for reg in regions: 
            os.mkdir('{}/weights_{}'.format(out_path, reg))
            os.mkdir('{}/logs_{}'.format(out_path, reg))
    except: 
        print("Result directories have already been created.")

## function: pretty-print time
def pretty_time(tpass):
    hr = int(tpass//3600)
    mn = int(tpass%3600//60)
    sc = int(tpass%3600%60)
    return '{:d} hr, {:d} mn, {:d} sc'.format(hr, mn, sc)

## function: simulated annealing 
def simulated_annealing(reg): 
    
    ## data 
    phen_array = phens[reg] ## (N subjects * 1)
    gene_array = genes[reg] ## (N subjects * M genes) 
    expr_matrx = exprs[reg] ## (N subjects * M genes) 
    n_genes = gene_array.shape[0]
    
    ## initialize weights (constrain number of genes to 32)
    selected_idx = np.random.choice(n_genes, 32, replace=False)
    mask = np.ones(n_genes, dtype=bool)
    mask[selected_idx] = False
    unselected_idx = np.arange(n_genes)[mask]

    W = np.zeros(n_genes, dtype=int) ## (M genes * 1)
    W[selected_idx] = 1

    ## initial correlation 
    expr_array = np.matmul(expr_matrx, W) ## (N subjects * 1)  
    RHO, pval = pearsonr(expr_array, phen_array)

    ## decay rate 
    decay_rate = DECAY_RATE
    if DECAY_RATE == -1: decay_rate = n_genes

    ## start 
    tstart = time()
    step = 0 
    same_RHO = 0
    temp = TEMP 

    while (same_RHO < SAME_RHO): ## i.e. stop if RHO hasn't changed in k steps 
        
        ## gradual temp decay 
        if step%(decay_rate) == 0:
            temp *= DECAY 

        if step%(10*decay_rate) == 0: 
            rho_ = round(RHO,5)
            tmp_ = round(temp,5)
            wgt_ = np.sum(W)
            print('[STEP {:6d}] RHO: {:.5f} / TEMP: {:.5f} / NSEL: {:3d}'.format(step, rho_, tmp_, wgt_))

        ## randomly select two weights to switch 
        sel2unsel = np.random.choice(selected_idx)
        unsel2sel = np.random.choice(unselected_idx)

        ## delta calcs 
        d_expr_array = expr_array \
            - np.squeeze(expr_matrx[:,sel2unsel]) \
            + np.squeeze(expr_matrx[:,unsel2sel])

        new_rho, new_pval = pearsonr(d_expr_array, phen_array) 

        ## update by chance or if good step  
        old_rho = RHO
        rho_diff = RHO - new_rho
        if (RHO < new_rho) or (np.exp(-rho_diff/temp) > np.random.rand()):
            W[sel2unsel] = 0
            W[unsel2sel] = 1 
            tmp = sel2unsel 
            selected_idx[selected_idx==sel2unsel] = unsel2sel
            unselected_idx[unselected_idx==unsel2sel] = tmp
            RHO = new_rho
            pval = new_pval
            expr_array = d_expr_array
        
        ## check if RHO has changed much
        if abs(old_rho - RHO) < TOLERANCE: same_RHO += 1 
        else: same_RHO = 0 

        ## avoid running forever
        step += 1 
        if step > STEP_LIMIT:
            rho_ = round(RHO,5)
            tmp_ = round(temp,5)
            wgt_ = np.sum(W)
            tpass = pretty_time(time()-tstart)
            print('[STEP {:6d}] RHO: {:.5f} / TEMP: {:.5f} / NSEL: {:3d}'.format(step, rho_, tmp_, wgt_))
            print('NONC - {}\n'.format(tpass))
            return RHO, pval, W, tpass, 'NONC'

    ## print details of final step before convergence
    rho_ = round(RHO,5)
    tmp_ = round(temp,5)
    wgt_ = np.sum(W)
    tpass = pretty_time(time()-tstart)
    print('[STEP {:6d}] RHO: {:.5f} / TEMP: {:.5f} / NSEL: {:3d}'.format(step, rho_, tmp_, wgt_))
    print('CONV ({} / {} genes) - {}\n'.format(wgt_, n_genes, tpass))
    return RHO, pval, W, tpass, 'CONV - step {}'.format(step)

## run 
shuffle(regions)
for rdx, reg in enumerate(regions): 

    ## original correlation, non-weighted 
    n_genes = genes[reg].shape[0]  
    expr_array = np.sum(exprs[reg], axis=1)
    phen_array = phens[reg] 
    nw_rho, nw_pval = pearsonr(expr_array, phen_array) 

    info = '{:3d} genes: r = {:.5f} (p ≤ {:.3f})'.format(n_genes, nw_rho, nw_pval)
    print('<{:>2d}> {} | {}'.format(rdx+1, reg, info))

    ## simulated annealing 
    w_rho, w_pval, W, tpass, conv_stat = simulated_annealing(reg) 
    
    ## save results
    wgt = '{}/weights_{}/{}.hdf5'.format(out_path, reg, FID)
    log = '{}/logs_{}/{}.log'.format(out_path, reg, FID)

    with h5py.File(wgt, 'w') as f:
        f['W'] = W

    with open(log, 'w') as f:
        uw_line = '{:3d} genes: r = {:.10f} (p ≤ {:.3f})'.format(n_genes, nw_rho, nw_pval)
        wt_line = '{:3d} genes: r = {:.10f} (p ≤ {:.3f})'.format(np.sum(W), w_rho, w_pval)
        ct_line = '{} ({})\n'.format(conv_stat, tpass)
        f.write('\n'.join([uw_line, wt_line, ct_line]))
