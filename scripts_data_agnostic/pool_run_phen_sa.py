'''
Regional Phenotype Simulated Annealing (binary gene selection) 

This version searches for genes using the training cohort, 
and then applies the solution to the test cohort. The number 
of genes is unconstrained. 

This script uses Pool to run in parallel. 

- Nhung, updated Oct 2021

Notes: 
pass in phen and reg 
and fid index 

'''

import numpy as np
from scipy.stats import pearsonr, spearmanr
import h5py
from time import time
import sys
from random import shuffle 
import os
from multiprocessing import Pool

## job array params
ATL = 'hoacer_sn_hth' ## atlas used in calculating phenotype 
SEL = -1 ## number of genes to select (-1 means unrestrained) 
#PHN = sys.argv[1] ## phenotype 
#REG = sys.argv[2] ## region 
#IDX = int(sys.argv[3]) ## index wrt to FID order 

## input paths 
#phen_file = '/data1/rubinov_lab/brain_genomics/analyses_HCP/DATA_OUTPUT/phen_regress/{}.hdf5'.format(PHN)
expr_dir = '/data1/rubinov_lab/brain_genomics/analyses_HCP/DATA_OUTPUT/expr_regress' 

## simulated annealing params
TEMP = 0.01
DECAY = 0.99
DECAY_RATE = 1000 ## num of iterations before temp decay (-1 if it should be n_genes)
SAME_RHO = 1e5 ## converges after rho has stayed the same for this number of steps
TOLERANCE = 1e-7 ## rho is considered the same if difference is within this range
STEP_LIMIT = 1e9 ## converge or max out on this number of steps

## function: pretty-print time
def pretty_time(tpass):
    hr = int(tpass//3600)
    mn = int(tpass%3600//60)
    sc = int(tpass%3600%60)
    return '{:d} hr, {:d} mn, {:d} sc'.format(hr, mn, sc)

########################################################################################

## function: simulated annealing 
def simulated_annealing(phen_array, gene_array, expr_matrx): 
    
    ## data - NO TRAIN/TEST SPLIT 
    n_genes = gene_array.shape[0]
    
    ## initialize weights 
    W = np.random.choice([0,1], n_genes) 

    ## initial correlation 
    expr_array = np.matmul(expr_matrx, W) ## (N subjects * 1)  
    rho, pval = pearsonr(expr_array, phen_array)

    ## decay rate 
    decay_rate = DECAY_RATE
    if DECAY_RATE == -1: decay_rate = n_genes

    ## start 
    tstart = time()
    step = 0 
    same_rho = 0
    temp = TEMP

    while (same_rho < SAME_RHO): ## i.e. stop if rho hasn't changed in k steps 
        
        ## gradual temp decay 
        if step%(decay_rate) == 0:
            temp *= DECAY 

        ## print search status 
        #if step%(10*decay_rate) == 0: 
        #    rho_ = round(rho,5)
        #    tmp_ = round(temp,5)
        #    wgt_ = np.sum(W)
        #    print('[STEP {:6d}] RHO: {:.5f} / TEMP: {:.5f} / NSEL: {:3d}'.format(step, rho_, tmp_, wgt_))

        ## randomly select one weight to flip 
        w_idx = np.random.choice(n_genes) 

        ## delta calcs 
        if W[w_idx] == 0: d_expr_array = expr_array + np.squeeze(expr_matrx[:,w_idx])
        if W[w_idx] == 1: d_expr_array = expr_array - np.squeeze(expr_matrx[:,w_idx])

        new_rho, new_pval = pearsonr(d_expr_array, phen_array)

        ## update by chance or if good step  
        old_rho = rho
        rho_diff = rho - new_rho
        if (rho < new_rho) or (np.exp(-rho_diff/temp) > np.random.rand()):
            W[w_idx] = 1 - W[w_idx] 
            rho = new_rho
            pval = new_pval
            expr_array = d_expr_array
        
        ## check if rho has changed much
        if abs(old_rho - rho) < TOLERANCE: same_rho += 1 
        else: same_rho = 0 

        ## avoid running forever
        step += 1 
        if step > STEP_LIMIT:
            rho_ = round(rho,5)
            tmp_ = round(temp,5)
            wgt_ = np.sum(W)
            tpass = pretty_time(time()-tstart)
            #print('[STEP {:6d}] RHO: {:.5f} / TEMP: {:.5f} / NSEL: {:3d}'.format(step, rho_, tmp_, wgt_))
            #print('NONC - {}\n'.format(tpass))
            return rho, pval, W, tpass, 'NONC'

    ## print details of final step before convergence
    #rho_ = round(rho,5)
    #tmp_ = round(temp,5)
    #wgt_ = np.sum(W)
    tpass = pretty_time(time()-tstart)
    #print('[STEP {:6d}] RHO: {:.5f} / TEMP: {:.5f} / NSEL: {:3d}'.format(step, rho_, tmp_, wgt_))
    #print('CONV ({} / {} genes) - {}'.format(wgt_, n_genes, tpass))

    return rho, pval, W, tpass, 'CONV - step {}'.format(step)

########################################################################################

## ONE RUN: lofo on one reg-phen  
def simann(data_dict):

    REG = data_dict['REG']
    PHN = data_dict['PHN']

    ## parse single-assoc for (p <= 0.05) gene indices 
    pval_file = '/data1/rubinov_lab/brain_genomics/analyses_HCP/assoc/pvals_{}/{}.txt'.format(PHN, REG) 
    pval_data = np.loadtxt(pval_file, delimiter='\t', skiprows=1, usecols=[2])
    pval_mask = np.zeros_like(pval_data, dtype=bool) 
    pval_mask[pval_data <= 0.05] = True 

    ## input data: reg phenotype array  
    phen_file = '/data1/rubinov_lab/brain_genomics/analyses_HCP/DATA_OUTPUT/phen_regress/{}.hdf5'.format(PHN)
    with h5py.File(phen_file, 'r') as f:
        phen_array = np.array(f[REG]) ## (subjects,) 

    ## input data: expression (only for significant genes)  
    with h5py.File('{}/{}.hdf5'.format(expr_dir, REG), 'r') as f: 
        expr_matrx = np.array(f['pred_expr'])[pval_mask].T ## (subjects * genes) 

    ## genes to search over 
    gene_list = np.loadtxt(pval_file, delimiter='\t', skiprows=1, usecols=[0], dtype=str)
    genes = gene_list[pval_mask] ## (genes,)
    n_genes = genes.size 

    print('{:>20s} | {:>21}: {} genes'.format(PHN, REG, expr_matrx.shape[1]))

    ## non-weighted correlations  
    expr_array = np.sum(expr_matrx, axis=1) ## (N subjects,) 
    rho, pval = pearsonr(expr_array, phen_array)

    ## simulated annealing  
    w_rho, w_pval, W, tpass, conv_stat = simulated_annealing(phen_array, genes, expr_matrx)

    ## save all results
    status_line = '{} ({}) [{}/{} genes selected]'.format(conv_stat, tpass, np.sum(W).astype(int), np.size(W))
    all_line = '{:.10f} (p ≤ {:.3f}) | {:.10f} (p ≤ {:.3f})'.format(rho, pval, w_rho, w_pval)
    placeholder = 'N/A (p ≤ N/A) | N/A (p ≤ N/A)'
    weights = '\t'.join(W.astype(str)) + '\n' 

    out_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/assoc/simann'
    log = '{}_{}/{}.log'.format(out_path, PHN, REG)
    with open(log, 'w') as f:
        f.write('\n'.join([status_line, all_line, placeholder, weights]))

    ## print results 
    n_selc_weights = np.sum(W).astype(int) 
    all_line1 = '{:.3f} (p ≤ {:.3f}) | {:.3f} (p ≤ {:.3f})'.format(rho, pval, w_rho, w_pval)
    print(status_line)
    print(all_line1) 

########################################################################################

regs = os.listdir('/data1/rubinov_lab/brain_genomics/analyses_HCP/lofo_assoc/train_pvals_alff')
phens = ['alff', 'regional_homogeneity', 'gm_volume']
reg_phen = [{'REG':reg, 'PHN':phen} for reg in regs for phen in phens] 
pool = Pool(processes=10)
pool.map(simann, reg_phen)  

