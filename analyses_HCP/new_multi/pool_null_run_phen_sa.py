'''
(Null) Regional Phenotype Simulated Annealing (binary gene selection) 

This version searches for genes using the training cohort,
and then applies the solution to the train cohort. The number
of genes is unconstrained.

This script uses Pool to run in parallel.

- Nhung, updated Oct 2021
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
#ATL = sys.argv[1] ## atlas used in calculating phenotype 
#SEL = int(sys.argv[2]) ## number of genes to select
#PHN = sys.argv[3] ## phenotype 
ATL = 'hoacer_sn_hth' 
SEL = -1
PHN = sys.argv[1]

perm_file = '/data1/rubinov_lab/brain_genomics/analyses_HCP/null_permutations.hdf5' 
#perm_file = '/data/rubinov_lab/brain_genomics_project/platypus2/null_permutations.hdf5' 

#shuffle_run = int(sys.argv[4]) 
#SID = int(shuffle_run/10) ## shuffle version 
#RID = shuffle_run%10 ## run of this shuffle 

## output path 
#out_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/multi_gene_assoc/obsv_{}'.format(PHN)
#out_path = '/data/rubinov_lab/brain_genomics_project/platypus2/null_{}'.format(PHN)
out_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/new_multi/null_results_{}'.format(PHN)
if not os.path.exists(out_path): os.mkdir(out_path) 

## simulated annealing params
TEMP = 0.01
DECAY = 0.99
DECAY_RATE = 1000 ## num of iterations before temp decay (-1 if it should be n_genes)
SAME_RHO = 1e5 ## converges after rho has stayed the same for this number of steps
TOLERANCE = 1e-7 ## rho is considered the same if difference is within this range
STEP_LIMIT = 1e9 ## converge or max out on this number of steps

## input data: phenotypes 
phen_file = '/data1/rubinov_lab/brain_genomics/data_HCP/{}/phenotypes/{}.hdf5'.format(ATL, PHN)
#phen_file = '/data/rubinov_lab/brain_genomics_project/platypus2/{}/phenotypes/{}.hdf5'.format(ATL, PHN)
with h5py.File(phen_file, 'r') as f:
    phens = {reg: np.array(f[reg]) for reg in f.keys()} ## k: reg, v: (subjects,)

regions = np.array(list(phens.keys()))

## input data: expression 
#expr_dir = '/data1/rubinov_lab/brain_genomics/data_HCP/expression/filtered_quality_r0.3_p0.01'
#samp_file = '/data1/rubinov_lab/brain_genomics/data_HCP/predixcan_v8/sample_ids.txt'
#expr_dir = '/data/rubinov_lab/brain_genomics_project/platypus2/filtered_quality_r0.3_p0.01'
#samp_file = '/data/rubinov_lab/brain_genomics_project/platypus2/sample_ids.txt'
expr_dir = '/data1/rubinov_lab/brain_genomics/data_HCP/expression/filtered_quality_r0.3_p0.01/cohort890_single-assoc'

#genes = {} ## k: reg, v: (genes,)
#exprs = {} ## k: reg, v: (subjects * genes)
#for expr_file in os.listdir(expr_dir):
#    if expr_file[-5:] != '.hdf5': continue
#    with h5py.File('{}/{}'.format(expr_dir, expr_file), 'r') as f:
#        reg = expr_file.split('.')[0]
#        genes[reg] = np.array([g.decode("utf-8") for g in np.array(f['genes'])])
#        exprs[reg] = np.array(f['pred_expr']).T[samp_order] 

genes = {} ## k: reg, v: (genes,)
exprs = {} ## k: reg, v: (subjects * genes)
for reg in regions:
    with h5py.File('{}/{}_{}.hdf5'.format(expr_dir, PHN, reg), 'r') as f:
        genes[reg] = np.array([g.decode('utf-8') for g in np.array(f['genes'])])
        exprs[reg] = np.array(f['expr']).T

## input data: order  
#order_file = '/data1/rubinov_lab/brain_genomics/analyses_HCP/subj_samp_assoc_order.hdf5'
#order_file = '/data/rubinov_lab/brain_genomics_project/platypus2/subj_samp_assoc_order.hdf5'
#with h5py.File(order_file, 'r') as f:
#    samp_order = np.array(f['sample_idx_1142'])

## function: pretty-print time
def pretty_time(tpass):
    hr = int(tpass//3600)
    mn = int(tpass%3600//60)
    sc = int(tpass%3600%60)
    return '{:d} hr, {:d} mn, {:d} sc'.format(hr, mn, sc)

########################################################################################

## function: simulated annealing 
def simulated_annealing(reg, train_subj_idx, train_samp_idx): 
    
    ## data - TRAIN GROUP ONLY 
    phen_array = phens[reg][train_subj_idx] ## (N subjects * 1)
    gene_array = genes[reg] ## (M genes * 1) 
    expr_matrx = exprs[reg][train_samp_idx] ## (N subjects * M genes) 
    n_genes = gene_array.shape[0]
    
    ## initialize weights 
    W = np.random.choice([0,1], n_genes)
    #selected_idx = np.random.choice(n_genes, SEL, replace=False)
    #mask = np.ones(n_genes, dtype=bool)
    #mask[selected_idx] = False
    #unselected_idx = np.arange(n_genes)[mask]

    #W = np.zeros(n_genes, dtype=int) ## (M genes * 1)
    #W[selected_idx] = 1

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

        '''
        if step%(10*decay_rate) == 0: 
            rho_ = round(rho,5)
            tmp_ = round(temp,5)
            wgt_ = np.sum(W)
            print('[STEP {:6d}] RHO: {:.5f} / TEMP: {:.5f} / NSEL: {:3d}'.format(step, rho_, tmp_, wgt_))
        '''

        '''
        ## randomly select two weights to switch 
        sel2unsel = np.random.choice(selected_idx)
        unsel2sel = np.random.choice(unselected_idx)

        ## delta calcs 
        d_expr_array = expr_array \
            - np.squeeze(expr_matrx[:,sel2unsel]) \
            + np.squeeze(expr_matrx[:,unsel2sel])

        new_rho, new_pval = pearsonr(d_expr_array, phen_array)
        '''

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
            #W[sel2unsel] = 0
            #W[unsel2sel] = 1 
            #tmp = sel2unsel 
            #selected_idx[selected_idx==sel2unsel] = unsel2sel
            #unselected_idx[unselected_idx==unsel2sel] = tmp
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
    rho_ = round(rho,5)
    tmp_ = round(temp,5)
    wgt_ = np.sum(W)
    tpass = pretty_time(time()-tstart)
    #print('[STEP {:6d}] RHO: {:.5f} / TEMP: {:.5f} / NSEL: {:3d}'.format(step, rho_, tmp_, wgt_))
    #print('CONV ({} / {} genes) - {}'.format(wgt_, n_genes, tpass))

    return rho, pval, W, tpass, 'CONV - step {}'.format(step)

########################################################################################

## one shuffle run: one phenotype, all regions, 
def gene_search(shuffle_run):
    for rdx, reg in enumerate(regions):

        SID = int(shuffle_run/10) ## shuffle version 
        RID = shuffle_run%10 ## run of this shuffle 
    
        ## permutation indices
        with h5py.File(perm_file, 'r') as f: 
            train_subj_idx = np.array(f['subj_train_idx'])
            train_samp_idx = np.array(f['samp_train_idx_{}'.format(SID)])
            test_subj_idx = np.array(f['subj_test_idx'])
            test_samp_idx = np.array(f['samp_test_idx_{}'.format(SID)])

        ## train/test split 
        n_genes = genes[reg].shape[0]
        expr_matrx = exprs[reg] ## (N subjects * M genes) 
        phen_array = phens[reg] ## (N subjects,)
    
        train_expr = expr_matrx[train_samp_idx] 
        train_phen = phen_array[train_subj_idx] 
        
        test_expr = expr_matrx[test_samp_idx]
        test_phen = phen_array[test_subj_idx]
    
        #print('<{:>2d}> {} | {}'.format(rdx+1, reg, expr_matrx.shape[1]))
    
        ## non-weighted correlations  
        train_array = np.sum(train_expr, axis=1) ## (N subjects,)
        ntrain_rho, ntrain_pval = pearsonr(train_array, train_phen)
    
        test_array = np.sum(test_expr, axis=1) ## (N subjects,)
        ntest_rho, ntest_pval = pearsonr(test_array, test_phen)
    
        ## simulated annealing - train results 
        wtrain_rho, wtrain_pval, W, tpass, conv_stat = simulated_annealing(reg, train_subj_idx, train_samp_idx)
    
        ## test results 
        wtest_array = np.matmul(test_expr, W) ## (N subjects,)
        wtest_rho, wtest_pval = pearsonr(wtest_array, test_phen)
    
        ## save all results
        status_line = '{} ({}) [{}/{} genes selected]'.format(conv_stat, tpass, np.sum(W).astype(int), np.size(W))
        train_line = '{:.10f} (p ≤ {:.3f}) | {:.10f} (p ≤ {:.3f})'.format(ntrain_rho, ntrain_pval, wtrain_rho, wtrain_pval)
        test_line = '{:.10f} (p ≤ {:.3f}) | {:.10f} (p ≤ {:.3f})'.format(ntest_rho, ntest_pval, wtest_rho, wtest_pval)
        weights = '\t'.join(W.astype(str)) + '\n' 
    
        log_dir = '{}/{}'.format(out_path, reg)
        if not os.path.exists(log_dir): os.mkdir(log_dir) 
    
        log = '{}/unconstrained_{}_{}.log'.format(log_dir, SID, RID)
        with open(log, 'w') as f:
            f.write('\n'.join([status_line, train_line, test_line, weights]))
    
        ## print results 
        train_line = '{:.3f} (p ≤ {:.3f}) | {:.3f} (p ≤ {:.3f})'.format(ntrain_rho, ntrain_pval, wtrain_rho, wtrain_pval)
        test_line = '{:.3f} (p ≤ {:.3f}) | {:.3f} (p ≤ {:.3f})'.format(ntest_rho, ntest_pval, wtest_rho, wtest_pval)
        #print(status_line)
        #print(' train: {}'.format(train_line))
        #print(' test: {}'.format(test_line))

########################################################################################

shuffle_runs = [f for f in range(1000)]
pool = Pool(processes=17)
pool.map(gene_search, shuffle_runs) 
