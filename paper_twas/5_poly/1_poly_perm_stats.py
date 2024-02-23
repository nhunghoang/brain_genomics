'''
Permutation version of polygenic analysis. 

Per regional phenotype, a permutation is regression 
fit to a random set of genes (size equal to 
TWAS p < 0.001 genes). 

Run N permutations per regional phenotype. 

Save matrix of test r2s with shape (num_nulls, num_splits).

- Nhung, Oct 2023
'''

import numpy as np
import pandas as pd
import h5py
import sys

from scipy.stats import pearsonr

from statsmodels.regression.linear_model import OLS as OLS
from statsmodels.tools.tools import add_constant

from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score

from multiprocessing import Pool 

## params 
group = sys.argv[1] ## HCP, HCP/nonTwin, ...

num_split = None ## num of unique HCP fams
regs = ['dlpfc', 'anterior-cingulate', 'amygdala', 'hippocampus', \
        'caudate', 'putamen', 'nucleus-accumbens', 'cerebellar-hemisphere']
phens = ['vol_mean', 'alff_mean', 'reho_noGS_mean', 'connmean_noGS_mean']
reg_phens = [(r,p) for r in regs for p in phens]

## paths
main_path = '/data1/rubinov_lab/brain_genomics/paper_twas/'

subj_path = main_path + f'inputs_{group}/cohort.txt'
fams_path = main_path + 'inputs_HCP/demographics.csv'

twas_path = main_path + 'outputs_{}/twas_JTI/{}/{}.txt' ## group, phen, reg
grex_path = main_path + 'inputs_{}/grex_JTI/{}.hdf5' ## group, reg
phen_path = main_path + 'inputs_{}/phenotypes/{}_regr.csv' ## group, phen

## out paths
ntmp_path = main_path + f'outputs_{group}/polygenic_models/null_tmp' ## /phen_reg_itr.hdf5
null_path = main_path + f'outputs_{group}/polygenic_models/nulls.hdf5' ## keys: (reg, phen)

## load subject array
subj_arr = pd.read_csv(subj_path, sep=' ', header=None)[0].values

## load family IDs
df = pd.read_csv(fams_path, usecols=['subject', 'family_id'], index_col='subject')
df = df.reindex(subj_arr)
subj_fids = df['family_id'].values
uniq_fids = np.unique(subj_fids)

num_split = uniq_fids.size * 2

## split group into train and test sets
train_masks = {} ## k: itr, v: mask
for i in range(num_split):
    train_fids, _ = train_test_split(uniq_fids, \
                    test_size=0.2, \
                    random_state=i*150)
    train_masks[i] = np.isin(subj_fids, train_fids)

## load grex and phen data
grex_all = {} ## k: reg, v: (subj * gene)
phen_all = {} ## k: (reg, phen), v: (subj,)
gene_num = {} ## k: (reg, phen), v: num of sig genes 

for phen in phens: 
    pfile = phen_path.format(group, phen)
    pdf = pd.read_table(pfile, sep='\t')

    for reg in regs: 

        ## store regional phen values 
        phen_all[(reg,phen)] = pdf[reg].values

        ## store num of genes with TWAS p < 0.001
        tfile = twas_path.format(group, phen, reg)
        tdf = pd.read_table(tfile) 
        num = np.sum(tdf['pvalue'] < 0.001)

        gene_num[(reg,phen)] = num 
        print(f'{num} genes for {reg} {phen}')

        ## store regional grex
        if phen != phens[0]: continue 
        gfile = grex_path.format(group, reg)
        with h5py.File(gfile, 'r') as f:
            grex_mat = f['pred_expr'][()]
        grex_all[reg] = grex_mat.T

## func: fit regression model on training cohort
##       and evaluate on test cohort
def regression(itr, grex_mat, phen_arr):

    ## get train/test split 
    train_mask = train_masks[itr]
    testt_mask = ~train_mask

    ## get phen data
    ytrain = phen_arr[train_mask]
    ytestt = phen_arr[testt_mask]

    ## get grex data
    xtrain = grex_mat[train_mask]
    xtestt = grex_mat[testt_mask]

    ## for OLS: add intercept constant
    xtrain = add_constant(xtrain)
    xtestt = add_constant(xtestt)

    ## fit multi-variate regression to
    ## training data then get test data outputs
    model = OLS(ytrain, xtrain).fit()
    ptestt = model.predict(xtestt)

    ## fit regression to test data output
    ptestt = add_constant(ptestt)
    model2 = OLS(ytestt, ptestt).fit()
    r2 = model2.rsquared

    return r2

## func: permutation pool 
def pool_run(reg_phen_null_xtimes): 

    ## pool params 
    (reg, phen, null, xtimes) = reg_phen_null_xtimes 
    
    perm_rs = np.zeros((xtimes, num_split))
    
    for pidx in range(xtimes): 
        
        ## get a rand list of genes
        rng = np.random.RandomState()
        grex = grex_all[reg]
        idxs = np.arange(grex.shape[1])
        gidx = rng.choice(idxs, gene_num[(reg,phen)], replace=False)

        ## get grex and phen
        grex_mat = grex[:, gidx]
        phen_arr = phen_all[(reg,phen)]

        ## get pearsons based on splits
        perm_rs[pidx] = [regression(i, grex_mat, phen_arr) for i in range(num_split)]

    ## tmp save 
    tfile = f'{ntmp_path}/{phen}_{reg}_{null}.hdf5'
    with h5py.File(tfile, 'w') as f: 
        f['split_r2s'] = perm_rs 
    print(phen, reg, null, '\n')

## pool calls 
xtimes = 10
num_nulls = 100
params = [(r,p,i,xtimes) for r in regs for p in phens for i in range(num_nulls)]

pool = Pool(processes=250)
pool.map(pool_run, params)

## concat results  
stats = {} ## k: (regXphen), v: (nulls, splits)
for (reg, phen) in reg_phens:

    rstats = np.zeros((num_nulls * xtimes, num_split))
    for n in range(num_nulls): 
        rfile = f'{ntmp_path}/{phen}_{reg}_{n}.hdf5'
        with h5py.File(rfile, 'r') as f: 
            s = n*xtimes; e = s+xtimes 
            rstats[s:e] = f['split_r2s'][()]

    stats[(reg,phen)] = rstats

with h5py.File(null_path, 'w') as f:
    for (reg, phen), rstats in stats.items():
        f[f'{reg}X{phen}'] = rstats

