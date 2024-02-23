'''
Polygenic analysis (observed stats). 

Per regional phenotype, fit a regression 
to the genes with TWAS p < 0.001 using 
an 80/20 train/test split on the family units. 

The number of train/test splits equals twice 
the number of HCP families.  

Store the following stats per regional phenotype: 
- matrix with dimensions (test subjs * [ytrue, ypred]) 
  for the median model (median based on test r2) 
- array of test r2s with dimensions (num_split) 
- array of test pvals with dimensions (num_split) 

- Nhung, Oct 2023
'''

import numpy as np
import pandas as pd
import h5py
import sys 

from statsmodels.regression.linear_model import OLS as OLS
from statsmodels.tools.tools import add_constant
from statsmodels.stats.multitest import multipletests as sm

from sklearn.model_selection import train_test_split

group = sys.argv[1] ## HCP, HCP/nonTwin, ...

num_split = None ## num of unique HCP fams
regions = ['dlpfc', 'anterior-cingulate', 'amygdala', 'hippocampus', \
           'caudate', 'putamen', 'nucleus-accumbens', 'cerebellar-hemisphere']
phens = ['vol_mean', 'alff_mean', 'reho_noGS_mean', 'connmean_noGS_mean']
reg_phens = [(r,p) for r in regions for p in phens]

## paths
main_path = '/data1/rubinov_lab/brain_genomics/paper_twas/'

subj_path = main_path + f'inputs_{group}/cohort.txt'
fams_path = main_path + 'inputs_HCP/demographics.csv'

twas_path = main_path + 'outputs_{}/twas_JTI/{}/{}.txt' ## group, phen, reg
grex_path = main_path + 'inputs_{}/grex_JTI/{}.hdf5' ## group, reg
phen_path = main_path + 'inputs_{}/phenotypes/{}_regr.csv' ## group, phen

## out paths 
afile = main_path + f'outputs_{group}/polygenic_models/median_ytrue_ypred.hdf5'
bfile = main_path + f'outputs_{group}/polygenic_models/split_r2s.hdf5'
cfile = main_path + f'outputs_{group}/polygenic_models/split_pvs.hdf5'

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
grex_all = {} ## k: (reg, phen), v: (subj * gene)
phen_all = {} ## k: (reg, phen), v: (subj,)

for (reg, phen) in reg_phens: 

    ## get significant TWAS genes 
    tfile = twas_path.format(group, phen, reg)
    df = pd.read_table(tfile, sep='\t', usecols=['gene', 'pvalue'])

    keep = (df['pvalue'] < 0.001)
    df_keep = df.loc[keep]
    twas_genes = df_keep['gene'].values

    print(keep.sum(), reg, phen)

    ## get grex matrices and slice 
    gpath = grex_path.format(group, reg)
    with h5py.File(gpath, 'r') as f:
        reg_gene = f['genes'][()].astype(str)
        grex_mat = f['pred_expr'][()]

    gset, idx1, idx2 = np.intersect1d(reg_gene, twas_genes, return_indices=True)
    grex_all[(reg, phen)] = grex_mat[idx1].T

    ## get phenotype arrays 
    ppath = phen_path.format(group, phen)
    df = pd.read_table(ppath, sep='\t', usecols=[reg])
    phen_all[(reg,phen)] = df[reg].values

## func: fit regression model on training cohort
##       and evaluate on test cohort
def regression(itr, reg, phen):

    train_mask = train_masks[itr]
    testt_mask = ~train_mask

    ## get phen data
    ytrain = phen_all[(reg,phen)][train_mask]
    ytestt = phen_all[(reg,phen)][testt_mask]

    ## get grex data
    xtrain = grex_all[(reg,phen)][train_mask]
    xtestt = grex_all[(reg,phen)][testt_mask]

    ## for OLS: add intercept constant
    xtrain = add_constant(xtrain)
    xtestt = add_constant(xtestt)

    ## fit multi-variate regression to 
    ## training data then get test data outputs
    model = OLS(ytrain, xtrain).fit()
    ptestt = model.predict(xtestt) 

    ## evaluate test data output  
    ptestt = add_constant(ptestt)
    model2 = OLS(ytestt, ptestt).fit()
    r2 = model2.rsquared 
    pv = model2.pvalues[1] 

    results = {'r2': r2, 'pv': pv, \
               'ytest': ytestt, \
               'ypred': ptestt[:,1]}

    return results 

## main loop 
panel_a = {} ## k: (reg, phen), v: test subjs * [ytrue, ypred]
panel_b = {} ## k: (reg, phen), v: split r2s 
panel_c = {} ## k: (reg, phen), v: split pvs (will plot median)

for (reg, phen) in reg_phens: 

    split_r2s = np.zeros(num_split)
    split_pvs = np.zeros(num_split)
    split_ytests = {} ## k: split num, v: [ytests]
    split_ypreds = {} ## k: split num, v: [ypreds]

    ## regression on splits 
    results = [regression(i, reg, phen) for i in range(num_split)]
    for i in range(num_split): 
        split_r2s[i] = results[i]['r2']
        split_pvs[i] = results[i]['pv']
        split_ytests[i] = results[i]['ytest']
        split_ypreds[i] = results[i]['ypred']
    
    ## get median model 
    midx = np.argsort(split_r2s)[321]
    med_r2 = split_r2s[midx]
    med_pv = split_pvs[midx]
    med_ytest = split_ytests[midx][:,None]
    med_ypred = split_ypreds[midx][:,None]

    print('med r2: {:.3f}, med pv: {} ({} {})'.format(med_r2, med_pv, reg, phen))

    ## store data 
    panel_a[(reg,phen)] = np.concatenate([med_ytest, med_ypred], axis=1)
    panel_b[(reg,phen)] = split_r2s 
    panel_c[(reg,phen)] = split_pvs 

a = h5py.File(afile, 'w')
b = h5py.File(bfile, 'w')
c = h5py.File(cfile, 'w')

for (reg, phen) in reg_phens: 
    key = '{}X{}'.format(reg,phen)
    a[key] = panel_a[(reg,phen)]
    b[key] = panel_b[(reg,phen)]
    c[key] = panel_c[(reg,phen)]

_ = [f.close() for f in [a,b,c]]

