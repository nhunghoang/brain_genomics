'''
For each regional phenotype, grab the genes 
at single-gene significance and fit a linear 
regression model. Split the cohort into 
train/test sets for model evaluation. 

The test set evaluation is highly dependent on 
which samples are in the test set (because the 
number of samples overall is so low). To counter 
this, run 10k iterations of the regression model, 
based on a different random split each time. This
is not to be confused with permutation testing -- 
we are only looking at the observed associations   

- Nhung, updated Feb 2023
'''

import h5py
import numpy as np 
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error
from scipy.stats import pearsonr
import matplotlib.pyplot as plt

## regions and phenotypes
phens = ['gm_volume', 'myelination', 'alff', 'reho_noGS', 'connmean_noGS']
regns = ['frontal-pole', 'anterior-cingulate', 'caudate', 'putamen', \
         'nucleus-accumbens', 'hippocampus', 'amygdala', 'hypothalamus', \
         'substantia-nigra', 'cerebellar-hemisphere'] 
reg_phens = [(r,p) for p in phens for r in regns]

## paths 
base_path = '/data1/rubinov_lab/brain_genomics/scripts_twas'
expr_path = base_path + '/inputs_HCP/expr_regress'
phen_path = base_path + '/inputs_HCP/expr_regress'

demo_file = '/data1/rubinov_lab/brain_genomics/data_HCP/' + \
            'subject_demographics/cohort890_demographics.txt'
assc_path = base_path + '/outputs_HCP/assoc_1M'

out_file = base_path + '/outputs_HCP/multigene_linear_sets/test_eval_corrs.hdf5'

## scale expr and phen values  
expr_regr = {} ## k: reg, v: expr (genes, subjs) 
phen_regr = {} ## k: (reg, phen), v: phen (subjs,) 
for reg in regns: 
    rfile = '{}/{}.hdf5'.format(expr_path, reg)
    with h5py.File(rfile, 'r') as f: 
        evals = f['pred_expr'][()]
        emin = evals.min(axis=1)[:, None]
        emax = evals.max(axis=1)[:, None] 
        escl = (evals - emin) / (emax - emin)
        expr_regr[reg] = escl

for phen in phens: 
    pfile = '{}/{}.hdf5'.format(expr_path, phen) 
    with h5py.File(pfile, 'r') as f: 
        for reg in regns: 
            phevals = f[reg][()]
            phe_min = phevals.min() 
            phe_max = phevals.max()
            phe_scl = (phevals - phe_min) / (phe_max - phe_min)
            phen_regr[(reg, phen)] = phe_scl

## split cohort into train and test sets by family 
## generate a random split for each iteration 
subj_fids = np.loadtxt(demo_file, \
            delimiter='\t', \
            skiprows=1, \
            usecols=[2], \
            dtype=str)

fids = np.unique(subj_fids)
train_masks = {} 
testt_masks = {} 

ITR = 10000
for i in range(ITR): 
    train_fids, testt_fids = train_test_split(fids, \
                             test_size=0.25, \
                             random_state=i*150, \
                             shuffle=True)
    train_mask = np.isin(subj_fids, train_fids)
    testt_mask = ~train_mask

    train_masks[i] = train_mask 
    testt_masks[i] = testt_mask 

## gather observed corrs and pvals 
assocs = {} ## k: (reg, phen), v: [corrs, pvals]
for (reg, phen) in reg_phens: 
    afile = '{}/pvals_{}/{}.hdf5'.format(assc_path, phen, reg) 
    with h5py.File(afile, 'r') as f: 
        data = f['pearson'][()][:,:-1]
        #data[:,0] = np.where(np.isnan(data[:,0]), 0, data[:,0])
        #data[:,1] = np.where(np.isnan(data[:,1]), 100, data[:,1])
    assocs[(reg, phen)] = data 

## function: call the best-gene and gene-set models for a given reg phen 
def call_models(reg, phen, itr): 

    train_mask = train_masks[itr]
    testt_mask = testt_masks[itr]

    ## select gene and gene sets for models 
    midx = np.argmin(assocs[(reg, phen)][:,1])
    mask = (assocs[(reg, phen)][:,1] <= 0.05)

    ## split train and test sets 
    s_train = expr_regr[reg][midx][train_mask].reshape(-1,1)
    s_testt = expr_regr[reg][midx][testt_mask].reshape(-1,1)
    m_train = expr_regr[reg][mask][:,train_mask].T
    m_testt = expr_regr[reg][mask][:,testt_mask].T

    y_train = phen_regr[(reg, phen)][train_mask].reshape(-1,1)
    y_testt = phen_regr[(reg, phen)][testt_mask].reshape(-1,1)
    
    ## run regression without CV
    s_model = LinearRegression().fit(s_train, y_train)
    s_preds = s_model.predict(s_testt) 
    [s_rho, s_pvl] = pearsonr(y_testt[:,0], s_preds[:,0])

    m_model = LinearRegression().fit(m_train, y_train)
    m_preds = m_model.predict(m_testt) 
    [m_rho, m_pvl] = pearsonr(y_testt[:,0], m_preds[:,0])

    ## compare multi-gene weights to their single-gene correlations 
    ## TODO

    return [s_rho, s_pvl, m_rho, m_pvl]

## run all regression models 
all_corrs = np.zeros((5, 10, ITR, 2)) ## (phen, reg, itr, single/multi) 

for p, phen in enumerate(phens): 
    for r, reg in enumerate(regns):
        for itr in range(ITR):
            
            ## run regression 
            res = call_models(reg, phen, itr) 

            all_corrs[p][r][itr][0] = res[0]
            all_corrs[p][r][itr][1] = res[2]

        print(phen, reg)

## save results 
dims = 'phenotype, regions, iterations, (best-gene, gene-set)'
with h5py.File(out_file, 'w') as f: 
    f['pred-vs-expect-pearson'] = all_corrs
    f['dimensions'] = dims.encode()
    f['phenotype-order'] = phns.astype(bytes) 
    f['region-order'] = regs.astype(bytes)

