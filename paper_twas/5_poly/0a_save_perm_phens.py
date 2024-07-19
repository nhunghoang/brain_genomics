'''
'''

import numpy as np 
import pandas as pd 
import h5py 

from scipy.stats import pearsonr

## params 
group = 'HCP/nonTwin' #sys.argv[1] ## HCP, HCP/nonTwin, ...

regions = ['dlpfc', 'anterior-cingulate', 'amygdala', 'hippocampus', \
           'caudate', 'putamen', 'nucleus-accumbens', 'cerebellar-hemisphere']
phens = ['vol_mean', 'alff_mean', 'reho_noGS_mean', 'connmean_noGS_mean']
reg_phens = [(r,p) for r in regions for p in phens]

## paths 
main_path = '/data1/rubinov_lab/brain_genomics/paper_twas'
phen_path = f'{main_path}/inputs_{group}/phenotypes'

## load phen data
phen_all = {} ## k: (reg, phen), v: (perm, subj)

for phen in phens:
    ppath = f'{phen_path}/{phen}.csv'
    df = pd.read_table(ppath, sep='\t', usecols=regions)
    for reg in regions:
        phen_all[(reg,phen)] = df[reg].values

############################################################

## DO ONCE: save 10k permutations of the phens 

nperms = 5000
subjs = phen_all[(regions[0], phens[0])].shape[0]
ncols = len(regions) * nperms
cols = ['{}_{}'.format(reg, itr+5000) for reg in regions for itr in range(nperms)] 

for phen in phens: 
    data = np.zeros((subjs, ncols))
    accu = 0 
    for reg in regions:
        vals = phen_all[(reg,phen)]
        for itr in range(nperms): 
            data[:,accu] = np.random.permutation(vals)
            accu += 1

    df = pd.DataFrame(data, columns=cols)
    ofile = f'{phen_path}/perm_{phen}2.csv'
    df.to_csv(ofile, sep='\t', index=False)
    print(phen)

