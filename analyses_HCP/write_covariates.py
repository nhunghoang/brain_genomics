'''
Script for parsing HCP demographic files to 
generate a covariate file in the format needed 
for PrediXcanAssociation. 

Continuous variables are left as is, 
categorical variables are one-hot encoded. 

- Nhung (Sept 2021) 
'''

import csv
import numpy as np 
import h5py
import sys 

## paths 
dem_file = '/data1/rubinov_lab/brain_genomics/data_HCP/subject_demographics/COMPLETE_DEMOGRAPHICS.csv'
gen_file = '/data1/rubinov_lab/brain_genomics/data_HCP/subject_demographics/HCP_unrestricted_4_23_2019_21_59_10.csv'
grp_file = '/data1/rubinov_lab/brain_genomics/analyses_HCP/DATA_OUTPUT/subj_samp_assoc_order.hdf5' 

spl_file = '/data1/rubinov_lab/brain_genomics/analyses_HCP/DATA_OUTPUT/train_test_assoc_split.hdf5'
pca_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/DATA_OUTPUT/eigen_results'
cov_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/DATA_OUTPUT'

## covariates dictionary 
covs = ['race', 'age', 'gender', 'PC1', 'PC2']
cdict = {c:[] for c in covs}

## train/test split 
with h5py.File(spl_file, 'r') as f: 
    train_idx = np.array(f['train_idx_890'])
    testt_idx = np.array(f['test_idx_890'])
with h5py.File(grp_file, 'r') as f: 
    all_subjs = np.array(f['subjects']).astype(str)
subj_dict = {'train':all_subjs[train_idx], 'test':all_subjs[testt_idx]}

## parse demographic files (which contains info for all HCP subjects)  
race_all = {}; age_all = {}; gender_all = {} ## k: subject ID, v: value as text 
with open(dem_file) as csvf: 
    reader = csv.reader(csvf, delimiter=',') 
    for line in reader: 
        subj = line[0]
        race_all[subj] = line[9]   
        age_all[subj] = line[1] 

with open(gen_file) as csvf: 
    reader = csv.reader(csvf, delimiter=',') 
    for line in reader: 
        subj = line[0]
        gender_all[subj] = line[3]

## parse PCA files 
train_pcs = np.loadtxt('{}/train_cohort.pca'.format(pca_path), skiprows=4)
testt_pcs = np.loadtxt('{}/test_cohort.pca'.format(pca_path), skiprows=4)
pc_dict = {'train':train_pcs, 'test':testt_pcs}

## write covariates to file  
## header: raceNA isNative isAsian isAfrican isEuropean isMulti isFemale age PC1 PC2
cols = ['raceNA', 'isNative', 'isAsian', 'isAfrican', 'isEuropean',\
        'isMulti', 'isFemale', 'age', 'PC1', 'PC2']
header = '\t'.join(cols)

race_num = {'Am. Indian/Alaskan Nat.': 1, \
            'Asian/Nat. Hawaiian/Othr Pacific Is.': 2, \
            'Black or African Am.': 3, \
            'More than one': 5, \
            'Unknown or Not Reported': 0, \
            'White': 4}

## loop through train/test 
for grp in ['train', 'test']:
    cov_file = '{}/{}_covariates.txt'.format(cov_path, grp)
    subjects = subj_dict[grp]
    pcs = pc_dict[grp]
    pc1 = pcs[:,0]
    pc2 = pcs[:,1]

    with open(cov_file, 'w') as f: 
        f.write(header + '\n')

        for i,s in enumerate(subjects): 
            data = np.zeros(len(cols), dtype=int).astype(str)
        
            ## one hot encoding for race 
            r = race_num[race_all[s]]
            data[r] = '1' 

            ## one hot encoding for gender 
            if gender_all[s] == 'F': data[6] = '1' 

            ## age and PCs
            data[7] = age_all[s] 
            data[8] = pc1[i] 
            data[9] = pc2[i]

            ## save subject line 
            line = '\t'.join(data)
            f.write(line + '\n')  

