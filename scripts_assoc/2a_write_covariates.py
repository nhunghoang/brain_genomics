'''
Script for parsing HCP/UKB demographic files to 
generate a covariate file in the format needed 
for PrediXcanAssociation. 

Continuous variables are left as is, 
categorical variables are one-hot encoded. 

- Nhung (Sept 2021), updated by Tim (March 2022)  
'''

import csv
import numpy as np 
import h5py
import sys 

dset = sys.argv[1] ## HCP or UKB 

## paths 
if dset == 'HCP':
    dem_file = '/data1/rubinov_lab/brain_genomics/data_HCP/subject_demographics/COMPLETE_DEMOGRAPHICS.csv'
    gen_file = '/data1/rubinov_lab/brain_genomics/data_HCP/subject_demographics/HCP_unrestricted_4_23_2019_21_59_10.csv'
    sub_file = '/data1/rubinov_lab/brain_genomics/analyses_HCP/DATA_OUTPUT/subj_samp_assoc_order.hdf5' 
if dset == 'UKB': 
    dem_file = '/data1/rubinov_lab/brain_genomics/data_UKB/subject_demographics/demo_UKB.txt'

pca_path = '/data1/rubinov_lab/brain_genomics/analyses_{}/DATA_OUTPUT/eigen_results'.format(dset)
cov_file = '/data1/rubinov_lab/brain_genomics/analyses_{}/DATA_OUTPUT/covariates.txt'.format(dset) 

## covariates dictionary 
covs = ['race', 'age', 'gender', 'PC1', 'PC2']
cdict = {c:[] for c in covs}

###################################################################################################

## get subject IDs 
if dset == 'HCP':
    with h5py.File(sub_file, 'r') as f: 
        subjects = np.array(f['subjects']).astype(str) ## (890,)  

if dset == 'UKB': 
    subjects = np.loadtxt(dem_file, \
                          delimiter='\t', \
                          skiprows=1, \
                          dtype=str, \
                          usecols=[1], \
                          dtype=str) ## (2773,) 

## parse race, age, and gender 
if dset == 'HCP': 
    d_race = {}; d_agee = {}; d_gend = {} 
    with open(dem_file) as csvf:
        reader = csv.DictReader(csvf) 
        for row in reader: 
            subj = row['Subject'] 
            d_race[subj] = row['Race']
            d_agee[subj] = row['Age_in_Yrs']
    with open(gen_file) as csvf: 
        reader = csv.DictReader(csvf) 
        for row in reader: 
            subj = row['Subject'] 
            d_gend[subj] = row['Gender']

    races = np.array([d_race[s] for s in subjects])
    agess = np.array([d_agee[s] for s in subjects])
    gends = np.array([d_gend[s] for s in subjects])

if dset == 'UKB':
    races = np.loadtxt(dem_file, delimiter='\t', skiprows=1, dtype=str, usecols=[5])
    agees = np.loadtxt(dem_file, delimiter='\t', skiprows=1, dtype=str, usecols=[4])
    gends = np.loadtxt(dem_file, delimiter='\t', skiprows=1, dtype=str, usecols=[3])

## parse PCA 
pcs = np.loadtxt('{}/{}_cohort.pca'.format(pca_path, dset), skiprows=4)
pc1 = pcs[:,0]
pc2 = pcs[:,1]

###################################################################################################

## header: raceNA isNative isAsian isAfrican isEuropean isMulti isFemale age PC1 PC2
if dset == 'HCP':
    cols = ['raceNA', 'isNative', 'isAsian', 'isAfrican', 'isEuropean',\
            'isMulti', 'isFemale', 'age', 'PC1', 'PC2']
    header = '\t'.join(cols)

    race_num = {'Am. Indian/Alaskan Nat.': 1, \
                'Asian/Nat. Hawaiian/Othr Pacific Is.': 2, \
                'Black or African Am.': 3, \
                'More than one': 5, \
                'Unknown or Not Reported': 0, \
                'White': 4}

## header: raceNA isChinese isSouth_Asian isBlack isWhite isMulti isFemale age PC1 PC2
if dset == 'UKB':
    cols = ['raceNA', 'isChinese', 'isSouth_Asian', 'isBlack', 'isWhite',\
            'isMulti', 'isFemale', 'age', 'PC1', 'PC2']
    header = '\t'.join(cols)

    race_num = {'Chinese': 1, \
                'South_Asian': 2, \
                'Black': 3, \
                'Mixed': 5, \
                'Unknown': 0, \
                'White': 4}

## write covariates to file  
with open(cov_file, 'w') as f: 
    f.write(header + '\n')

    for s in range(subjects.size):
        data = np.zeros(len(cols), dtype=int).astype(str)
    
        ## one hot encoding for race 
        r = race_num[races[s]]
        data[r] = '1' 

        ## one hot encoding for gender 
        if gends[s] == 'F': data[-4] = '1' 

        ## age and PCs
        data[-3] = agees[s] 
        data[-2] = pc1[s] 
        data[-1] = pc2[s]

        ## save subject line 
        line = '\t'.join(data)
        f.write(line + '\n')  

