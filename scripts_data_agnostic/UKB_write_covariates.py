'''
Script for parsing UKB demographic files to 
generate a covariate file in the format needed 
for PrediXcanAssociation. 

Continuous variables are left as is, 
categorical variables are one-hot encoded. 

- Tim, Mar 13 2022
'''

import csv
import numpy as np 
import h5py
import sys 

## paths 
dem_file = '/data1/rubinov_lab/brain_genomics/data_UKB/subject_demographics/demo_UKB.txt'
pca_path = '/data1/rubinov_lab/brain_genomics/analyses_UKB/DATA_OUTPUT/eigen_results'
cov_file = '/data1/rubinov_lab/brain_genomics/analyses_UKB/DATA_OUTPUT/covariates.txt'

## covariates dictionary 
covs = ['race', 'age', 'gender', 'PC1', 'PC2']
cdict = {c:[] for c in covs}


## get subject IDs for UKB 2773 cohort
subjects = np.loadtxt(dem_file, delimiter='\t', skiprows=1, dtype=str, usecols=[1])
race_all = np.loadtxt(dem_file, delimiter='\t', skiprows=1, dtype=str, usecols=[5])
age_all = np.loadtxt(dem_file, delimiter='\t', skiprows=1, dtype=str, usecols=[4])
gender_all = np.loadtxt(dem_file, delimiter='\t', skiprows=1, dtype=str, usecols=[3])

## parse PCA file 
pcs = np.loadtxt('{}/UKB_cohort.pca'.format(pca_path), skiprows=4)
pc1 = pcs[:,0]
pc2 = pcs[:,1]

## write covariates to file  
## header: raceNA isChinese isSouth_Asian isBlack isWhite isMulti isFemale age PC1 PC2
cols = ['raceNA', 'isChinese', 'isSouth_Asian', 'isBlack', 'isWhite',\
        'isMulti', 'isFemale', 'age', 'PC1', 'PC2']
header = '\t'.join(cols)

race_num = {'Chinese': 1, \
            'South_Asian': 2, \
            'Black': 3, \
            'Mixed': 5, \
            'Unknown': 0, \
            'White': 4}

with open(cov_file, 'w') as f: 
    f.write(header + '\n')

    for i in range(subjects.shape[0]): 
        data = np.zeros(len(cols), dtype=int).astype(str)
    
        ## one hot encoding for race 
        r = race_num[race_all[i]]
        data[r] = '1' 

        ## one hot encoding for gender 
        if gender_all[i] == 'F': data[6] = '1' 

        ## age and PCs
        data[7] = age_all[i] 
        data[8] = pc1[i] 
        data[9] = pc2[i]

        ## save subject line 
        line = '\t'.join(data)
        f.write(line + '\n')  

