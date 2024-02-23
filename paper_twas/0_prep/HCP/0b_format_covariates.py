'''
Script for parsing UKB demographic files to 
generate a covariate file in the format needed 
for PrediXcanAssociation. 

Continuous variables are left as is, 
categorical variables are one-hot encoded. 

Columns: subject sample gender age PC1 PC2 ...  

- Nhung (updated Feb 2024)  
'''

import numpy as np 
import pandas as pd 
import sys 

## paths 
main_path = '/data1/rubinov_lab/brain_genomics/paper_twas/inputs_HCP'
demo_path = f'{main_path}/demographics.csv'
eign_path = f'{main_path}/eigen_output/jti2.pca.evec'
coho_path = f'{main_path}/cohort.txt' 

outs_path = f'{main_path}/covariates.csv'

## load demographics 
demo = pd.read_csv(demo_path, usecols=['subject', 'sample', 'age', 'gender'])
demo['isMale'] = demo['gender'].apply(lambda x: 1 if x=='M' else 0) 
demo = demo.drop('gender', axis=1)

## load eigendata 
with open(eign_path, 'r') as f: 
    _ = f.readline() 
    lines = f.readlines() 

info = np.array([i.split() for i in lines])[:,:-1]
cols = ['sample'] + [f'PC{i}' for i in range(1,41)]
eign = pd.DataFrame(info, columns=cols)

## merge 
df = demo.merge(eign, on='sample', how='inner')
df = df.set_index('sample')

## reorder 
coho = pd.read_csv(coho_path, sep=' ', header=None)
coho.columns = ['subject', 'sample']
df = df.reindex(coho['sample'].values)

## save 
df.to_csv(outs_path, sep='\t')
