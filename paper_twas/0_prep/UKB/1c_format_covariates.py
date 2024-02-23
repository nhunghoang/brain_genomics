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

group = sys.argv[1] ## UKB, UKB/nonBrit, ...

## paths 
main_path = '/data1/rubinov_lab/brain_genomics/'
demo_path = f'{main_path}/data_UKB/downloads/all_covariates.csv'
coho_path = f'{main_path}/paper_twas/inputs_{group}/cohort.txt' 

outs_path = f'{main_path}/paper_twas/inputs_{group}/covariates.csv'

## format covariate names
cols_old = ['21003-2.0', '22001-0.0'] + ['22009-0.{}'.format(i) for i in range(1,41)]
cols_new = ['age', 'isMale']  + ['PC{}'.format(i) for i in range(1,41)]
cdict = {co: cn for co, cn in zip(cols_old, cols_new)}

## read demographics data 
df = pd.read_csv(demo_path, index_col='eid')
df = df[cols_old].rename(columns=cdict)

## load cohort data and slice demographics  
coho = pd.read_csv(coho_path, header=None, sep=' ')[0]
covs = df.reindex(coho)

## save 
covs.to_csv(outs_path, sep='\t')
