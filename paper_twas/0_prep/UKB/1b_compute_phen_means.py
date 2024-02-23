'''
Get phenotypes for UKB subset of interest.

- Nhung, Feb 2024
'''

import numpy as np 
import pandas as pd 
import sys 

group = sys.argv[1] ## UKB, UKB/nonBrit, ... 

## paths
main_path = '/data1/rubinov_lab/brain_genomics'
pids_path = f'{main_path}/paper_twas/inputs_UKB/phen_download_info.csv'

vols_path = f'{main_path}/data_UKB/downloads/structural/data'
coho_path = f'{main_path}/paper_twas/inputs_{group}/cohort.txt' 

outs_path = f'{main_path}/paper_twas/inputs_{group}/phenotypes/vol_mean.csv'

## load phenotype IDs (from UKB downloads) 
pids = pd.read_csv(pids_path) 

## load cohort EIDs 
coho = pd.read_csv(coho_path, header=None, sep=' ')[0].values

## parse original phen files 
pdfs = [] 
for pf, info in pids.groupby('ukb_file'): 

    pfile = f'{vols_path}/{pf}.csv'
    pcols = info['ukb_field_id'].values 

    ## read phens
    pdf = pd.read_csv(pfile, index_col='eid')

    ## slice cohort 
    pdf = pdf.reindex(coho)

    ## slice volumes
    pdf = pdf[pcols] 
    pdfs.append(pdf)

## merge phens 
pdf = pd.concat(pdfs, axis=1)

## get mean of every L/R pair 
for ureg, info in pids.groupby('ukb_region'): 
    rcols = info['ukb_field_id']
    pdf[ureg] = pdf[rcols].mean(axis=1)

## get mean of every TWAS region 
for treg, info in pids.groupby('twas_region'):  
    uregs = info['ukb_region'].unique()

    pdf[treg] = pdf[uregs].mean(axis=1)

## save vol means 
tregs = pids['twas_region'].unique()
pdf.to_csv(outs_path, columns=tregs, sep='\t', index=True) 

