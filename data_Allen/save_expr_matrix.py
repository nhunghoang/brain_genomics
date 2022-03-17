'''
Save (n_subjects * n_regions * n_probes) Allen matrix 
with nans as needed. 

- Nhung, March 2022 
'''

import numpy as np 
import os 
import h5py 

## input / output paths 
abi_path = '/data1/rubinov_lab/brain_genomics/data_Allen'
probe_file = abi_path + '/probes_as_entrez.txt'
region_file = '/data1/rubinov_lab/brain_genomics/data_HCP/hoacer_sn_hth/regions-115.txt' 
abbrev_file = '/data1/rubinov_lab/brain_genomics/data_HCP/hoacer_sn_hth/naming-115.txt'

symbols_file = '/data1/rubinov_lab/brain_genomics/analyses_HCP/predixcan_validation/compare_allen_hcp/2022_ensembl_entrez.txt'

out_file = abi_path + '/allen_expr_matrix.hdf5'

## matrix dimensions 
abi_subjs = np.array([i for i in os.listdir(abi_path) if i.split('_')[0] == 'normalized'])
abi_probes = np.loadtxt(probe_file, dtype=str)
abi_regions = np.loadtxt(region_file, delimiter='\n', dtype=str) 
region_idx = {reg:i for i,reg in enumerate(abi_regions)} 

## map biomaRt Entrez gene IDs to Ensembl IDs, skip if no mapping
with open(symbols_file, 'r') as f:
    f.readline()
    lines = [i.strip().split(',') for i in f.readlines()]
gene_dict = {j[0]:j[1] for j in lines if j[0] != ''} ## k: Entrez, v: Ensembl

## convert Allen genes from Entrez to Ensembl 
ensembl = [] 
for entrez in abi_probes: 
    try: ensembl.append(gene_dict[entrez].encode()) 
    except KeyError: ensembl.append('n/a'.encode()) 
abi_probes = np.array(ensembl) 

## initial NaN matrix 
data_mat = np.empty((abi_subjs.size, abi_regions.size, abi_probes.size)) 
data_mat[:] = np.nan 

for s,subj in enumerate(abi_subjs):

    ## load subject region dimension 
    subj_reg_path = '{}/{}/hoacer_sn_hth_regions.txt'.format(abi_path, subj)
    subj_regs = np.loadtxt(subj_reg_path, delimiter='\t', usecols=[1], dtype=str)

    ## gather index mask of unique regions 
    reg_masks = {} ## k: reg, v: (all sites) mask  
    for reg in np.unique(subj_regs): 
        reg_masks[reg] = np.in1d(subj_regs, [reg])
    
    ## load subject expression matrix 
    expr_path = '{}/{}/MicroarrayExpression.csv'.format(abi_path, subj)
    expr_all = np.loadtxt(expr_path, delimiter=',') ## (all probes * all sites)
    expr_all = expr_all[:,1:] ## first column is probe ID 

    ## normalize expression of each probe (across sites) 
    expr_mean = expr_all.mean(axis=1)[:,None]
    expr_stdv = expr_all.std(axis=1)[:,None] 
    expr_norm = (expr_all-expr_mean) / expr_stdv 
    #expr_norm = expr_all

    ## average site expression per unique region 
    for reg, mask in reg_masks.items():
        expr_reg = expr_norm[:, mask] ## (all probes * region sites) 
        expr_avg = expr_reg.mean(axis=1) ## (all probes) 
    
        ## add average regional expression to main matrix 
        reg_idx = region_idx[reg]
        data_mat[s][reg_idx] = expr_avg 

## map relevant atlas labels to PrediXcan abbreviations 
with open(abbrev_file, 'r') as f:
    f.readline()
    lines = [i.split('\t') for i in f.readlines()]
    reg_map = [(j[0], j[1]) for j in lines] 

abbrev_dict = {} ## k: atlas label (left hemisphere only), v: PrediXcan abbreviation
for (label,abbrev) in reg_map: 
    if label.split(' ')[0] != 'Right': 
        if abbrev.split('-')[-1] != 'LH': 
            abbrev_dict[label] = abbrev 
        else: 
            abbrev_dict[label] = abbrev[:-3] 

## change Allen atlas names to PrediXcan abbreviations if possible 
encoded_regs = [] 
for reg in abi_regions: 
    try: 
        abbrev = abbrev_dict[reg] 
        encoded_regs.append(abbrev.encode())  
    except KeyError: 
        encoded_regs.append(reg.encode()) 

## save main matrix and its dimensions 
with h5py.File(out_file, 'w') as f: 
    f['subjects'] = [s.split('_')[-1].encode() for s in abi_subjs] 
    f['probe_genes'] = abi_probes 
    f['regions'] = encoded_regs 
    f['expression'] = data_mat 

print('Allen matrix saved with dimensions {}'.format(data_mat.shape))

