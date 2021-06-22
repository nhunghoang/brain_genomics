'''
Parse the PrediXcan expression imputation logs, and filter out all gene 
models that HCP did not have all the snps for. Record those genes in a 
collective file. Per tissue, save the list of valid genes.  

- Nhung, May 2021 
'''

import os 
import numpy as np 
import h5py 

predix_dir = '/data1/rubinov_lab/brain_genomics/data_HCP/expression'
filter_dir = '/data1/rubinov_lab/brain_genomics/data_HCP/expression/filtered_completeness'
gmodel_dir = '/data1/rubinov_lab/brain_genomics/data_HCP/predixcan_v8/gene_models_completeness'

no_keep_file = '/data1/rubinov_lab/brain_genomics/data_HCP/expression/not_kept_model_complete.txt'
no_keep = open(no_keep_file, 'w') 

for log in os.listdir(predix_dir): 
    if log[-4:] != '.log': continue 
    tissue = log.split('.')[0] 

    genes_to_keep = [] 
    index_to_keep = [] 
    with open('{}/{}'.format(predix_dir, log), 'r') as f: 
        f.readline() ## header 
        lines = f.readlines() 

    ## compare num snps in gene model to num snps used from HCP 
    for i,line in enumerate(lines): 
        [gene, gene_name, n_snps_in_model, n_snps_used] = line.split('\t')[:4]
        if n_snps_used == 'NA':
            no_keep_line = '\t'.join([tissue, gene, n_snps_in_model, n_snps_used])
            no_keep.write(no_keep_line + '\n') 
        elif float(n_snps_in_model) != float(n_snps_used):
            no_keep_line = '\t'.join([tissue, gene, n_snps_in_model, n_snps_used])
            no_keep.write(no_keep_line + '\n') 
        else:
            genes_to_keep.append(gene) 
            index_to_keep.append(i)  

    ## save tissue-specific file of gene models kept 
    with open('{}/{}.txt'.format(gmodel_dir, tissue), 'w') as f1: 
        for gene in genes_to_keep: 
            f1.write(gene + '\n') 

    ## write new expression hdf5 files with just the valid genes 
    expr_data = {} 
    with h5py.File('{}/{}.hdf5'.format(predix_dir, tissue), 'r') as f2: 
        for key in f2.keys(): 
            expr_data[key] = np.array(f2[key])
    with h5py.File('{}/{}.hdf5'.format(filter_dir, tissue), 'w') as f3: 
        f3['genes'] = expr_data['genes'][index_to_keep]
        f3['pred_expr'] = expr_data['pred_expr'][index_to_keep]
        f3['samples'] = expr_data['samples']

    print('{:22}: {} / {} genes kept'.format(tissue, len(genes_to_keep), len(lines)))

no_keep.close()
        
