'''
Per phenotype permutation, gather top N assocs based 
on pvals and the number of observed instances. 

- Nhung, March 2022 
'''

import os 
import numpy as np 
import requests 
import json 

phens = ['alff', 'regional_homogeneity', 'gm_volume', \
        'connectivity_mean', 'connectivity_variance', \
        'falff', 'gradient', 'myelination', \
        'timeseries_variance']#, 'fa', 'md']
regs = os.listdir('/data1/rubinov_lab/brain_genomics/analyses_HCP/assoc/null_pvals_alff') 

## paths 
obsv_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/predixvu/symbol_ensembl_region'
null_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/predixvu/top_null/sorted_pvals'
out_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/predixvu/top_null/symbol_ensembl_region'
if not os.path.exists(out_path): os.mkdir(out_path) 

for phen in phens: 

    ## count number of observed instances 
    with open('{}/{}.txt'.format(obsv_path, phen), 'r') as f: 
        obsv_count = len(f.readlines())

    for perm in range(100): 

        lines = [] 

        ## gather the top N null assocs based on observed count 
        path = '{}/{}_{}.txt'.format(null_path, phen, perm)
        data = np.loadtxt(path, delimiter='\t', dtype=str) 
        data = data[:obsv_count] 
        regs = data[:,2] 
        sel_genes = data[:,0] 

        ## get gene symbols 
        url = 'https://biotools.fr/human/ensembl_symbol_converter/' 
        json_genes = json.dumps(sel_genes.tolist())          
        in_data = {'api':1, 'ids':json_genes} 
        out_data = requests.post(url, data=in_data) 
        symbols = json.loads(out_data.text) ## {ensembl:symbol} 

        if len(symbols) == 0: continue 

        ## format lines 
        for reg, ens in zip(regs, sel_genes): 
            line = '{} {} {}\n'.format(symbols[ens], ens, reg)
            lines.append(line) 

        ## write to file 
        with open('{}/{}_{}.txt'.format(out_path, phen, perm), 'w') as f: 
            for line in lines: 
                f.write(line) 

    print('done with {}'.format(phen))
        

