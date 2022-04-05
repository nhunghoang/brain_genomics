'''
Prepare gene data for PrediXVU query. Gather all regional genes 
per phenotype *permutation* and store alongside gene symbols.

- Nhung, updated April 2022 
'''

import os 
import numpy as np 
import requests 
import json 
from multiprocessing import Pool

## paths 
in_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/assoc' ## null_pvals_{phen}/{region}/{perm}.txt 
out_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/predixvu/null/symbol_ensembl_region'
if not os.path.exists(out_path): os.mkdir(out_path) 

## phenotypes and regions 
phens = ['alff', 'regional_homogeneity', 'gm_volume', \
        'connectivity_mean', 'connectivity_variance', \
        'falff', 'gradient', 'myelination', \
        'timeseries_variance'] #, 'fa', 'md']
regs = os.listdir('/data1/rubinov_lab/brain_genomics/analyses_HCP/assoc/null_pvals_alff') 

## function: gather regional genes by phenotype permutation 
def gather_null_genes(data):  
    phen = data['phen'] 
    perm = data['perm'] 
    lines = [] 
    for reg in regs: 

        ## single genes (FDR <= 0.05) 
        data = np.loadtxt('{}/null_pvals_{}/{}/{}.txt'.format(in_path, phen, reg, perm), \
            dtype=str, delimiter='\t', skiprows=1, usecols=[0,3]) 
        genes = np.array([g.split('.')[0] for g in data[:,0]])
        pvals = data[:,1].astype(float) 
        sel_genes = genes[pvals <= 0.05] 

        ## get gene symbols 
        url = 'https://biotools.fr/human/ensembl_symbol_converter/' 
        json_genes = json.dumps(sel_genes.tolist())          
        in_data = {'api':1, 'ids':json_genes} 
        out_data = requests.post(url, data=in_data) 
        symbols = json.loads(out_data.text) ## {ensembl:symbol} 

        if len(symbols) == 0: continue 

        ## format lines 
        for ens,sym in symbols.items(): 
            line = '{} {} {}\n'.format(sym, ens, reg) 
            lines.append(line) 

    ## write to file 
    with open('{}/{}_{}.txt'.format(out_path, phen, perm), 'w') as f: 
        for line in lines: 
            f.write(line) 
    
    print('{}-{}: {}'.format(phen, perm, len(lines)))

phen_perms = [{'phen':phen, 'perm':perm} for phen in phens for perm in range(100)]
pool = Pool(processes=60)
pool.map(gather_null_genes, phen_perms) 
