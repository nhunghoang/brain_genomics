'''
Search for clinical traits from PrediXVU that associate 
with selected genes per phenotype. 

This script uses multi-threading. Any PrediXVU gene 
filenames that can't be found are noted in 'need_run' 
so they can be identified in a different way (search_predixvu2.py). 

- Nhung, updated April 2022 
'''

import os 
import numpy as np 
import pandas as pd 
import subprocess 
from statsmodels.stats.multitest import fdrcorrection as FDR
from multiprocessing import Pool

## paths 
pdx_path = '/dors/capra_lab/data/predixVU/byGene' 
gene_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/predixvu/symbol_ensembl_region' 
brain_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/predixvu/brain_phecode_defns.csv'

out_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/predixvu/phecodes_genes' 
if not os.path.exists(out_path): os.mkdir(out_path) 

except_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/predixvu/need_run' 
if not os.path.exists(except_path): os.mkdir(except_path) 

## run once: 
## map phecodes to phenotypes (only for brain-related phenotypes)  
brain_phens = ['mental disorders', 'neurological'] 
#pfile = '/dors/capra_lab/data/predixVU/phecode_definitions1.2.tsv' 
#df = pd.read_csv(pfile, sep='\t')
#df = df.loc[df['category'].isin(brain_phens)]
#df.to_csv(brain_path)

## load phecode-phenotype mappings (clinical phenotypes)  
df = pd.read_csv(brain_path, sep=',')
phen_dict = df.set_index('phecode').to_dict()['phenotype'] ## k: float(phecode), v: clinical phenotype
brain_phecodes = list(phen_dict.keys())

## function: query phecodes that associate 
## with genes for the neural phenotype  
def query_predixvu(phen): 

    ## significant phecodes (across all regions) 
    phe_dict = {} ## k: phecode, v: [gene symbols]

    ## parse PredixVU gene filenames 
    with open('{}/{}.txt'.format(gene_path, phen), 'r') as f: 
        names = f.readlines() 
    prefixes = ['_'.join(n.split(' ')[:2]) for n in names] ## [symbol]_[ensembl]
    unique_prefixes, counts = np.unique(prefixes, return_counts=True)

    ## parse PredixVU gene file 
    for prefix in unique_prefixes: 
        fname = '{}/{}_predixVU.csv.gz'.format(pdx_path, prefix)

        ## try to read gene file
        if prefix.split('_')[0] == 'None': 
            print('n/a symbol ({})'.format(prefix)) 
            continue 
        try:
            df = pd.read_csv(fname) 
        except: 
            ## add gene to 'need run' list if filename not found 
            with open('{}/list_{}.txt'.format(except_path, phen), 'a') as f: 
                f.write(prefix + '\n') 
            print('n/a filename ({})'.format(prefix)) 
            continue  

        ## query significant PrediXVU associations 
        ## consider brain-related phecodes only 
        df = df[df['gene'].notnull()] 
        pcodes = df['phecode'].str[1:].astype(float) ## example format: X008.52
        df = df.loc[pcodes.isin(brain_phecodes)] 
        
        ## FDR correction 
        pvals0 = df['pvalue']
        _, pvals = FDR(pvals0)
        df = df.loc[pvals <= 0.05] 

        ## record phecode-gene associations  
        for index, row in df.iterrows(): 
            key = float(row['phecode'][1:]) ## example format: X008.52
            val = row['gene_name'] ## gene symbol 
            try: phe_dict[key].append(val)
            except KeyError: phe_dict[key] = [val]

    ## convert phecodes to clinical phenotypes 
    pkeys = list(phe_dict.keys())
    medtypes = [phen_dict[p] for p in pkeys]
    print('{}: {}'.format(phen, len(medtypes)))

    ## save to file 
    with open('{}/{}.txt'.format(out_path, phen), 'w') as f: 
        for pk in pkeys: 
            genes = ' '.join(np.unique(phe_dict[pk]))
            line = '{:>6}\t{}\t{}\n'.format(pk, phen_dict[pk], genes)
            f.write(line)

## neural phenotypes 
phens = ['alff', 'regional_homogeneity', 'gm_volume', \
        'connectivity_mean', 'connectivity_variance', \
        'falff', 'gradient', 'myelination', \
        'timeseries_variance'] #, 'fa', 'md']

## run the query 
pool = Pool(processes=9) 
pool.map(query_predixvu, phens) 
