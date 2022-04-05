'''
Search for clinical traits from PrediXVU that associate
with selected genes per phenotype *permutation*.

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
gene_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/predixvu/null/symbol_ensembl_region' 
brain_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/predixvu/brain_phecode_defns.csv'

out_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/predixvu/null/phecodes_genes'
if not os.path.exists(out_path): os.mkdir(out_path)

except_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/predixvu/null/need_run'
if not os.path.exists(except_path): os.mkdir(except_path)

## read phecode-phenotype mapping 
brain_phens = ['mental disorders', 'neurological'] 
df = pd.read_csv(brain_path, sep=',')
phen_dict = df.set_index('phecode').to_dict()['phenotype'] ## k: float(phecode), v: phenotype
brain_phecodes = list(phen_dict.keys())

## function: query phecodes that associate
## with genes for the neural phenotype permutation 
def query_predixvu(data):
    phen = data['phen']
    perm = data['perm'] 

    ## significant phecodes (across all regions)
    phe_dict = {} ## k: phecode, v: [gene symbols]

    ## parse PredixVU gene filenames
    with open('{}/{}_{}.txt'.format(gene_path, phen, perm), 'r') as f:
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
            with open('{}/list_{}_{}.txt'.format(except_path, phen, perm), 'a') as f:
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
    print('{}-{}: {}'.format(phen, perm, len(medtypes)))

    ## save to file
    with open('{}/{}_{}.txt'.format(out_path, phen, perm), 'w') as f:
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
phens_perms = [{'phen':phen, 'perm':i} for phen in phens for i in range(100)]
pool = Pool(processes=60)
pool.map(query_predixvu, phens_perms) 
