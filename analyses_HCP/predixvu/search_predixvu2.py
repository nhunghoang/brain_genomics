'''
Search for clinical traits from PrediXVU that associate 
with selected genes per phenotype. 

This script is for genes whose expected filenames weren't 
found (i.e. gene symbol and ensembl don't match). 
 
If another file is found for the gene, update the 
corresponding results. 

- Nhung, updated April 2022 
'''

import os 
import numpy as np 
import pandas as pd 
import subprocess 
from statsmodels.stats.multitest import fdrcorrection as FDR

## paths 
pdx_path = '/dors/capra_lab/data/predixVU/byGene' 
gene_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/predixvu/symbol_ensembl_region' 
brain_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/predixvu/brain_phecode_defns.csv'

out_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/predixvu/phecodes_genes' 
except_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/predixvu/need_run' 

## load phecode-phenotype mappings (clinical phenotypes)  
brain_phens = ['mental disorders', 'neurological'] 
df = pd.read_csv(brain_path, sep=',')
phen_dict = df.set_index('phecode').to_dict()['phenotype'] ## k: float(phecode), v: clinical phenotype
brain_phecodes = list(phen_dict.keys())

## neural phenotypes
phens = ['alff', 'regional_homogeneity', 'gm_volume', \
        'connectivity_mean', 'connectivity_variance', \
        'falff', 'gradient', 'myelination', \
        'timeseries_variance'] #, 'fa', 'md']

## loop through phenotypes 
for phen in phens: 
    phe_dict = {} ## k: phecode, v: [gene symbols]

    ## gather genes to query 
    plist = '{}/list_{}.txt'.format(except_path, phen) 
    if not os.path.exists(plist): continue 

    with open(plist, 'r') as f: 
        prefixes = [i.strip() for i in f.readlines()]  

    for prefix in prefixes: 
        fname = '{}/{}_predixVU.csv.gz'.format(pdx_path, prefix) 

        ## try to find PrediXVU gene file
        cmd = 'ls {}/*{}*'.format(pdx_path, prefix.split('_')[1])
        try: sp_out = subprocess.check_output(cmd, shell=True, stderr=subprocess.DEVNULL)
        except: print('DNE: {}'.format(prefix.split('_')[1])); continue
        fname = sp_out.decode().strip()
        df = pd.read_csv(fname)

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

    ## parse existing results 
    with open('{}/{}.txt'.format(out_path, phen), 'r') as f: 
        old_lines = f.readlines()

    results = {} ## k: (phecode, medphen), v: [gene symbols]
    for line in old_lines:
        [pcode, mphen, genes] = line.strip().split('\t')
        results[(pcode, mphen)] = genes.split(' ')

    ## add new results 
    old_keys = results.keys() 
    for pcode, genes in phe_dict.items(): 
        key = (pcode, phen_dict[pcode]) 
        if key in old_keys: 
            results[key].extend(genes) 
            print('{}: added gene'.format(phen))
        else: 
            results[key] = genes 
            print('{}: added new med phen'.format(phen))

    ## rewrite file with updated results 
    with open('{}/{}.txt'.format(out_path, phen), 'w') as f: 
        for (pcode, mphen), genes in results.items(): 
            genes_str = ' '.join(np.unique(genes))
            line = '{:>6}\t{}\t{}\n'.format(pcode, mphen, genes_str)
            f.write(line)
