'''
Identify clinical traits from PrediXVU that
associate with our neuroimaging genes of
interest, per set of phenotype genes.

Per clinical trait:

- (PrediXVU table look-up)
  Record and count the neuroimaging genes that 
  are associated with phecodes in PrediXVU. 

- (Enrichment analysis)
  Count the number of observed neuroimaging genes 
  associated with the trait and compare to the 
  average number of corresponding permutated genes.

The FDR correction is computed over all 156 
phecodes available and all phenotypes. 

- Nhung, March 2023

'''

import numpy as np 
import h5py
import pandas as pd 
from statsmodels.stats.multitest import fdrcorrection as FDR

num_nulls = 10000
phenotypes = ['gm_volume', 'alff', 'myelination', 'reho_noGS', 'connmean_noGS'] 
regions = ['caudate', 'nucleus-accumbens', 'hippocampus', 'amygdala', 'anterior-cingulate', \
    'putamen', 'hypothalamus', 'substantia-nigra', 'frontal-pole', 'cerebellar-hemisphere']

## paths 
main_path = '/data1/rubinov_lab/brain_genomics'  

table_file = main_path + '/models_PrediXcan_v8/predixvu_assocs.csv'
assoc_path = main_path + '/scripts_assoc_clean/outputs_HCP/assoc_1M' 

out_path = main_path + '/scripts_assoc_clean/outputs_HCP/predixvu' 

## read PrediXVU table 
## ensembl, phecode, p-value, fdr-val, effsize, tissue 
df = pd.read_csv(table_file, \
     usecols=['ensembl', 'phecode', 'p-value', 'fdr-val'])

## select PrediXVU gene-phecode assocs that pass 
## threshold, ignoring tissue-specificity of genes
df = df.loc[df['p-value'] <= 0.01]
df = df[['ensembl', 'phecode']].drop_duplicates()
phecodes = df['phecode'].drop_duplicates()

## gather selected genes per phenotype  
obsv_genes = {} ## k: phen, v: gene array 
null_genes = {} ## k: phen, v: list of gene arrays 
for phen in phenotypes: 
    obsv_set = [] 
    null_set = [[] for _ in range(num_nulls)] 

    for reg in regions: 

        ## load genes and select observed data 
        ofile = '{}/pvals_{}/{}.hdf5'.format(assoc_path, phen, reg)
        with h5py.File(ofile, 'r') as f: 
            reg_genes = f['genes'][()].astype(str) 
            obsv_pvals = f['pearson'][()][:,1]
        obsv_set.extend(reg_genes[obsv_pvals <= 0.005])

        ## select null data 
        nfile = '{}/nulls/pvals_{}/{}.hdf5'.format(assoc_path, phen, reg)
        with h5py.File(nfile, 'r') as f: 
            null_pvals = f['null_pearson'][()][:,:,1] ## (perms, genes) 

        mask = (null_pvals <= 0.005)
        for i in range(num_nulls): 
            null_set[i].extend(reg_genes[mask[i]])            

    obsv_genes[phen] = np.unique(obsv_set)
    null_genes[phen] = [np.unique(n) for n in null_set]

## count/record genes per phenotype group of phecodes 
obsv_counts = {phen: {pc:0 for pc in phecodes} for phen in phenotypes}
null_counts = {phen: {pc:np.zeros(num_nulls) for pc in phecodes} for phen in phenotypes}

cloud_genes = {phen: {pc:[] for pc in phecodes} for phen in phenotypes} 

for phen in phenotypes: 
    obsv_set = obsv_genes[phen]
    pdx = df.loc[df['ensembl'].isin(obsv_set)] 
    counts = pdx['phecode'].value_counts() 
    for phecode in phecodes: 
        obsv_counts[phen][phecode] = counts[phecode]

        pgenes = pdx.loc[pdx['phecode']==phecode]['ensembl']
        cloud_genes[phen][phecode] = pgenes

    ## TODO: double check  
    null_set = null_genes[phen]
    for n, nset in enumerate(null_set):
        pdx = df.loc[df['ensembl'].isin(nset)] 
        counts = pdx['phecode'].value_counts() 
        for phecode in phecodes: 
            null_counts[phen][phecode][n] = counts[phecode]

## compute phecode enrichment p-values 
pdata = np.zeros((len(phenotypes), phecodes.size)) 
for i, phen in enumerate(phenotypes): 
    for j, code in enumerate(phecodes): 
        pdata[i][j] = (null_counts[phen][code] >= obsv_counts[phen][code]).mean() 

_, fdata = FDR(pdata) 

## save results 
with h5py.File(out_path + '/cloud_results.hdf5', 'w') as f: 
    for phen, codes in cloud_genes.items(): 
        for code, glist in codes.items(): 
            key = '{}-{}'.format(code, phen)
            f[key] = glist.to_numpy(dtype-bytes)

with h5py.File(out_path + '/set_enrichments.hdf5', 'w') as f: 
    f['phenotypes'] = np.array(phenotypes).astype(bytes) 
    f['phecodes'] = phecodes.to_numpy(dtype=bytes)
    f['p-values'] = pdata 
    f['fdr-values'] = fdata
     
