'''
For each phenotype, gather all regional instances. 
Sort FDR p-value of these regional phenotypes from lowest to highest. 
Save phenotype table as gene, FDR, region. 

- Nhung, March 2022 
'''

import numpy as np 
import os 

## paths  
pval_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/assoc' ## null_pvals_{phen}/{region}/{perm}.txt 
out_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/predixvu/top_null/sorted_pvals'
if not os.path.exists(out_path): os.mkdir(out_path) 

## regions and phenotypes 
regs = os.listdir('/data1/rubinov_lab/brain_genomics/analyses_HCP/assoc/null_pvals_alff')
phens = ['alff', 'regional_homogeneity', 'gm_volume', \
        'connectivity_mean', 'connectivity_variance', \
        'falff', 'gradient', 'myelination', \
        'timeseries_variance']

## loop through phenotype permutations 
for perm in range(100): 
    for phen in phens: 
        phen_data = [] ## np.array([gene, FDR, region])  

        for reg in regs: 
            ## gather assoc data for region-specific genes 
            reg_data = np.loadtxt('{}/null_pvals_{}/{}/{}.txt'\
                .format(pval_path, phen, reg, perm), \
                dtype=str, delimiter='\t', skiprows=1, usecols=[0,3])

            ## append all regional data to phenotype matrix 
            for data in reg_data:
                gene = data[0].split('.')[0] 
                pval = data[1] 
                phen_data.append(np.array([gene, pval, reg])) 

        ## sort regional assocs by FDR  
        phen_data = np.array(phen_data) 
        pvals = phen_data[:,1].astype(float) 
        idx = np.argsort(pvals) ## lowest to highest 
        phen_data = phen_data[idx] 

        ## save sorted phenotype permutation assocs 
        with open('{}/{}_{}.txt'.format(out_path, phen, perm), 'w') as f: 
            for p in phen_data: 
                line = '\t'.join(p) + '\n'
                f.write(line) 

        

            
