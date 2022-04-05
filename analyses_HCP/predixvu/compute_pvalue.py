'''
Count the number of times the observed PrediXVU results 
appear in the permutations. This script can be ran for 
either permutation type. 

- Nhung, March 2022 
'''

import os 
import sys

ptype = sys.argv[1] 

## paths 
opath = '/data1/rubinov_lab/brain_genomics/analyses_HCP/predixvu/phecodes_genes' 
npath = '/data1/rubinov_lab/brain_genomics/analyses_HCP/predixvu/{}/phecodes_genes'.format(ptype) 

## load observed data 
observed = [] ## [(phenotype, gene, clinical)] 
for phen_file in os.listdir(opath): 
    phen = phen_file.split('.')[0] 
    with open('{}/{}'.format(opath,phen_file), 'r') as f: 
        lines = [line.strip().split('\t') for line in f.readlines()] 
    for line in lines: 
        genes = line[2].split(' ') 
        if len(genes) == 1: 
            observed.append((phen, genes[0], line[1]))
        else: 
            for gene in genes: 
                observed.append((phen, gene, line[1])) 

## load null data 
null = {} ## k: phen, v: [100 searches]  
for phen_file in os.listdir(opath): 
    phen = phen_file.split('.')[0] 
    all_searches = [] ## len 100  
    counts = 0 
    for i in range(100): 
        searches = [] ## num of PrediXVU results  
        nfile = '{}/{}_{}.txt'.format(npath, phen, i)
        with open(nfile, 'r') as f: 
            lines = [line.strip().split('\t') for line in f.readlines()] 
        for line in lines: 
            genes = line[2].split(' ') 
            if len(genes) == 1: 
                searches.append((phen, genes[0], line[1]))
            else: 
                for gene in genes: 
                    searches.append((phen, gene, line[1])) 
        all_searches.append(searches) 
        counts += len(searches)
    null[phen] = all_searches
    print('{}: {} avg results'.format(phen, counts/100)) 
        
## count null instances, compute p-value 
for obsv in observed: 
    count = 0 
    searches = null[obsv[0]]
    for search in searches: 
        if obsv in search: count += 1 
    print('[{:.2f}] {}'.format(count/100, obsv))

