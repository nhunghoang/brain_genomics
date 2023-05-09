'''
Script for parsing full VCF files to generate 
.geno, .snp, and .ind files as defined by 
Eigenstrat documentation. 

doc: https://github.com/DReichLab/EIG/tree/master/CONVERTF

- Nhung, updated May 2023 
'''

import os 
import sys 
import numpy as np 
import pandas as pd 
from time import time 
from multiprocessing import Pool

# ind: sampleID, gender, ancestry 
# snp: snp-name chr Morgan=0 pos ref alt
# geno: 0/1/2/9

## paths 
main_path = '/data1/rubinov_lab/brain_genomics/'
demo_path = main_path + 'data_HCP/subject_demographics/cohort1142_demographics.txt'
geno_path = main_path + 'data_HCP/Marchini_phg000989/snps_by_chrom_hg37'

eign_path = main_path + 'scripts_twas/inputs_HCP/eigen_input/1142' ## .ind, .snp, .geno

## write .ind data
## one line per sample: sampleID gender ancestry 
ind_data = pd.read_table(demo_path, sep='\t', dtype=str, \
                        usecols=['SAMPLE_ID', 'GENDER', 'Race'])
ind_data.to_csv(eign_path + '.ind', sep='\t', \
                index=False, header=False)  

## func: parse lines of chromosome file 
def parse_chrom_subset(data): 

    ch = data['chrom']
    start = data['start']
    end = data['end']
    version = data['version']

    ## line per snp: varID chrom Morgan=0 pos ref alt 
    ## line per gen: 0/1/2/9 for each individual  
    snp_data = [] 
    gen_data = [] 
    cpath = '{}/chr{}.filtered.sampid.vcf'.format(geno_path, ch)

    with open(cpath, 'r') as f: 
        i = 0 
        for line in f: 
            if line[0] == '#': continue     
            
            ## versioning 
            if i < start: continue 
            if i == end: break  

            ## parse 
            info = line.strip().split('\t')  
            name = info[2]
            posn = info[1]
            refa = line[3]
            alta = line[4] 

            desc = '\t'.join([name, ch, '0.0', posn, refa, alta])
            snp_data.append(desc + '\n') 

            ## dosages 
            prob = info[9:] 
            dosg = [np.array(p.split(','), dtype=float) \
                    .argmax() for p in prob] 
            dosg = np.array(dosg, dtype=str) 
            vals = ''.join(dosg)
            gen_data.append(vals + '\n') 
    
            i += 1 

    ## write subset to file 
    snp_file = '{}_c{}_v{}.snp'.format(eign_path, ch, version)
    with open(snp_file, 'w') as f: 
        f.writelines(snp_data)
    
    gen_file = '{}_c{}_v{}.genp'.format(eign_path, ch, version)
    with open(gen_file, 'w') as f: 
        f.writelines(gen_data)

chr_ver = [] 
for i, s in enumerate(np.arange(0, 2286000, 10000)): 
    for c in np.arange(22, 0, -1): 
        chr_ver.append({'chrom':str(c), \
                        'start':s, \
                        'end':s+10000, \
                        'version':str(i) \
                        })

pool = Pool(processes=10) 
pool.map(parse_chrom_subset, chr_ver)

