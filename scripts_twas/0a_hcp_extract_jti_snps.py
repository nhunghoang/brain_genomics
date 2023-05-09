'''
Save the VCF subset of HCP SNPs that are needed in the JTI models. 

Chromosomal positions from HCP were first mapped to the
NCBI rsids (i.e., the rsids in the HCP VCFs were ignored).

NCBI rsID positions were downloaded on May 8, 2023 from: 
https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/common_all_20180423.vcf.gz 

This file was subsequently filtered to only include rsIDs used in JTI models, 
saved as variable `JTI_path` below. 

- Nhung, May 2023
'''

import os 
import numpy as np 
import pandas as pd
import subprocess 
from multiprocessing import Pool

## paths 
top_path = '/data1/rubinov_lab/brain_genomics/'
JTI_path = top_path + 'scripts_twas/aux_files/JTI_commonSNPs_GRCh37_b151.csv'
vcf_path = top_path + 'data_HCP/Marchini_phg000989/snps_by_chrom_hg37'

out_path = top_path + 'scripts_twas/inputs_HCP/vcf_JTI'

## read the JTI-specific SNP reference 
## k: (chr, pos), v: rsid 
pos_map = pd.read_csv(JTI_path, \
                      sep='\t', \
                      usecols=['ID', 'CHROM', 'POS'], \
                      index_col=['CHROM', 'POS']) \
            .to_dict()['ID']

## func: slice chromosomes
def extract(chrom): 

    path = '{}/chr{}.filtered.sampid.vcf'.format(vcf_path, chrom)
    cpos = subprocess.check_output('cut -f 2 {}'.format(path), shell=True)
    cpos = cpos.decode().split('\n')

    cidx = []; dne = [] 
    for i, pos in enumerate(cpos): 
        try:
            rs = pos_map[(chrom, int(pos))]
            cidx.append(i) 
        except: 
            dne.append(pos) 
         
    print('chr {}: {} total, {} dne, {} in JTI'\
        .format(chrom, len(cidx) + len(dne), len(dne), len(cidx)))

    new_lines = [] 
    with open(path, 'r') as f: 
        for i, line in enumerate(f): 
            if i not in cidx: continue 

            info = line.split('\t')
            [chrom, pos] = info[:2]

            info[2] = pos_map[(int(chrom), int(pos))]
            new_line = '\t'.join(info)
            new_lines.append(new_line)
    
    opath = '{}/c{}.vcf'.format(out_path, chrom)
    with open(opath, 'w') as f: 
        f.writelines(new_lines)

pool = Pool(processes=15)
pool.map(extract, [c for c in np.arange(22,0,-1)])
