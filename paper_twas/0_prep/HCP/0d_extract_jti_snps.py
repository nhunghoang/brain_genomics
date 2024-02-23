'''
Save the VCF subset of HCP SNPs that are needed in the JTI models.  

* Chromosomal positions from HCP were first mapped to the
  NCBI rsids (i.e., the rsids in the HCP VCFs were ignored).

* NCBI rsID positions were downloaded on May 8, 2023 from: 
  https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/common_all_20180423.vcf.gz 

* This file was subsequently filtered to only include rsIDs used in JTI models 
  (which are saved as variable `JTI_path` below). 

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
phn_path = top_path + 'paper_twas/inputs_HCP/phenotypes/all_phens_ids.txt'

## out paths
out_path = top_path + 'paper_twas/inputs_HCP/vcf_JTI'

## read the JTI-specific SNP reference 
## k: (chr, pos), v: rsid 
pos_map = pd.read_csv(JTI_path, \
                      sep='\t', \
                      usecols=['ID', 'CHROM', 'POS'], \
                      index_col=['CHROM', 'POS']) \
            .to_dict()['ID']

## read phenotype order (save SNPs in the same order) 
pids = pd.read_csv(phn_path, sep='\t')['sample'].values

## read vcf order and get corresponding indices 
path = vcf_path + '/chr22.filtered.sampid.vcf'
with open(path, 'r') as f: 
    _ = [f.readline() for i in range(2)] 
    header = f.readline() 

gids = np.array(header.strip().split('\t')[9:])
pidx = [np.where(gids == p)[0][0] for p in pids]

## func: slice chromosomes
def extract(chrom): 

    ## gather all (SNP) positions in this chrom
    path = '{}/chr{}.filtered.sampid.vcf'.format(vcf_path, chrom)
    cpos = subprocess.check_output('cut -f 2 {}'.format(path), shell=True)
    cpos = cpos.decode().split('\n')

    ## check if SNP is in JTI 
    cidx = []; dne = [] 
    for i, pos in enumerate(cpos): 
        try:
            rs = pos_map[(chrom, int(pos))]
            cidx.append(i) 
        except: 
            dne.append(pos) 

    print('chr {}: {} total, {} dne, {} in JTI'\
        .format(chrom, len(cidx) + len(dne), len(dne), len(cidx)))

    ## write JTI SNPs to file 
    new_lines = [] 
    with open(path, 'r') as f: 
        for i, line in enumerate(f): 
            if i not in cidx: continue 

            ## parse chrom, pos, and geno probs
            info = line.strip().split('\t')
            [chrom, pos] = info[:2]
            geno = info[9:][pidx]

            ## update rsid based on chrom and pos 
            info[2] = pos_map[(int(chrom), int(pos))]
            new_info = info[:9]
            new_info.extend(geno)

            new_line = '\t'.join(new_info) + '\n'
            new_lines.append(new_line)

    opath = '{}/c{}.vcf'.format(out_path, chrom)
    with open(opath, 'w') as f: 
        f.writelines(new_lines)
    print('wrote chr {}'.format(chrom))

pool = Pool(processes=22)
pool.map(extract, [c for c in np.arange(22,0,-1)])
