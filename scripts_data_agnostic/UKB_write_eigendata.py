'''
Script for parsing full VCF files to generate 
new .geno, .snp, and .ind files as defined by 
Eigenstrat documentation. 

Note: this script needs to be ran per chromosome. 

doc: https://github.com/DReichLab/EIG/tree/master/CONVERTF

- Nhung (Sept 2021)

Updated for UKB

- Tim (Mar 7 2022)
'''

import os 
import sys 
import numpy as np 
from time import time 

ch = int(sys.argv[1]) ## chromosome 
dset = sys.argv[2] ## dataset, e.g., HCP or UKB 

chr_counts = {10:4381390, 11:4447616, 12:4246710, 13:3130222, 14:2906455, 15:2654777, \
    16:2982874, 17:2546620, 18:2495574, 19:1995346, 1:7081482, 20:2001787, 21:1206469, \
    22:1200863, 2:7807335, 3:6428176, 4:6287665, 5:5826087, 6:5508700, 7:5178847, \
    8:5086660, 9:3909140} ## UKB counts

## input / output paths 
if dset == 'UKB': ##TODO 
    header_dir = '/data1/rubinov_lab/brain_genomics/data_UKB/snps_by_chrom_hg38/lifted_clean_vcf_format'
    id_file = '/data1/rubinov_lab/brain_genomics/data_UKB/subject_demographics/demo_UKB.txt'
    dos_dir = '/data1/rubinov_lab/Yiting/UKB_preprocessing/dosage_format'

eig_dir = '/data1/rubinov_lab/brain_genomics/analyses_{}/DATA_OUTPUT/eigendata'.format(dset)
if not os.path.exists(eig_dir): os.mkdir(eig_dir) 

## load cohort sampleIDs 
if dset == 'UKB': ##TODO
    cohort = np.loadtxt(id_file, delimiter='\t', skiprows=1, dtype=str, usecols=[1])

## load vcf header (which contains all sampleIDs)
header_list = None
vcf = open(header_dir + '/chr22.vcf') 
for i, line in enumerate(vcf): 
    if i == 27: header_list = np.array(line.strip().split('\t'), dtype='<U7')
    if i > 27: break 
vcf.close()


## find the indices of cohort IDs in vcf header (in cohort order) 
_, _, vcf_idx = np.intersect1d(cohort, header_list, assume_unique=True, return_indices=True) 
sorted_idx = np.argsort(cohort) 
cohort_idx = [] ## indices of cohort in vcf (in cohort order)
for i in range(cohort.size): 
    sidx = sorted_idx[i] ## position of sample in sorted order 
    vidx = vcf_idx[sidx] - 3 ## position of sample in dos files ** dos has 3 fewer cols than vcf, hence -3 
    cohort_idx.append(vidx) 
## assert(np.array_equal(cohort, header_list[cohort_idx])) ** this would fail b/o previous lines

## chr, pos, rsid, ref, alt, cohort 
header_idx = [1,2,3,4]
header_idx.extend(cohort_idx)

## ind 
# one line per indiv; sampID gender ancestry 
data = np.loadtxt(id_file, delimiter='\t', skiprows=1, dtype=str, usecols=[1,3,5])
with open('{}/{}_cohort.ind'.format(eig_dir, dset), 'w') as f: 
    for d in data: 
        line = '\t'.join(d) 
        f.write(line + '\n')

## parse chr-specific vcf files 
gen_pos = '0.0' 
chr_file = '{}/chr{}.dosage.txt'.format(dos_dir, ch) 

acc = -1 
p = -1 
num_snps = chr_counts[ch]
num_idvs = cohort.shape[0]
gen_arr = np.empty((num_snps, num_idvs), dtype=str)
snp_arr = np.empty((num_snps, 6), dtype='<U11') 

tstart = time()
with open(chr_file) as f: 
    for line in f: 
        ## parse 
        p += 1 
        d = np.array(line.strip().split('\t'))[header_idx]

        snp = d[1]
        pos = d[0] 
        ref = d[2]
        alt = d[3] 
        dos = [np.array(d.split(','), dtype=float) for d in d[4:]]
        
        ## snp 
        ## one line per snp; snp-name chr Morgan=0 pos ref alt 
        snp_arr[p] = np.array([snp, str(ch), gen_pos, pos, ref, alt])

        ## geno
        ## one line per snp; 0/1/2/9 per indiv
        gen_arr[p] = dos 

        if (p%10000 == 0): 
            print('chr {}: {:.2f}% complete'.format(ch, (p/num_snps)*100))

tpass = time() - tstart 
hr = int(tpass//3600); mn = int(tpass%3600//60); sc = int(tpass%3600%60)
print('chr {}: 100% complete ({} hr, {} mn, {} sc - to parse)'.format(ch, hr, mn, sc))

## wrote gen and snp arrays to file 
gen_file = '{}/chr{}_{}_cohort.geno'.format(eig_dir, ch, dset)
snp_file = '{}/chr{}_{}_cohort.snp'.format(eig_dir, ch, dset)

tstart = time() 
np.savetxt(gen_file, gen_arr[:p+1], delimiter='', fmt='%s')
np.savetxt(snp_file, snp_arr[:p+1], delimiter='\t', fmt='%s')
tpass = time() - tstart 

hr = int(tpass//3600); mn = int(tpass%3600//60); sc = int(tpass%3600%60)
print('chr {}: {} hr, {} mn, {} sc (to write)'.format(ch, hr, mn, sc))



