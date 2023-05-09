'''
Script for parsing full VCF/dosage files to generate 
new .geno, .snp, and .ind files as defined by 
Eigenstrat documentation. This script needs to be 
ran per chromosome. Only the cohort of interest is 
included (i.e., not all subjects).   

For HCP, VCF files are parsed. 
For UKB, dosage files are parsed. 

doc: https://github.com/DReichLab/EIG/tree/master/CONVERTF

- Nhung (Sept 2021), updated by Tim (March 2022) 
'''

import os 
import sys 
import numpy as np 
from time import time 

ch = int(sys.argv[1]) ## chromosome 
dset = sys.argv[2] ## HCP or UKB 

## number of SNPs per chromosome in the VCF/dosage file  
hcp_counts = {10:1342893, 11:1329939, 12:1288862, 13:959601, 14:877972, 15:796319, \
    16:862214, 17:753613, 18:757090, 19:615502, 1:2110772, 20:593277, 21:372780, \
    22:138911, 2:2285131, 3:1910134, 4:1940830, 5:1750448, 6:1715824, 7:1574295, \
    8:1503483, 9:1169483}   
ukb_counts = {10:4381390, 11:4447616, 12:4246710, 13:3130222, 14:2906455, 15:2654777, \
    16:2982874, 17:2546620, 18:2495574, 19:1995346, 1:7081482, 20:2001787, 21:1206469, \
    22:1200863, 2:7807335, 3:6428176, 4:6287665, 5:5826087, 6:5508700, 7:5178847, \
    8:5086660, 9:3909140}
if dset == 'HCP': chr_counts = hcp_counts 
if dset == 'UKB': chr_counts = ukb_counts 

## paths
if dset == 'HCP':
    id_file = '/data1/rubinov_lab/brain_genomics/data_HCP/subject_demographics/cohort890_demographics.txt'
    vcf_dir = '/data1/rubinov_lab/brain_genomics/data_HCP/Marchini_phg000989/snps_by_chrom_hg38'
    gen_dir = vcf_dir 

if dset == 'UKB': 
    id_file = '/data1/rubinov_lab/brain_genomics/data_UKB/subject_demographics/demo_UKB.txt'
    vcf_dir = '/data1/rubinov_lab/brain_genomics/data_UKB/snps_by_chrom_hg38/lifted_clean_vcf_format' 
    gen_dir = '/data1/rubinov_lab/Yiting/UKB_preprocessing/dosage_format' 

eig_dir = '/data1/rubinov_lab/brain_genomics/scripts_assoc_clean/inputs_{}/eigen_input'.format(dset)
if not os.path.exists(eig_dir): os.mkdir(eig_dir) 

## load cohort sampleIDs 
cohort = np.loadtxt(id_file, delimiter='\t', skiprows=1, dtype=str, usecols=[1])

## load vcf header (which contains all sampleIDs)
header_list = None
if dset == 'HCP': x = 4 
if dset == 'UKB': x = 27 
vcf = open(vcf_dir + '/chr22.vcf') 
for i, line in enumerate(vcf): 
    if i == x: 
        header_list = np.array(line.strip().split('\t'))
    if i > x: 
        break 
vcf.close()

## find the indices of cohort IDs in vcf header (in cohort order) 
_, _, vcf_idx = np.intersect1d(cohort, header_list, assume_unique=True, return_indices=True) 
sorted_idx = np.argsort(cohort) 
cohort_idx = [] ## indices of cohort in vcf (in cohort order)
for i in range(cohort.size): 
    sidx = sorted_idx[i] ## position of sample in sorted order 
    if dset == 'HCP':
        vidx = vcf_idx[sidx] ## position of sample in vcf header  
    if dset == 'UKB':
        vidx = vcf_idx[sidx]-3 ## dos has 3 fewer cols than vcf, hence -3
    cohort_idx.append(vidx) 
#assert(np.array_equal(cohort, header_list[cohort_idx]))

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
if dset == 'HCP': chr_file = '{}/chr{}.vcf'.format(gen_dir, ch) 
if dset == 'UKB': chr_file = '{}/chr{}.dosage.txt'.format(gen_dir, ch)

acc = -1 
p = -1 
num_snps = chr_counts[ch]
num_idvs = cohort.shape[0]
gen_arr = np.empty((num_snps, num_idvs), dtype=str)
snp_arr = np.empty((num_snps, 6), dtype='<U11') 

tstart = time()
with open(chr_file) as f: 
    for line in f: 

        ## skip header info (HCP only)  
        if dset == 'HCP': 
            acc += 1 
            if acc < 5: continue 

        ## otherwise, parse 
        p += 1 
        d = np.array(line.strip().split('\t'))[header_idx] 

        if dset == 'HCP': snp = d[1].split(',')[0]
        if dset == 'UKB': snp = d[1]
        pos = d[0] 
        ref = d[2]
        alt = d[3] 
        if dset == 'HCP':
            gtp = [np.array(d.split(','), dtype=float) for d in d[4:]]
            dos = np.array([np.argmax(p) for p in gtp])
        if dset == 'UKB': 
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


