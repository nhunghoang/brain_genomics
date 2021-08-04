'''
For GCTA, save vcf containing only the relevant 
samples. 

-Nhung, July 2021 
'''

import h5py 
import numpy as np 
from time import time 
import sys 
from multiprocessing import Pool 

A = int(sys.argv[1])
B = int(sys.argv[2])

## pretty-print time
def pretty_time(tpass):
    hr = int(tpass//3600)
    mn = int(tpass%3600//60)
    sc = int(tpass%3600%60)
    return '{:d} hr, {:d} mn, {:d} sc'.format(hr, mn, sc)

## input paths 
all_dir = '/data1/rubinov_lab/brain_genomics/data_HCP/Marchini_phg000989/snps_by_chrom_hg38' 
pdx_dir = '/data1/rubinov_lab/brain_genomics/data_HCP/predixcan_v8/vcf_format' 

## output paths 
samp_all_dir = '/data1/rubinov_lab/brain_genomics/analyses_HCP/heritability/vcf/all'
samp_pdx_dir = '/data1/rubinov_lab/brain_genomics/analyses_HCP/heritability/vcf/pdx'

## gather corresponding IDs for analysis cohort (n = 891)  
map_file = '/data1/rubinov_lab/brain_genomics/analyses_HCP/multi_gene_assoc/ccc_constants.hdf5' 
with h5py.File(map_file, 'r') as f: IDs = np.array(f['IDs'], dtype=str)
samples = ['MH0'+s for s in IDs[0]]
print('{} samples'.format(len(samples)))

## prepare header (and note the sample columns to ignore) 
with open(all_dir + '/chr22.vcf', 'r') as f:
    lines = [f.readline() for line in range(5)]
header = ''.join(lines[:4])
cat_names = lines[4].strip().split('\t')
non_samp_cats = '\t'.join(cat_names[:9]) + '\t'
all_samples = cat_names[9:] 

samp_subset, all_idx, samp_idx = np.intersect1d(all_samples, samples, \
    return_indices=True, assume_unique=True)

new_cat_names = non_samp_cats + '\t'.join(samp_subset) + '\n' 

## write new vcf files containing just the samples being used 
#for ch in range(22, 0, -1): 
def one_run(ch):
    tstart = time() 

    log = open('chr{}.log'.format(ch), 'w') 

    ## all snps 
    all_vcf = '{}/chr{}.vcf'.format(all_dir, ch)
    with open(all_vcf, 'r') as f: all_lines = f.readlines() 
    log.write('[{:>2}] read ALL variants\n'.format(ch))

    f_samp_all = '{}/chr{}.vcf'.format(samp_all_dir, ch) 
    with open(f_samp_all, 'w') as f: 
        f.write(header)
        f.write(new_cat_names) 
        clen = len(all_lines[5:])
        for c, line in enumerate(all_lines[5:]):
            line_list = line.strip().split('\t')
            samp_values = np.array(line_list[9:])[all_idx]  
            othr_values = line_list[:9]
            new_line = '\t'.join(othr_values) + '\t' + '\t'.join(samp_values) + '\n'
            f.write(new_line) 
            if (c%100000 == 0): log.write('[{:>2}] {:.2f}% (all)\n'.format(ch, (c/clen)*100)); log.flush()
    log.write('{:>4} wrote ALL variants\n'.format(''))

    ## pdx snps 
    pdx_vcf = '{}/chr{}.vcf'.format(pdx_dir, ch)
    with open(pdx_vcf, 'r') as f: pdx_lines = f.readlines() 
    log.write('{:>4} read PDX variants\n'.format(''))

    f_samp_pdx = '{}/chr{}.vcf'.format(samp_pdx_dir, ch) 
    with open(f_samp_pdx, 'w') as f: 
        f.write(header)
        f.write(new_cat_names) 
        clen = len(all_lines[5:])
        for c, line in enumerate(all_lines[5:]):
            line_list = line.strip().split('\t')
            samp_values = np.array(line_list[9:])[all_idx]  
            othr_values = line_list[:9]
            new_line = '\t'.join(othr_values) + '\t' + '\t'.join(samp_values) + '\n'
            f.write(new_line) 
            if (c%100000 == 0): log.write('[{:>2}] {:.2f}% (pdx)\n'.format(ch, (c/clen)*100)); log.flush() 
    log.write('{:>4} wrote PDX variants\n'.format(''))

    tpass = pretty_time(time()-tstart)
    log.write('{:>4} finished in {}\n'.format('', tpass))
    log.close()

## implement parallel computing 
chrs = [c for c in range(A, B+1)]
pool = Pool()
pool.map(one_run, chrs) 
