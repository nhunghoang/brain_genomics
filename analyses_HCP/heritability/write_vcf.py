'''
For GCTA, save vcf containing only the relevant 
samples. 

-Nhung, July 2021 
'''

import h5py 
import numpy as np 

all_dir = '/data1/rubinov_lab/brain_genomics/data_HCP/Marchini_phg000989/snps_by_chrom_hg38' 
pdx_dir = '/data1/rubinov_lab/brain_genomics/data_HCP/predixcan_v8/vcf_format' 

samp_all_dir = '/data1/rubinov_lab/brain_genomics/analyses_HCP/heritability/vcf_hoacer/all'
samp_pdx_dir = '/data1/rubinov_lab/brain_genomics/analyses_HCP/heritability/vcf_hoacer/pdx'

ts_file = '/data1/rubinov_lab/brain_genomics/data_hoacer/timeseries_order.hdf5'
mp_file = '/data1/rubinov_lab/brain_genomics/analyses_HCP/multi_gene_assoc/ccc_constants.hdf5' 

## gather sample IDs for which we have timeseries for (in timeseries order)  
with h5py.File(mp_file, 'r') as f: IDs = np.array(f['IDs'], dtype=str)
id_map = {IDs[1][s]:'MH0'+IDs[0][s] for s in range(IDs.shape[1])} ## k: subj ID, v: samp ID 

with h5py.File(ts_file, 'r') as f: subjects = np.array(f['subjects'], dtype=str)
samples = [id_map[s] for s in subjects]

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

## write new vcf files containing just the hoacer samples 
for ch in range(21, 0, -1): 
    print(ch) 

    all_vcf = '{}/chr{}.vcf'.format(all_dir, ch)
    with open(all_vcf, 'r') as f: all_lines = f.readlines() 

    pdx_vcf = '{}/chr{}.vcf'.format(pdx_dir, ch)
    with open(pdx_vcf, 'r') as f: pdx_lines = f.readlines() 

    ## all snps 
    f_samp_all = '{}/chr{}.vcf'.format(samp_all_dir, ch) 
    with open(f_samp_all, 'w') as f: 
        f.write(header)
        f.write(new_cat_names) 
        for line in all_lines[5:]:
            line_list = line.strip().split('\t')
            samp_values = np.array(line_list[9:])[samp_idx]  
            othr_values = line_list[:9]
            new_line = '\t'.join(othr_values) + '\t' + '\t'.join(samp_values) + '\n'
            f.write(new_line) 

    ## pdx snps 
    f_samp_pdx = '{}/chr{}.vcf'.format(samp_pdx_dir, ch) 
    with open(f_samp_pdx, 'w') as f: 
        f.write(header)
        f.write(new_cat_names) 
        for line in pdx_lines[5:]:
            line_list = line.strip().split('\t')
            samp_values = np.array(line_list[9:])[samp_idx]  
            othr_values = line_list[:9]
            new_line = '\t'.join(othr_values) + '\t' + '\t'.join(samp_values) + '\n'
            f.write(new_line) 

