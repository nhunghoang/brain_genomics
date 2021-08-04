'''

'''

import numpy as np
import h5py 
import numpy as np 
import sys 
import os 

phenotype = sys.argv[1] 

## path to written phen files  
out_dir = '/data1/rubinov_lab/brain_genomics/analyses_HCP/heritability/phen_files_SM'
if not os.path.exists(out_dir): os.mkdir(out_dir) 

## cohort with both expression and phenotypes (n = 891)
samp_subj_file = '/data1/rubinov_lab/brain_genomics/analyses_HCP/subj_samp_assoc_order.hdf5'

## cohort with genotypes (and expression, consequently) (n = 1142) 
samp_demo_file = '/data1/rubinov_lab/brain_genomics_accre/data_HCP/sampleID_race_familyID_ethnicity.txt' 

## cohort with phenotypes (n = 939, out of 1206 possible) 
subj_file = '/data1/rubinov_lab/brain_genomics/data_HCP/timeseries_order.hdf5' 
phen_file = '/data1/rubinov_lab/brain_genomics/data_HCP/phenotypes/{}.hdf5'.format(phenotype)

## fetch order of (gen & phen) individuals (n = 891)
## use to filter and reorder regional phen arrays 
with h5py.File(samp_subj_file, 'r') as f: 
    samp_cohort = np.array(['MH0'+s for s in np.array(f['samples'], dtype=str)])
    subj_cohort_idx = np.array(f['subject_idx_939'], dtype=int)

## read region-specific phenotype values 
## order (gen & phen) individuals  
reg_phen = {} ## k: region, v: subject array of phen values 
with h5py.File(phen_file, 'r') as f: 
    for region in f.keys():
        phen_vals = np.array(f[region]).astype(str)
        reg_phen[region] = phen_vals[subj_cohort_idx] 

## parse subject family file (remove 'T' twin tag)  
fam = {} ## k: subject id, v: family id 
with open(samp_demo_file, 'r') as f: 
    f.readline()
    info = [line.split('\t') for line in f.readlines()] 
    for subj in info: 
        sid = subj[0]
        fid = subj[2]
        if fid[-1] == 'T': 
            fid = fid[:-1]  
        fam[sid] = fid 

## prepare FIDs in phenotype-subject order
fids = [fam[sid] for sid in samp_cohort]

## save region-specific .phen files 
## tabulated format: FID, SID, PHEN
for region in reg_phen.keys():
    phens = reg_phen[region]
    with open('{}/{}_{}.phen'.format(out_dir, phenotype, region), 'w') as f: 
        for i in range(len(samp_cohort)): 
            f.write('{}\t{}\t{}\n'.format(fids[i], samp_cohort[i], phens[i])) 
