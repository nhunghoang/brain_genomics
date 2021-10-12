'''

'''

import numpy as np
import h5py 
import numpy as np 
import sys 
import os 

phenotype = sys.argv[1] 
tag = sys.argv[2] 
data = sys.argv[3]

## path to written phen files  
out_dir = '/data1/rubinov_lab/brain_genomics/analyses_HCP/heritability/phen_files_{}'.format(tag)
if not os.path.exists(out_dir): os.mkdir(out_dir) 

## cohort with both expression and phenotypes (n = 890)
samp_subj_file = '/data1/rubinov_lab/brain_genomics/analyses_HCP/subj_samp_assoc_order.hdf5'

## cohort with genotypes (and expression, consequently) (n = 1142) 
samp_demo_file = '/data1/rubinov_lab/brain_genomics_accre/data_HCP/sampleID_race_familyID_ethnicity.txt' 

## cohort with phenotypes (n = 891) 
subj_file = '/data1/rubinov_lab/brain_genomics/data_HCP/{}/timeseries_order.hdf5'.format(data) 
phen_file = '/data1/rubinov_lab/brain_genomics/data_HCP/{}/phenotypes/{}.hdf5'.format(data, phenotype)

## assert that the phenotype array is already in the correct order 
with h5py.File(samp_subj_file, 'r') as f: 
    subj_cohort = np.array(f['subjects'])
    samp_cohort = np.array(['MH0'+s for s in np.array(f['samples'], dtype=str)])
with h5py.File(subj_file, 'r') as f: 
    phen_subjs = np.array(f['subjects']) 

## ignore subjects with incomplete data
X = np.setdiff1d(subj_cohort, phen_subjs)
mask = np.ones_like(subj_cohort, dtype=bool)
for x in X:
    idx = np.where(subj_cohort == x)[0][0]
    mask[idx] = False
subj_cohort = subj_cohort[mask]
samp_cohort = samp_cohort[mask]

if not (np.array_equal(subj_cohort, phen_subjs)):
    print('Error: phenotype subjects are not in the same order as samp/subj cohort.')
    sys.exit()

## FA/MD specific: some subject data missing 
if phenotype in ['fa','md']:
    miss_path = '/data1/rubinov_lab/brain_genomics/data_HCP/hoacer_sn_hth/phenotypes/famd_missing_idx.txt' 
    with open(miss_path, 'r') as f: lines = f.readlines()
    miss_idx = [int(i) for i in lines]

    if subj_cohort.size != 890: 
        print('Error: FA/MD subject array is not of expected size ({} instead of 890)'.format(subj_cohort.size))
        sys.exit()
    mask = np.ones_like(subj_cohort, dtype=bool)
    mask[miss_idx] = False 
    subj_cohort = subj_cohort[mask]
    samp_cohort = samp_cohort[mask]
    
## read region-specific phenotype values 
## order (gen & phen) individuals  
reg_phen = {} ## k: region, v: subject array of phen values 
with h5py.File(phen_file, 'r') as f: 
    for region in f.keys():
        phen_vals = np.array(f[region]).astype(str)
        reg_phen[region] = phen_vals 

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
