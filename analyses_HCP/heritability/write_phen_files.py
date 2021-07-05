'''

'''

import numpy as np
import h5py 
import numpy as np 
import sys 
import os 

phenotype = sys.argv[1] 

## cohort with both expression and phenotypes (n = 891)
samp_subj_file = '/data1/rubinov_lab/brain_genomics/analyses_HCP/multi_gene_assoc/ccc_constants.hdf5'  

## cohort with genotypes (and expression, consequently) (n = 1142) 
samp_demo_file = '/data1/rubinov_lab/brain_genomics_accre/data_HCP/sampleID_race_familyID_ethnicity.txt' 

## cohort with phenotypes (n = 939, out of 1206 possible) 
subj_file = '/data1/rubinov_lab/brain_genomics/data_HCP/timeseries_order.hdf5' 
phen_file = '/data1/rubinov_lab/brain_genomics/data_HCP/phenotypes/{}.hdf5'.format(phenotype)

## path to written phen files  
out_dir = '/data1/rubinov_lab/brain_genomics/analyses_HCP/heritability/phen_files'
if not os.path.exists(out_dir): os.mkdir(out_dir) 

## fetch order of (gen & phen) individuals (n = 891)
## use to filter and reorder regional phen arrays 
with h5py.File(samp_subj_file, 'r') as f: 
    all_people = np.array(f['IDs'])
    phen_subjs = all_people[1]  
    phen_samps = ['MH0' + str(s) for s in all_people[0]]
with h5py.File(subj_file, 'r') as f: 
    all_subjs = np.array(f['subjects'])
subj_idx = [np.where(all_subjs==s)[0][0] for s in phen_subjs]

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
fids = [fam[sid] for sid in phen_samps]

## read region-specific phenotype values and order (gen & phen) individuals  
reg_phen = {} ## k: region, v: subject array of phen values 
with h5py.File(phen_file, 'r') as f: 
    for region in f.keys():
        phen_vals = np.array(f[region]).astype(str)
        reg_phen[region] = phen_vals[subj_idx] 

## save region-specific .phen files 
## tabulated format: FID, SID, PHEN
for region in reg_phen.keys():
    phens = reg_phen[region]
    with open('{}/{}_{}.phen'.format(out_dir, phenotype, region), 'w') as f: 
        for i in range(len(phen_samps)): 
            f.write('{}\t{}\t{}\n'.format(fids[i], phen_samps[i], phens[i])) 
