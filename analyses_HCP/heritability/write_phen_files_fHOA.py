'''

'''

import numpy as np
import h5py 
import numpy as np 
import sys 
import os 

phenotype = sys.argv[1] 
# <KeysViewHDF5 ['alff', 'avgconn', 'connvar', 'falff', 'partcoef12', 'partcoef6', 'tsvar']>

out_dir = '/data1/rubinov_lab/brain_genomics/analyses_HCP/heritability/phen_files_fHOA'
if not os.path.exists(out_dir): os.mkdir(out_dir) 

regions = ['amygdala', 'hippocampus', 'putamenBG', 'aCingulateCortex', 'frontalCortex', 'nAccumbensBG', 'caudateBG', 'cerebellum'] 
regions_v8 = ['amygdala', 'hippocampus', 'putamen', 'anterior-cingulate', 'frontal-pole', 'nucleus-accumbens', 'caudate', 'cerebellar-hemisphere']
rmap = {r7:r8 for r7,r8 in zip(regions, regions_v8)}

subj_file = '/data1/rubinov_lab/brain_genomics_accre/scripts/coexpr_coactivity/sim_annealing/sa_variables.hdf5'  
sfam_file = '/data1/rubinov_lab/brain_genomics/data_HCP/subject_demographics/sampleID_race_familyID_ethnicity.txt' 
phen_file = '/data1/rubinov_lab/brain_genomics_accre/data_HCP/phenotypes.hdf5'

## read phenotype subject order (but the sample ID version)  
with h5py.File(subj_file, 'r') as f: 
    phen_subjs = ['MH0' + str(s) for s in np.array(f['gids'])][1:-1]

## parse subject family file (remove 'T' twin tag)  
fam = {} ## k: subject id, v: family id 
with open(sfam_file, 'r') as f: 
    f.readline()
    info = [line.split('\t') for line in f.readlines()] 
    for subj in info: 
        sid = subj[0]
        fid = subj[2]
        if fid[-1] == 'T': 
            fid = fid[:-1]  
        fam[sid] = fid 

## prepare FIDs in phenotype-subject order
fids = [fam[sid] for sid in phen_subjs]

## read region-specific phenotype values 
reg_phen = {} ## k: region, v: subject array of phen values 
with h5py.File(phen_file, 'r') as f: 
    phen_matrix = np.array(f[phenotype]).astype(str) ## (890, 8)
    for r, region in enumerate(regions):
        reg_phen[region] = phen_matrix[1:-1][:,r] 

## save region-specific .phen files 
## tabulated format: FID, SID, PHEN
if phenotype == 'tsvar': phenotype = 'timeseries_variance' ## BANDAID
if phenotype == 'connvar': phenotype = 'connectivity_variance' ## BANDAID
if phenotype == 'avgconn': phenotype = 'connectivity_mean' ## BANDAID
for region in reg_phen.keys():
    phens = reg_phen[region]
    region_v8 = rmap[region]
    with open('{}/{}_{}.phen'.format(out_dir, phenotype, region_v8), 'w') as f: 
        for i in range(len(phen_subjs)): 
            f.write('{}\t{}\t{}\n'.format(fids[i], phen_subjs[i], phens[i])) 
