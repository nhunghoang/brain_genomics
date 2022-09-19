'''
Format neuroimaging phenotypes and save as hdf5 files. 

keywords: region names 
values: ordered subject array 

- Nhung (updated Sept 2022)  
'''

import numpy as np 
import h5py 
import os 
import sys

dset = sys.argv[1] ## HCP/UKB 

## input paths 
SUBJ_ORD = '/data1/rubinov_lab/brain_genomics/analyses_{}/DATA_OUTPUT/subj_samp_assoc_order.hdf5'.format(dset)
PHEN_DIR = '/data1/rubinov_lab/brain_genomics/data_{}/hoacer_sn_hth/phenotypes'.format(dset)
MATS_DIR = '/data1/rubinov_lab/brain_genomics/data_{}/hoacer_sn_hth/timeseries'.format(dset) 
TS_ORDER = '/data1/rubinov_lab/brain_genomics/data_{}/hoacer_sn_hth/timeseries_order.hdf5'.format(dset)
REG_NAME = '/data1/rubinov_lab/brain_genomics/data_{}/hoacer_sn_hth/naming-115.txt'.format(dset)

## phenotypes 
GMVL_MAT = '/data1/rubinov_lab/brain_genomics/neuro_phenotypes/{}/GM_hoacer.mat'.format(dset) 
MYEL_MAT = '/data1/rubinov_lab/brain_genomics/neuro_phenotypes/{}/Myelin_hoacer.mat'.format(dset)
ALFF_MAT = '/data1/rubinov_lab/brain_genomics/neuro_phenotypes/{}/ALFF_hoacer.mat'.format(dset) 
REHO_MAT = '/data1/rubinov_lab/brain_genomics/neuro_phenotypes/{}/Reho_hoacer.mat'.format(dset)   
RHGS_MAT = '/data1/rubinov_lab/brain_genomics/neuro_phenotypes/{}/GSout/Reho_hoacer.mat'.format(dset) 
FCNN_MAT = '/data1/rubinov_lab/brain_genomics/neuro_phenotypes/{}/FCmean_hoacer.mat'.format(dset)   
FCGS_MAT = '/data1/rubinov_lab/brain_genomics/neuro_phenotypes/{}/GSout/FCmean_hoacer.mat'.format(dset)
COACT_MX = '/data1/rubinov_lab/brain_genomics/data_{}/hoacer_sn_hth/coactivity_matrices.hdf5'.format(dset)

##########################################################################

## Load order of available subject timeseries 
## hdf5 keys: subject-array, subject-to-timeseries-dictionary  
def load_timeseries_order(): 
    with h5py.File(TS_ORDER, 'r') as f: 
        subjects = np.array([str(s) for s in np.array(f['subjects'])]) 
        mapping = {s:np.array(f[s]) for s in subjects}
    print('- timeseries order loaded ({} subj files)'.format(subjects.size)) 
    return subjects, mapping 

##########################################################################

## Save Neda's mat files as hdf5 files (just hemisphere means) 
def save_data(in_file, key, out_name, reg_idx):
    with h5py.File(in_file, 'r') as f: 
        data = np.array(f[key]) 
    #assert(data.shape == (890,115))

    out_file = '{}/{}.hdf5'.format(PHEN_DIR, out_name)
    with h5py.File(out_file, 'w') as f: 
        for reg, idxs in reg_idx.items(): 
            f[reg] = data[:,idxs].mean(axis=1) 

    print('done: {}'.format(out_name)) 

##########################################################################

def write_gm_volume(reg_idx): 

    ## load GM volume  
    with h5py.File(GMVL_MAT, 'r') as f: 
        data = np.array(f['GM_115parc'])

    ## BANDAID - remove subject 173233 for having < 1200 timepoints
    bad_index = 446 ## index relative to 891 subject order
    data = np.delete(data, bad_index, axis=0)

    ## write hemisphere averages to file 
    with h5py.File(PHEN_DIR + '/gm_volume.hdf5', 'w') as ff: 
        for reg, idxs in reg_idx.items(): 
            ff[reg] = data[:,idxs].mean(axis=1) 

##########################################################################

def write_conn_mean(subjs, reg_idx): 

    ## load all subject coactivity matrices 
    with h5py.File(COACT_MX, 'r') as f: 
        all_coacts = {s:np.array(f[s]) for s in subjs}

    ## set up region-based arrays 
    all_means = {reg:np.empty(subjs.size) for reg in reg_idx} 

    ## order of averaging: each scan, all scans, hemispheres 
    for s,subj in enumerate(subjs): 
        coacts = all_coacts[subj] 
        cmeans = coacts.mean(axis=1).mean(axis=0) 
        for reg, idxs in reg_idx.items(): 
            all_means[reg][s] = cmeans[idxs].mean()

    ## save regional means as subject arrays 
    with h5py.File(PHEN_DIR + '/connectivity_mean.hdf5', 'w') as m: 
        for reg, means in all_means.items(): 
            m[reg] = means 

##########################################################################
##########################################################################

## Main 
def main(): 

    regions = ['hippocampus', 'amygdala', 'hypothalamus', 'substantia-nigra',\
               'caudate', 'putamen', 'nucleus-accumbens', 'anterior-cingulate',\
               'frontal-pole', 'cerebellar-hemisphere']

    ## k: reg abbrev, v: idxs relative to atlas (115) 
    reg_idx = {reg: [] for reg in regions}  

    with open(REG_NAME, 'r') as f: 
        f.readline()
        lines = f.readlines()  
    for line in lines: 
        atlas, abbrev, index = line.strip().split('\t') 
        if abbrev[-1] == 'H': reg = abbrev[:-3] 
        else: reg = abbrev 
        reg_idx[reg].append(int(index))

    subjects, ts_map = load_timeseries_order() 

    ## call phenotype functions 
    write_gm_volume(reg_idx)
    write_conn_mean(subjects, reg_idx)
    
    files = [MYEL_MAT, ALFF_MAT, REHO_MAT, RHGS_MAT, FCNN_MAT, FCGS_MAT] 
    keys = ['Myelin_115parc', 'ALFF_115parc', 'Reho_allsubjs', 'Reho_allsubjs', \
            'FCmean_parcelvalues', 'FCmean_GSout_parcelvalues'] 
    outs = ['myelination', 'alff', 'regional_homogeneity', 'reho_noGS', 'connectivity_mean', 'connmean_noGS'] 

    for i in range(len(files)): 
        save_data(files[i], keys[i], outs[i], reg_idx)

main() 
