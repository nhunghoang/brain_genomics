'''
Write neural phenotypes to hdf5 files. 

Phenotype List: 
- Interregional: CCC 
- Signals: FALFF, ALFF, Regional Homogeneity 
- Connectivity: Variance/Mean, Participation Coefficient  
- Structural: Thickness, Volume  
- Networks: e.g., DMN, Salience, etc.  

- Nhung, May 2021 
'''

import numpy as np 
import h5py 
import os 

PHEN_DIR = '/data1/rubinov_lab/brain_genomics/data_HCP/phenotypes'
MATS_DIR = '/data1/rubinov_lab/brain_genomics/data_HCP/timeseries_seq'
TS_ORDER = '/data1/rubinov_lab/brain_genomics/data_HCP/timeseries_order.hdf5'
COACT_MX = '/data1/rubinov_lab/brain_genomics/data_HCP/phenotypes/coactivity_matrices.hdf5'
PARC_MAT = '/data1/rubinov_lab/brain_genomics/data_HCP/parc-121.mat'

##########################################################################

## (RUN ONCE) Save order of available subject timeseries  
## hdf5 keys: subject-array, subject-to-timeseries-dictionary  
def save_timeseries_order():
    mat_files = os.listdir(MATS_DIR)
    subjects = [] 
    timeseries = {} 
    for mat_file in mat_files: 
        with h5py.File('{}/{}'.format(MATS_DIR, mat_file), 'r') as f: 
            objs = f.get('timeseries') 
            if objs is None: 
                print('{} returned a NoneType timeseries.'.format(mat_file))
                continue 
            ts = [] 
            for t in range(4): 
                obj = objs[0][t] 
                name = h5py.h5r.get_name(obj, f.id) 
                ts_size = np.array(f[name]).size
                if ts_size > 2: 
                    ts.append(t) 
            if len(ts) != 0:
                subj = mat_file.split('.')[0]
                subjects.append(subj)
                timeseries[subj] = ts 
    with h5py.File(TS_ORDER, 'w') as f1:
        f1['subjects'] = np.array(subjects, dtype=int)
        for subj in subjects: 
            f1[subj] = timeseries[subj]
    print('- finished saving timeseries order ({} subj files)'.format(len(subjects))) 

##########################################################################

## Load order of available subject timeseries 
## hdf5 keys: subject-array, subject-to-timeseries-dictionary  
def load_timeseries_order(): 
    with h5py.File(TS_ORDER, 'r') as f: 
        subjects = [str(s) for s in np.array(f['subjects'])] 
        mapping = {s:np.array(f[s]) for s in subjects}
    return subjects, mapping 
    print('- timeseries order loaded') 

##########################################################################

## Compute and write coactivity matrices 
## Note: normalize scan and zero out self-connections 
## hdf5 key: [subject-id] // val: (nscans,121,121) matrix 
def write_coactivity_matrices(subjects, ts_map): 
    cf = h5py.File(COACT_MX, 'w') 
    for subj in subjects: 
        coact_all = [] 
        scan_nums = ts_map[subj]
        with h5py.File('{}/{}.mat'.format(MATS_DIR, subj), 'r') as f: 
            objs = f.get('timeseries') 
            for sn in scan_nums: 
                obj = objs[0][sn] 
                name = h5py.h5r.get_name(obj, f.id) 
                ts = np.array(f[name]).T ## (121,1200) 
                ts_mean = np.mean(ts, axis=1)[:,None] 
                ts_norm = (ts-ts_mean)/ts_mean
                coact = np.corrcoef(ts_norm) 
                np.fill_diagonal(coact, 0) 
                coact_all.append(coact) 
        cf[subj] = np.array(coact_all)
    cf.close()
    print('- finished computing co-activty matrices') 

##########################################################################

## Compute and write connectivity mean and variance for regions of interest 
## hdf5 keys: regions // val: value in subject-scan order 
def write_conn_mean_var(subjs, reg_idx): 
    ## load all subject coact matrices 
    with h5py.File(COACT_MX, 'r') as f: 
        all_coacts = {s:np.array(f[s]) for s in subjs}
    ## set up arrays for subj means and variances 
    nsubjs = len(subjs) 
    all_means = {r:np.zeros(nsubjs, dtype=float) for r in reg_idx} 
    all_varis = {r:np.zeros(nsubjs, dtype=float) for r in reg_idx} 
    ## compute mean/var for each scan across regions THEN average scans together 
    for s,subj in enumerate(subjs): 
        coacts = all_coacts[subj]
        means = np.mean(np.mean(coacts, axis=1), axis=0)
        varis = np.mean(np.var(coacts, axis=1), axis=0)
        for reg,idx in reg_idx.items(): 
            all_means[reg][s] = means[idx]
            all_varis[reg][s] = varis[idx]     
    ## save regional means and variances as subject arrays 
    with h5py.File(PHEN_DIR + '/connectivity_mean.hdf5', 'w') as m: 
        for reg,arr in all_means.items(): 
            m[reg] = arr  
    with h5py.File(PHEN_DIR + '/connectivity_variance.hdf5', 'w') as v: 
        for reg,arr in all_varis.items(): 
            v[reg] = arr  
    print('- finished computing connectivity mean and variance') 

##########################################################################

## Source code copied directly from bctpy (with minor modifications) 
## add citations, couldn't get library to run 
## git: https://github.com/aestrivex/bctpy 
def participation_coefficient(W, ci): 
    _, ci = np.unique(ci, return_inverse=True)
    ci += 1
    n = len(W)  # number of vertices
    Ko = np.sum(W, axis=1)  # (out) degree
    Gc = np.dot((W != 0), np.diag(ci))  # neighbor community affiliation
    Kc2 = np.zeros((n,))  # community-specific neighbors
    for i in range(1, int(np.max(ci)) + 1):
        Kc2 += np.square(np.sum(W * (Gc == i), axis=1))
    P = np.ones((n,)) - Kc2 / np.square(Ko)
    P[np.where(np.logical_not(Ko))] = 0
    return P

def write_part_coef(subjs, ci, reg_idx): 
    ## load all subject coact matrices 
    with h5py.File(COACT_MX, 'r') as f: 
        all_coacts = {s:np.array(f[s]) for s in subjs} 
    subj_pcs = []  
    ## compute PC per subject, for all their coact matrices 
    for subj,cmat in all_coacts.items(): 
        n_mats = cmat.shape[0]
        pcs = np.zeros((n_mats, 121), dtype=float)
        for i in range(n_mats): 
            W = cmat[i] 
            ## zero out negatives
            W = np.where(W<0, W, 0)
            pc = participation_coefficient(W, ci)
            pcs[i] = pc 
        ## represent subj by their avg PC array 
        avg_pc = np.mean(pcs, axis=0)
        subj_pcs.append(avg_pc) 
    ## save subject array of PCs per region 
    subj_pcs = np.array(subj_pcs) ## (n_subjs, n_regs)  
    with h5py.File(PHEN_DIR + '/participation_coefficient.hdf5', 'w') as p:
        for reg,idx in reg_idx.items(): 
            p[reg] = subj_pcs[:,idx]      

            m0 = subj_pcs[:,idx].min()
            m1 = subj_pcs[:,idx].max()
            print('{:>25s}: {:.3f} - {:.3f}'.format(reg, m0, m1))

    print('- finished computing participation coefficient')
         
##########################################################################

## Regional Homogeneity 
def write_reg_homogeneity(): 
    pass 
    print('- finished computing regional homogeneity') 

##########################################################################
##########################################################################

## Main 
def main(): 

    names_file = '/data1/rubinov_lab/brain_genomics/data_HCP/naming-121.txt'
    with open(names_file, 'r') as f: 
        f.readline() ## header  
        region_idx = {} ## k: reg abbrev // v: idx relative to 121 regs  
        for i in f.readlines():
            line = i.strip().split('\t')
            region_idx[line[1]] = int(line[2])
    print('- region naming/indices loaded') 

    #save_timeseries_order() 
    subjects, ts_map = load_timeseries_order() 

    #write_coactivity_matrices(subjects, ts_map) 

    #write_conn_mean_var(subjects, region_idx) 
    #write_reg_homogeneity() 

    clusters_file = '/data1/rubinov_lab/brain_genomics/data_HCP/clusters-subcort-yeo.txt'
    with open(clusters_file, 'r') as f:
        clusters = np.array([int(i.strip()) for i in f.readlines()])
    write_part_coef(subjects, clusters, region_idx)

main() 


