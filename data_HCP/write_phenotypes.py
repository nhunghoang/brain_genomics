'''
Write neural phenotypes to hdf5 files. 

Phenotype List: 
- Interregional: CCC 
- Signals: FALFF, ALFF, Regional Homogeneity 
- Connectivity: Variance/Mean, Participation Coefficient  
- Structural: Thickness, Volume  
- Networks: e.g., DMN, Salience, etc.  

Nhung, May 2021 (updated June 2021)  
'''

import numpy as np 
import h5py 
import os 
from scipy.stats import pearsonr, spearmanr 

PHEN_DIR = '/data1/rubinov_lab/brain_genomics/data_HCP/phenotypes'
MATS_DIR = '/data1/rubinov_lab/brain_genomics/data_HCP/timeseries'
TS_ORDER = '/data1/rubinov_lab/brain_genomics/data_HCP/timeseries_order.hdf5'
COACT_MX = '/data1/rubinov_lab/brain_genomics/data_HCP/coactivity_matrices.hdf5'
PARC_MAT = '/data1/rubinov_lab/brain_genomics/data_HCP/parc-121.mat'
REHO_MAT = '/data1/rubinov_lab/brain_genomics/data_HCP/phenotypes/regional_homogeneity.mat' 

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

## Write coactivity arrays for region pairs of interest 
## Order of operations: LH/RH specific coactivity, average between hemis, average across scans  
def write_coactivity(subjs, reg_idx): 

    ## define region pairs alphabetically
    ## ignore hemisphere tag
    regions = []
    for reg in reg_idx.keys(): 
        if reg[-1] == 'H': regions.append(reg[:-3])
        else: regions.append(reg)
    regions = np.unique(regions) 
    region_pairs = [(r1,r2) for i,r1 in enumerate(regions) for r2 in regions[i+1:]]

    ## gather LH & RH indices for every region pair 
    region_pair_idx = {} ## k: (reg1,reg2) // v: [(LH idxs), (RH idxs)] 
    for (areg,breg) in region_pairs: 
        aLH = 0; bLH = 0
        aRH = 0; bRH = 0 
        
        if areg in ['hypothalamus', 'substantia_nigra']: 
            aLH = reg_idx[areg] 
            aRH = reg_idx[areg]
        else: 
            aLH = reg_idx[areg + '-LH']
            aRH = reg_idx[areg + '-RH']

        if breg in ['hypothalamus', 'substantia-nigra']: 
            bLH = reg_idx[breg]
            bRH = reg_idx[breg]
        else: 
            bLH = reg_idx[breg + '-LH']
            bRH = reg_idx[breg + '-RH']

        region_pair_idx[(areg,breg)] = [(aLH,bLH), (aRH, bRH)]
        
    ## load all subject coactivity matrices 
    with h5py.File(COACT_MX, 'r') as f: 
        all_coacts = {s:np.array(f[s]) for s in subjs}
    nsubjs = len(subjs) 

    ## compute (for each region pair) an array of coactivity across subjects   
    coactivity = {rp:np.zeros(nsubjs, dtype=float) for rp in region_pairs} 
    for s,subj in enumerate(subjs): 
        coacts = all_coacts[subj] ## (scans, 121, 121) 
        for reg_pair in region_pairs: 
            [LH, RH] = region_pair_idx[reg_pair]

            ## A) get corr(A,B) in LH and RH separately (and per scan)
            coacts_LH = coacts[:, LH[0], LH[1]]
            coacts_RH = coacts[:, RH[0], RH[1]]

            ## B) average LH-corr and RH-corr together (still per scan)  
            coacts_hm = (coacts_LH + coacts_RH) / 2

            ## C) average corr across scans and save  
            coactivity[reg_pair][s] = np.mean(coacts_hm) 
        
    ## save coactivity per region pair as subject arrays  
    with h5py.File(PHEN_DIR + '/coactivity.hdf5', 'w') as c: 
        for reg_pair, subj_arr in coactivity.items():
            key = '_'.join(reg_pair)
            c[key] = subj_arr 
    print('- finished computing coactivity') 

##########################################################################

## Compute and write connectivity mean and variance for regions of interest 
## Check phenotype correlation between LH & RH, save average 
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

    ## check correlation between LH & RH values, store average     
    h_means = {}; h_varis = {} 
    c_means = {}; c_varis = {} 
    for reg in reg_idx.keys(): 
        if (reg=='hypothalamus') or (reg=='substantia-nigra'): 
            h_means[reg] = all_means[reg]   
            h_varis[reg] = all_varis[reg] 
        else: 
            reg_full = reg[:-3]
            m1 = all_means[reg] 
            v1 = all_varis[reg]
            try:
                m0 = h_means[reg_full] 
                v0 = h_varis[reg_full]
                h_means[reg_full] = (m0 + m1) / 2 
                h_varis[reg_full] = (v0 + v1) / 2 
                mr, mp = spearmanr(m0, m1)
                vr, vp = spearmanr(v0, v1)
                c_means[reg_full] = mr
                c_varis[reg_full] = vr
            except KeyError:
                h_means[reg_full] = m1 
                h_varis[reg_full] = v1 

    ## print correlation info 
    print('{:>25s}  {:7s} {:7s}'.format('REGION', 'MEANCORR', 'VARCORR'))
    for reg in c_means.keys(): 
        print('{:>25s}  {:.3f}    {:.3f}'.format(reg, c_means[reg], c_varis[reg]))

    ## save regional means and variances as subject arrays 
    with h5py.File(PHEN_DIR + '/connectivity_mean.hdf5', 'w') as m: 
        for reg,arr in h_means.items(): 
            m[reg] = arr  
    with h5py.File(PHEN_DIR + '/connectivity_variance.hdf5', 'w') as v: 
        for reg,arr in h_varis.items(): 
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
            '''
            Ko = np.sum(W, axis=1)
            if np.any(Ko == 0): 
                line = '{}-{}: {}\n'.format(subj, i, np.where(Ko==0)[0])
                print(line)            
            '''
            pc = participation_coefficient(W, ci)
            pcs[i] = pc 

        ## represent subj by their avg PC array 
        avg_pc = np.mean(pcs, axis=0)
        subj_pcs.append(avg_pc) 

    subj_pcs = np.array(subj_pcs) ## (n_subjs, n_regs)  

    ## check correlation between LH & RH, store average 
    print('{:>25s}  {:7s} {:7s} {:7s}'.format('REGION', 'CORRPC', 'MINPC', 'MAXPC'))
    h_pcs = {} 
    for reg,idx in reg_idx.items(): 
        reg_pcs = subj_pcs[:,idx] 
        if (reg=='hypothalamus') or (reg=='substantia-nigra'): 
            h_pcs[reg] = reg_pcs 
        else: 
            reg_full = reg[:-3] 
            if reg_full not in h_pcs.keys(): 
                h_pcs[reg_full] = reg_pcs 
            else: 
                h0 = h_pcs[reg_full] 
                vals = (h0 + reg_pcs)/2
                h_pcs[reg_full] = vals 
                rho, pval = spearmanr(h0, reg_pcs) 
                print('{:>25s}   {:.3f}  {:.3f}   {:.3f}'.\
                    format(reg_full, rho, vals.min(), vals.max()))

    ## save subject array of PCs per region 
    with h5py.File(PHEN_DIR + '/participation_coefficient.hdf5', 'w') as p:
        for reg,vals in h_pcs.items(): 
            p[reg] = vals      

    print('- finished computing participation coefficient')
         
##########################################################################

## Read ReHo values (computed in matlab, per regional hemisphere)
## Check phenotype correlation between LH & RH, save average 
## script: /data1/rubinov_lab/brain_genomics/data_HCP/compute_ReHo.m
def write_regional_homogeneity(): 
    print('{:>25s}  {:7s} {:7s} {:7s}'.format('REGION', 'CORRRH', 'MINRH', 'MAXRH'))
    rehos = {} ## k: region, v: subject array  
    with h5py.File(REHO_MAT, 'r') as f: 
        for reg in f.keys(): 
            val = np.array(f[reg])[0]
            if (reg=='hypothalamus') or (reg=='substantia_nigra'): 
                rehos[reg] = val 
            else: 
                reg_full = reg[:-3] 
                try: 
                    r0 = rehos[reg_full] 
                    avg = (r0 + val)/2
                    rehos[reg_full] = avg 
                    rho, pval = spearmanr(r0, val) 
                    print('{:>25s}   {:.3f}  {:.3f}   {:.3f}'.\
                        format(reg_full, rho, avg.min(), avg.max()))
                except KeyError: 
                    rehos[reg_full] = val 

    ## save subject array of ReHos per region 
    with h5py.File(PHEN_DIR + '/regional_homogeneity.hdf5', 'w') as rh:
        for reg0,arr in rehos.items(): 
            reg = reg0.replace('_', '-')
            rh[reg] = arr

    print('- finished saving regional homogeneity')

##########################################################################

## Read ALFF and FALFF values (computed by Neda, script path below)
## Check phenotype correlation between LH & LR, save average 
## script: /data1/rubinov_lab/Neda/ALFF_fALFF_updated.m
def write_alff_falff(subjs, reg_idx): 
    alff_path = '/data1/rubinov_lab/brain_genomics/neuro_phenotypes/HCP_phenotypes/ALFF_121regions.mat'
    falff_path = '/data1/rubinov_lab/brain_genomics/neuro_phenotypes/HCP_phenotypes/fALFF_121regions.mat'
    subj_path = '/data1/rubinov_lab/brain_genomics/neuro_phenotypes/HCP_phenotypes/subjs_list2.mat' 

    ## read subject order of phenotype files 
    with h5py.File(subj_path, 'r') as f: 
        phen_subjs = np.array(f['list'])[0].astype(int).astype('str')

    ## get indices of relevant subjects in phen array 
    subj_idx = [np.where(phen_subjs==s)[0][0] for s in subjs]

    ## read data - shape (1206, 121)  
    with h5py.File(alff_path, 'r') as f: all_alff = np.array(f['ALFF_subjcs'])  
    with h5py.File(falff_path, 'r') as f: all_falff = np.array(f['fALFF_subjcs']) 

    ## slice the parts of interest, compute hemisphere averages   
    print('{:>25s}  {:7s} {:7s} {:7s} {:7s} {:7s} {:7s}'.\
        format('REGION', 'CORR-a', 'MIN-a', 'MAX-a', 'CORR-f', 'MIN-f', 'MAX-f'))
    alff = {}; falff = {} 
    for reg0, idx in reg_idx.items(): 
        aval = all_alff[subj_idx][:,idx] 
        fval = all_falff[subj_idx][:,idx]  
        if (reg0=='hypothalamus') or (reg0=='substantia-nigra'): 
            alff[reg0] = aval; falff[reg0] = fval 
        else: 
            reg = reg0[:-3]
            try: 
                h_alff = alff[reg] 
                h_falff = falff[reg] 
                alff[reg] = (h_alff + aval) / 2
                falff[reg] = (h_falff + fval) / 2
                arho, apval = spearmanr(h_alff, aval) 
                frho, fpval = spearmanr(h_falff, fval) 
                print('{:>25s}   {:.3f}  {:.3f}   {:.3f}   {:.3f}   {:.3f}   {:.3f}'.\
                    format(reg, arho, alff[reg].min(), alff[reg].max(),\
                    frho, falff[reg].min(), falff[reg].max()))
            except KeyError: 
                alff[reg] = aval; falff[reg] = fval 

    ## save subject array of ALFF/FALFF per region 
    with h5py.File(PHEN_DIR + '/alff.hdf5', 'w') as aa:
        for reg,val in alff.items(): 
            aa[reg] = val 
    with h5py.File(PHEN_DIR + '/falff.hdf5', 'w') as ff:
        for reg,val in falff.items(): 
            ff[reg] = val 

    print('- finished saving ALFF/FALFF')

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

    #write_coactivity(subjects, region_idx)

    #write_conn_mean_var(subjects, region_idx) 

    #clusters_file = '/data1/rubinov_lab/brain_genomics/data_HCP/clusters-subcort-yeo.txt'
    #with open(clusters_file, 'r') as f:
    #    clusters = np.array([int(i.strip()) for i in f.readlines()])
    #write_part_coef(subjects, clusters, region_idx)

    #write_regional_homogeneity() 

    write_alff_falff(subjects, region_idx) 

main() 


