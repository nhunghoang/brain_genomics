'''
Write neural phenotypes to hdf5 files. 

Phenotype List: 
- Interregional: CCC 
- Signals: FALFF, ALFF, Regional Homogeneity, Variance  
- Connectivity: Variance/Mean, Participation Coefficient  
- Structural: Thickness, Volume  
- Networks: e.g., DMN, Salience, etc.  

Nhung, May 2021 (updated July 2021)  

'''

import numpy as np 
import h5py 
import os 
from scipy.stats import pearsonr, spearmanr 
import sys

data_dir = sys.argv[1]

## input data 
SUBJ_ORD = '/data1/rubinov_lab/brain_genomics/analyses_HCP/subj_samp_assoc_order.hdf5'
PHEN_DIR = '/data1/rubinov_lab/brain_genomics/data_HCP/{}/phenotypes'.format(data_dir)
MATS_DIR = '/data1/rubinov_lab/brain_genomics/data_HCP/{}/timeseries'.format(data_dir) 
TS_ORDER = '/data1/rubinov_lab/brain_genomics/data_HCP/{}/timeseries_order.hdf5'.format(data_dir)
REG_NAME = '/data1/rubinov_lab/brain_genomics/data_HCP/{}/naming-115.txt'.format(data_dir)
#PARC_MAT = '/data1/rubinov_lab/brain_genomics/data_HCP/parc-121.mat'

## phenotypes precomputed by Nhung
COACT_MX = '/data1/rubinov_lab/brain_genomics/data_HCP/{}/coactivity_matrices.hdf5'.format(data_dir)
REHO_DIR = '/data1/rubinov_lab/brain_genomics/data_HCP/{}/phenotypes/reho_by_subj'.format(data_dir)  

## phenotypes precomputed by Neda 
FAMD_MAT = '/data1/rubinov_lab/Neda/Diffusion/Diffusion_FA_MD_hoacer_sn_hth.mat'
GMVL_MAT = '/data1/rubinov_lab/Neda/GM/GM_vol_hoacer_sn_hth.mat'
MYEL_MAT = '/data1/rubinov_lab/Neda/Myelin/Myelination_hoacer_sn_hth.mat' 
ALFF_MAT = '/data1/rubinov_lab/Neda/ALFF/ALFF_hoacer_sn_hth_vox_parcelled.mat'
FALF_MAT = '/data1/rubinov_lab/Neda/ALFF/fALFF_hoacer_sn_hth_vox_parcelled.mat'

##########################################################################

## (RUN ONCE) Save order of available subject timeseries  
## hdf5 keys: subject-array, subject-to-timeseries-dictionary  
def save_timeseries_order():
    with h5py.File(SUBJ_ORD, 'r') as f: subj_list = np.array(f['subjects'], dtype=str)

    subjects = [] 
    timeseries = {} 
    for subj in subj_list: 
        if data_dir in ['hoacer_sn_hth', 'hoacer_sn_hth_nofilt']: mat_file = subj.strip() + '.mat'
        else: mat_file = subj.strip() + '/processed_hoacer.mat'
        mat_path = MATS_DIR + '/' + mat_file
        if not os.path.exists(mat_path): continue 

        with h5py.File(mat_path, 'r') as f: 
            objs = f.get('Vp_clean') 
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
                if data_dir in ['hoacer_sn_hth', 'hoacer_sn_hth_nofilt']:
                    subj = mat_file.split('.')[0]
                else: 
                    subj = mat_file.split('/')[0]
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
    print('- timeseries order loaded ({} subj files)'.format(len(subjects))) 
    return subjects, mapping 

##########################################################################

## Compute and write coactivity matrices 
## Note: normalize scan and zero out self-connections 
## hdf5 key: [subject-id] // val: (nscans,121,121) matrix 
def write_coactivity_matrices(subjects, ts_map): 
    cf = h5py.File(COACT_MX, 'w') 
    for subj in subjects: 
        coact_all = [] 
        scan_nums = ts_map[subj]
        if data_dir in ['hoacer_sn_hth', 'hoacer_sn_hth_nofilt']: mat_file = subj.strip() + '.mat'
        else: mat_file = subj.strip() + '/processed_hoacer.mat'
        mat_path = MATS_DIR + '/' + mat_file
        with h5py.File(mat_path, 'r') as f: 
            objs = f.get('Vp_clean') 
            for sn in scan_nums: 
                obj = objs[0][sn] 
                name = h5py.h5r.get_name(obj, f.id) 
                ts = np.array(f[name]).T ## (115,1140) 
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
        pcs = np.zeros((n_mats, 115), dtype=float)
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
## script: /data1/rubinov_lab/brain_genomics/data_HCP/compute_reho.m
def write_regional_homogeneity(subjs): 

    ## get region order that values are saved in 
    with open(REG_NAME, 'r') as f:
        f.readline() ## header 
        regs = [line.split('\t')[1] for line in f.readlines()] 

    ## get subject values and organize by region 
    subj_rehos = np.zeros((len(subjs), len(regs))) ## (num subjs, num regs) 
    for s,subj in enumerate(subjs):
        with h5py.File('{}/{}.mat'.format(REHO_DIR, subj), 'r') as f: 
            subj_reho = np.array(f['reho'])
        subj_rehos[s] = np.squeeze(subj_reho) 
            
    ## store average values across hemispheres per region 
    print('{:>25s}  {:7s} {:7s} {:7s}'.format('REGION', 'CORRRH', 'MINRH', 'MAXRH'))
    rehos = {} ## k: region, v: subject array  
    for r,reg in enumerate(regs): 
        vals = subj_rehos[:,r]

        if reg in ['hypothalamus', 'substantia-nigra']:
            rehos[reg] = vals  

        else:
            name = reg[:-3]
            try: 
                vals0 = rehos[name] 
                avg = (vals0 + vals) / 2 
                rehos[name] = avg 
                rho, pval = spearmanr(vals0, vals) 
                print('{:>25s}   {:.3f}  {:.3f}   {:.3f}'.\
                    format(name, rho, avg.min(), avg.max()))
            except KeyError: 
                rehos[name] = vals  

    ## save subject array of ReHos per region 
    with h5py.File(PHEN_DIR + '/regional_homogeneity.hdf5', 'w') as rh:
        for reg,arr in rehos.items(): 
            rh[reg] = arr

    print('- finished saving regional homogeneity')

##########################################################################

## Read ALFF and FALFF values (computed by Mika, script path below)
## Check phenotype correlation between LH & LR, save average 
## script: /data1/rubinov_lab/brain_genomics/data_HCP/hoacer_filt/preprocess_hcp_hoacer.m 
def write_alff_falff(subjs, ts_map, reg_idx): 

    ## set up region-specific subject arrays  
    nsubjs = len(subjs)
    all_alff = {r:np.zeros(nsubjs, dtype=float) for r in reg_idx}
    all_falff = {r:np.zeros(nsubjs, dtype=float) for r in reg_idx}

    ## read computed values from mat, take subject averages  
    for s,subj in enumerate(subjs):
        with h5py.File('{}/{}/processed_hoacer.mat'.format(MATS_DIR, subj), 'r') as f:
            subj_alff = []; a_objs = f.get('Vp_alff')
            subj_falff = []; f_objs = f.get('Vp_falff')
            for scan_num in ts_map[subj]: 
                a_obj = a_objs[0][scan_num] 
                f_obj = f_objs[0][scan_num] 
                a_name = h5py.h5r.get_name(a_obj, f.id)
                f_name = h5py.h5r.get_name(f_obj, f.id)
            
                alff = np.array(f[a_name]).T ## (115,1140) 
                falff = np.array(f[f_name]).T ## (115,1140) 
                if alff.size < 3: continue  
                subj_alff.append(alff)
                subj_falff.append(falff)
        avg_subj_alff = np.mean(subj_alff, axis=0)
        avg_subj_falff = np.mean(subj_falff, axis=0)
        for reg,idx in reg_idx.items():
            all_alff[reg][s] = avg_subj_alff[idx]
            all_falff[reg][s] = avg_subj_falff[idx]

    ## check correlation between LH & RH values
    print('{:>25s}   {:7s} {:7s}   {:7s} {:7s} {:7s} {:7s}'.\
        format('REGION', 'CORR-a', 'MIN-a', 'MAX-a', 'CORR-f', 'MIN-f', 'MAX-f'))
    alff = {}; falff = {} 
    for reg0, idx in reg_idx.items(): 
        aval = all_alff[reg0]; fval = all_falff[reg0]  
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
                print('{:>25s}   {:.3f}   {:.3f}   {:.3f}   {:.3f}   {:.3f}   {:.3f}'.\
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

## Compute and write signal variance across (normalized) timepoints per region
## hdf5 keys: regional hemispheres // val: value in subject-scan order
def write_ts_variance(subjs, ts_map, reg_idx):

    ## set up array for subject timeseries variances
    nsubjs = len(subjs)
    all_tsvars = {r:np.zeros(nsubjs, dtype=float) for r in reg_idx}

    ## compute tsvar for each scan across regions THEN average scans together
    for s,subj in enumerate(subjs):
        scan_nums = ts_map[subj]
        scan_tsvars = np.zeros((scan_nums.size, 115), dtype=float)

        if data_dir in ['hoacer_sn_hth', 'hoacer_sn_hth_nofilt']: mat_file = subj.strip() + '.mat'
        else: mat_file = subj.strip() + '/processed_hoacer.mat'
        mat_path = MATS_DIR + '/' + mat_file

        with h5py.File(mat_path, 'r') as f:
            objs = f.get('Vp_clean')
            for n,sn in enumerate(scan_nums):
                obj = objs[0][sn]
                name = h5py.h5r.get_name(obj, f.id)
                ts = np.array(f[name]).T ## (121,1200)
                ts_mean = np.mean(ts, axis=1)[:,None]
                ts_norm = (ts-ts_mean)/ts_mean
                tsvars = np.var(ts_norm, axis=1)
                scan_tsvars[n] = tsvars
        avg_scan_tsvars = np.mean(scan_tsvars, axis=0)
        for reg,idx in reg_idx.items():
            all_tsvars[reg][s] = avg_scan_tsvars[idx]

    ## check correlation between LH & RH values
    h_tsvars = {}; c_tsvars = {}
    for reg in reg_idx.keys():
        if (reg=='hypothalamus') or (reg=='substantia-nigra'):
            h_tsvars[reg] = all_tsvars[reg]
        else:
            reg_full = reg[:-3]
            t1 = all_tsvars[reg]
            try:
                t0 = h_tsvars[reg_full]
            except KeyError:
                h_tsvars[reg_full] = t1
                continue
            h_tsvars[reg_full] = (t0 + t1) / 2
            tr, tp = spearmanr(t0, t1)
            c_tsvars[reg_full] = tr

    ## print correlation info
    print('{:>25s}  {}'.format('REGION', 'TSVAR_CORR'))
    for reg in c_tsvars.keys():
        print('{:>25s}  {:.3f}'.format(reg, c_tsvars[reg]))

    ## save regional means and variances as subject arrays
    with h5py.File(PHEN_DIR + '/timeseries_variance.hdf5', 'w') as v:
        for reg,arr in h_tsvars.items():
            v[reg] = arr 
    with h5py.File(PHEN_DIR + '/h_timeseries_variance.hdf5', 'w') as v:
        for reg,arr in all_tsvars.items():
            v[reg] = arr  
    print('- finished computing timeseries variance')

##########################################################################

## Read ALFF and FALFF values (computed by Neda)
## Check phenotype correlation between LH & LR, save average 
def NEDA_write_alff_falff(reg_idx): 

    ## load ALFF and FALFF values  
    with h5py.File(ALFF_MAT, 'r') as f: a_data = np.array(f['ALFF_115parc'])
    with h5py.File(FALF_MAT, 'r') as f: f_data = np.array(f['fALFF_115parc'])

    ## BANDAID - remove subject 173233 for having < 1200 timepoints 
    bad_index = 446 ## index relative to 891 subject order 
    a_data = np.delete(a_data, bad_index, axis=0) 
    f_data = np.delete(f_data, bad_index, axis=0)  

    ## gather values from relevant hemispheric regions  
    a_hem = {}; f_hem = {} ## k: reg hem, v: subj array  
    a_avg = {}; f_avg = {} ## k: reg, v: subj array  
    print('{:>25s}  {:7s}  {:7s}'.format('REGION', 'CORR-AA', 'CORR-FF'))
    for reg, idx in reg_idx.items(): 
        a_val = a_data[:,idx] 
        f_val = f_data[:,idx] 

        a_hem[reg] = a_val
        f_hem[reg] = f_val
    
        ## compute region averages between hemispheres 
        if reg in ['hypothalamus', 'substantia-nigra']:
            a_avg[reg] = a_val 
            f_avg[reg] = f_val 
        else: 
            try: 
                a_avg[reg[:-3]] = (a_avg[reg[:-3]] + a_val) / 2
                f_avg[reg[:-3]] = (f_avg[reg[:-3]] + f_val) / 2
                ## check correlation between LH & RH values
                a_corr, _ = spearmanr(a_avg[reg[:-3]], a_val)
                f_corr, _ = spearmanr(f_avg[reg[:-3]], f_val)
                print('{:>25s}  {:>.5f}  {:>.5f}'.format(reg[:-3], a_corr, f_corr))
            except KeyError: 
                a_avg[reg[:-3]] = a_val 
                f_avg[reg[:-3]] = f_val

    ## save subject array of ALFF/FALFF per region 
    with h5py.File(PHEN_DIR + '/alff.hdf5', 'w') as aa: 
        for reg, val in a_avg.items(): 
            aa[reg] = val 
    with h5py.File(PHEN_DIR + '/falff.hdf5', 'w') as ff: 
        for reg, val in f_avg.items(): 
            ff[reg] = val 

    print('- finished saving ALFF/FALFF')

##########################################################################

## Read FA and MD values (computed by Neda)
## Check phenotype correlation between LH & LR, save average 
def NEDA_write_diffusion_fa_md(reg_idx): 

    ## load FA and MD values  
    with h5py.File(FAMD_MAT, 'r') as f: 
        fa_data = np.array(f['FA_115parc'])
        md_data = np.array(f['MD_115parc'])

    ## BANDAID - remove subject 173233 for having < 1200 timepoints 
    bad_index = 446 ## index relative to 891 subject order 
    fa_data = np.delete(fa_data, bad_index, axis=0)  
    md_data = np.delete(md_data, bad_index, axis=0) 

    ## record indices of subjects with missing data 
    nans = np.sum(np.isnan(fa_data), axis=1)
    miss = np.argwhere(nans > 0)
    good = np.squeeze(np.argwhere(nans == 0))
    fa_data = np.squeeze(fa_data[good])
    md_data = np.squeeze(md_data[good])

    ## gather values from relevant hemispheric regions, removing missing data  
    fa_hem = {}; md_hem = {} ## k: reg hem, v: subj array  
    fa_avg = {}; md_avg = {} ## k: reg, v: subj array  
    print('{:>25s}  {:7s}  {:7s}'.format('REGION', 'CORR-FA', 'CORR-MD'))
    for reg, idx in reg_idx.items(): 
        fa_val = fa_data[:,idx] 
        md_val = md_data[:,idx] 

        fa_hem[reg] = fa_val
        md_hem[reg] = md_val
    
        ## compute region averages between hemispheres 
        if reg in ['hypothalamus', 'substantia-nigra']:
            fa_avg[reg] = fa_val 
            md_avg[reg] = md_val 
        else: 
            try: 
                fa_avg[reg[:-3]] = (fa_avg[reg[:-3]] + fa_val) / 2
                md_avg[reg[:-3]] = (md_avg[reg[:-3]] + md_val) / 2
                ## check correlation between LH & RH values
                fa_corr, _ = spearmanr(fa_avg[reg[:-3]], fa_val)
                md_corr, _ = spearmanr(md_avg[reg[:-3]], md_val)
                print('{:>25s}  {:>.5f}  {:>.5f}'.format(reg[:-3], fa_corr, md_corr))
            except KeyError: 
                fa_avg[reg[:-3]] = fa_val 
                md_avg[reg[:-3]] = md_val

    ## save subject array of FA/MD per region 
    with h5py.File(PHEN_DIR + '/fa.hdf5', 'w') as ff: 
        for reg, val in fa_avg.items(): 
            ff[reg] = val 
    with h5py.File(PHEN_DIR + '/md.hdf5', 'w') as mm: 
        for reg, val in md_avg.items(): 
            mm[reg] = val 

    ## save missing subject data 
    with open(PHEN_DIR + '/famd_missing_idx.txt', 'w') as ss:
        for idx in miss: 
            ss.write('{}\n'.format(idx[0]))

    print('- finished saving FA/MD')

##########################################################################

## Read gray matter volume values (computed by Neda)
## Check phenotype correlation between LH & LR, save average 
def NEDA_write_gm_volume(reg_idx): 

    ## load GM volume  
    with h5py.File(GMVL_MAT, 'r') as f: 
        data = np.array(f['GM_115parc'])

    ## BANDAID - remove subject 173233 for having < 1200 timepoints 
    bad_index = 446 ## index relative to 891 subject order 
    data = np.delete(data, bad_index, axis=0)  

    ## gather values from relevant hemispheric regions  
    gm_hem = {} ## k: reg hem, v: subj array  
    gm_avg = {} ## k: reg, v: subj array  
    print('{:>25s}  {}'.format('REGION', 'CORR-GM-VOL'))
    for reg, idx in reg_idx.items(): 
        gm_val = data[:,idx] 
        gm_hem[reg] = gm_val
    
        ## compute region averages between hemispheres 
        if reg in ['hypothalamus', 'substantia-nigra']:
            gm_avg[reg] = gm_val 
        else: 
            try: 
                gm_avg[reg[:-3]] = (gm_avg[reg[:-3]] + gm_val) / 2
                ## check correlation between LH & RH values
                gm_corr, _ = spearmanr(gm_avg[reg[:-3]], gm_val)
                print('{:>25s}  {:>.5f}'.format(reg[:-3], gm_corr))
            except KeyError: 
                gm_avg[reg[:-3]] = gm_val 

    ## save subject array of GM volume per region 
    with h5py.File(PHEN_DIR + '/gm_volume.hdf5', 'w') as ff: 
        for reg, val in gm_avg.items(): 
            ff[reg] = val 

    print('- finished saving GM volume')

##########################################################################

## Read myelination values (computed by Neda)
## Check phenotype correlation between LH & LR, save average 
def NEDA_write_myelination(reg_idx): 

    ## load myelination volume  
    with h5py.File(MYEL_MAT, 'r') as f: 
        data = np.array(f['Myelin_115parc'])

    ## BANDAID - remove subject 173233 for having < 1200 timepoints 
    bad_index = 446 ## index relative to 891 subject order 
    data = np.delete(data, bad_index, axis=0)  

    ## gather values from relevant hemispheric regions  
    my_hem = {} ## k: reg hem, v: subj array  
    my_avg = {} ## k: reg, v: subj array  
    print('{:>25s}  {}'.format('REGION', 'CORR-MYELINATION'))
    for reg, idx in reg_idx.items(): 
        my_val = data[:,idx] 
        my_hem[reg] = my_val
    
        ## compute region averages between hemispheres 
        if reg in ['hypothalamus', 'substantia-nigra']:
            my_avg[reg] = my_val 
        else: 
            try: 
                my_avg[reg[:-3]] = (my_avg[reg[:-3]] + my_val) / 2
                ## check correlation between LH & RH values
                my_corr, _ = spearmanr(my_avg[reg[:-3]], my_val)
                print('{:>25s}  {:>.5f}'.format(reg[:-3], my_corr))
            except KeyError: 
                my_avg[reg[:-3]] = my_val 

    ## save subject array of GM volume per region 
    with h5py.File(PHEN_DIR + '/myelination.hdf5', 'w') as ff: 
        for reg, val in my_avg.items(): 
            ff[reg] = val 

    print('- finished saving GM volume')

##########################################################################
##########################################################################

## Main 
def main(): 

    with open(REG_NAME, 'r') as f: 
        f.readline() ## header  
        region_idx = {} ## k: reg abbrev // v: idx relative to 115 regs  
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

    #write_regional_homogeneity(subjects) 

    #write_alff_falff(subjects, ts_map, region_idx) 

    #write_ts_variance(subjects, ts_map, region_idx)     

    #NEDA_write_diffusion_fa_md(region_idx) 
    #NEDA_write_gm_volume(region_idx) 
    #NEDA_write_myelination(region_idx)
    NEDA_write_alff_falff(region_idx)

main() 


