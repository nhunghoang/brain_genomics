'''
Functions for reading data files 

- Nhung, July 2020 
'''

import numpy as np 
import os 
import h5py

## HCP expressions (PrediXcan v7) 
## expr_dict = {key: region abbreviation, value: gene list, expression matrix (subjects x genes)}
## return subj_ids, expr_dict
def HCP_v7_expr(region_subset=None): 
    expr_dir = '/data/rubinov_lab/brain_genomics_project/data_MarchiniPrediXcan_v7/v3_expr/'
    if region_subset is None: 
        expr_files = [f for f in os.listdir(expr_dir) if f.split('_')[-1] == 'expression.txt']
    else: 
        expr_files = ['{}_predicted_expression.txt'.format(r) for r in region_subset]
    subj_ids = None
    expr_dict = {} 
    for xfile in expr_files: 
        region = xfile.split('_')[0]
        with open(expr_dir+xfile, 'r') as f: 
            lines = [l.strip().split('\t') for l in f.readlines()]
        if subj_ids is None: subj_ids = np.array([line[0] for line in lines[1:]])
        genes = np.array(lines[0][2:])
        exprs = [line[2:] for line in lines[1:]]
        exprs = np.array(exprs, dtype=float)
        expr_dict[region] = {'genes':genes, 'exprs':exprs}
    return subj_ids, expr_dict

## HCP expressions (PrediXcan v7) 
## (regions x subjects x genes) matrix
## using just genes that all regions have in common
## return {regions, subjects, common genes, expr matrix} 
def HCP_v7_expr_intersect(region_subset=None):
    subjs, expr_dict = HCP_v7_expr(region_subset)
    regions = np.array(list(expr_dict.keys()))
    common_genes = None 
    for reg in regions:
        genes = expr_dict[reg]['genes']
        if common_genes is None: common_genes = genes 
        else: common_genes = np.intersect1d(genes, common_genes) 
    expr_mat = np.zeros((regions.shape[0], subjs.shape[0], common_genes.shape[0]))
    for i, reg in enumerate(regions): 
        gene_idx = np.nonzero(np.isin(expr_dict[reg]['genes'], common_genes))[0]
        expr_mat[i] = expr_dict[reg]['exprs'][:,gene_idx]
    return {'regions':regions, 'subjects':subjs, 'genes':common_genes, 'exprs':expr_mat}

## Allen expressions (for probes with Entrez ids) 
## expr_dict: {key: subject, value: region list, expression matrix (probes x regions)}
## return probes, corresponding genes, expr_dict
def Allen_expr():
    gene_file = '/data/rubinov_lab/brain_genomics_project/scripts/allen/genes_in_probe_order.txt'
    with open(gene_file, 'r') as f: 
        genes = np.array([l.strip() for l in f.readlines()])
    prefix = '/data/rubinov_lab/brain_genomics_project/scripts/allen/allen_probes_regions'
    with open(prefix + '.txt', 'r') as f: 
        lines = [line.strip().split('\t') for line in f.readlines()]
        subjs = np.array(lines[0])
        probes = np.array(lines[1])
        regions = {subjs[i]:np.array(lines[2+i]) for i in range(subjs.shape[0])}
    expr_dict = {}
    with h5py.File(prefix + '.hdf5', 'r') as f: 
        for s in subjs: 
            expr_dict[s] = {'regions':regions[s], 'exprs':np.array(f[s])} 
    return probes, genes, expr_dict
    
## Allen expressions (for probes with Entrez ids)
## (subjects x probes x regions) matrix 
## using just regions that all subjects have in common 
## return {subjects, probes, genes, regions, exprs}
def Allen_expr_intersect():
    probes, genes, expr_dict = Allen_expr()
    subjs = np.array(list(expr_dict.keys()))
    common_regions = None 
    for s in subjs: 
        regs = expr_dict[s]['regions']
        if common_regions is None: common_regions = regs 
        else: common_regions = np.intersect1d(regs, common_regions)
    common_regions = common_regions[:-1] # ignore Right Cerebellum (noise?)
    expr_mat = np.zeros((subjs.shape[0], probes.shape[0], common_regions.shape[0]))
    for i, subj in enumerate(subjs):
        region_idx = np.nonzero(np.isin(expr_dict[subj]['regions'], common_regions))[0]
        expr_mat[i] = expr_dict[subj]['exprs'][:,region_idx]
    return {'subjects':subjs, 'probes':probes, 'genes':genes, 'regions':common_regions, 'exprs':expr_mat}

## region names: map abbrev to HOA, and vice versa 
def map_HOA_region_names(): 
    reg_file = '/data/rubinov_lab/brain_genomics_project/data_MarchiniPrediXcan_v7/updated_brain_regions.txt'
    with open(reg_file, 'r') as f: lines = [line.strip().split('\t') for line in f.readlines()]
    abbrev_to_hoa = {}
    hoa_to_abbrev = {}
    for line in lines[1:]: 
        if line[2] != 'N/A':
            hoa = line[2].split('|')[0] # left region
            abbrev_to_hoa[line[0]] = hoa
            hoa_to_abbrev[hoa] = line[0]
    return abbrev_to_hoa, hoa_to_abbrev

## subject ids: map expression ids to imaging ids, and vice versa
## imaging ids are ints, expression ids are strings 
def map_HCP_subject_ids():
    ifile = '/data/rubinov_lab/brain_genomics_project/data_HCP/neuro_gen_ids.txt'
    gid_to_nid = {}
    nid_to_gid = {}
    with open(ifile, 'r') as f: lines = [i.strip().split('\t') for i in f.readlines()[1:]]
    for line in lines: 
        gid_to_nid[line[1]] = int(line[0])
        nid_to_gid[int(line[0])] = line[1]
    return gid_to_nid, nid_to_gid

## HCP functional connectivity using the HOA 
## for subjects with at least one scan 
## (subjects x regions x regions) FC matrix 
## return subject ids, FC matrix
def HCP_func_conn(regions=None, eqsize=False):
    fc_file = '/data/rubinov_lab/brain_genomics_project/data_HCP/func_conn.hdf5'
    if eqsize:
        fc_file = '/data/rubinov_lab/brain_genomics_project/data_HCP/func_conn_eqsize.hdf5'
    with h5py.File(fc_file, 'r') as f: 
        func_conn = np.array(f['func_conn'])
        subj_ids = np.array(f['subj_ids'])
    if regions is not None: 
        rfile = '/data/rubinov_lab/brain_genomics_project/scripts/HOA_regions.txt'
        with open(rfile, 'r') as f: hoa = [i.strip() for i in f.readlines()]
        idx = [hoa.index(reg) for reg in regions]
        func_conn = func_conn[:,idx,:][:,:,idx]
    return subj_ids, func_conn

## ONLY CALL ONCE (saves files) 
## variables needed specifically for the simulated annealing runs 
def save_simulated_annealing_variables(regions):

    ## gather data
    abbrev2hoa, hoa2abbrev = map_HOA_region_names()
    gid2nid, nid2gid = map_HCP_subject_ids()
    gids0, expr_dict = HCP_v7_expr(regions)
    hoa_regions = [abbrev2hoa[r] for r in regions]
    nids0, func_conn = HCP_func_conn(hoa_regions)

    ## get subjects with both expressions and scans
    nids_expr = [gid2nid[g] for g in gids0]
    nids = np.intersect1d(nids0, nids_expr)
    gids = np.array([nid2gid[n] for n in nids])
    nid_order = [np.where(nids0==n)[0][0] for n in nids]
    gid_order = [np.where(gids0==g)[0][0] for g in gids]

    ## gather exprs and CCC constants for all region pairs
    pairwise_genes = {}
    ccc_consts = {}
    for i, reg1 in enumerate(regions):
        for ii, reg2 in enumerate(regions[i+1:]):
            j = i + 1 + ii

            ## gene intersect
            r1_order = np.argsort(expr_dict[reg1]['genes'])
            r1_genes = expr_dict[reg1]['genes'][r1_order]
            r1_exprs = expr_dict[reg1]['exprs'][:,r1_order]

            r2_order = np.argsort(expr_dict[reg2]['genes'])
            r2_genes = expr_dict[reg2]['genes'][r2_order]
            r2_exprs = expr_dict[reg2]['exprs'][:,r2_order]

            common_genes = np.intersect1d(r1_genes, r2_genes)
            r1_idx = np.nonzero(np.isin(r1_genes, common_genes))[0]
            r2_idx = np.nonzero(np.isin(r2_genes, common_genes))[0]

            Gi = r1_exprs[gid_order][:,r1_idx]
            Gj = r2_exprs[gid_order][:,r2_idx]

            ## expression mean PER SUBJECT > (s,1) vectors
            Mi = np.mean(Gi, axis=1)[:,None]
            Mj = np.mean(Gj, axis=1)[:,None]
            ## centered expressions > (s,g) matrix
            cGi = Gi - Mi
            cGj = Gj - Mj
            ## Gi*Gj, Gi^2, Gj^2 > (s,g) matrix
            GiGj = np.multiply(cGi, cGj)
            Gi2 = np.square(cGi)
            Gj2 = np.square(cGj)

            pairwise_genes[(reg1,reg2)] = (Gi,Gj)
            ccc_consts[(reg1,reg2)] = (GiGj,Gi2,Gj2)    

    ## remove 'MH' prefix from gids, to save in h5py format 
    gids = [int(gid[2:]) for gid in gids]

    ## save Gi, Gj, GiGj, Gi2, Gj2  
    fdir = '/data/rubinov_lab/brain_genomics_project/scripts/coexpr_coactivity/sim_annealing/'
    fname = fdir + 'sa_variables.hdf5'
    f = h5py.File(fname, 'w') 
    f['nid_order'] = np.array(nid_order)
    f['nids'] = np.array(nids)
    f['gid_order'] = np.array(gid_order)
    f['gids'] = np.array(gids)
    for i, reg1 in enumerate(regions):
        for reg2 in regions[i+1:]: 
            (Gi,Gj) = pairwise_genes[(reg1,reg2)]
            (GiGj,Gi2,Gj2) = ccc_consts[(reg1,reg2)]
            mat = np.zeros( (5,Gi.shape[0],Gi.shape[1]), dtype=float ) 
            mat[0] = Gi
            mat[1] = Gj
            mat[2] = GiGj 
            mat[3] = Gi2 
            mat[4] = Gj2
            f['{}_{}'.format(reg1,reg2)] = mat 
    f.close()
    
    

