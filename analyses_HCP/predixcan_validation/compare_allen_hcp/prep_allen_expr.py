'''
Prep Allen expression for co-expression analysis against HCP. 
Per subject, reduce (probes * regional sites) expression matrix 
into (probes-as-genes * regions) expression matrix for genes with 
Ensembl IDs and regions with PrediXcan models. Normalize the expression 
of each probe across all sampling sites per subject prior to reducing.  
Also, only consider the left-hemisphere. 

- Nhung, updated March 2022
'''

import os 
import numpy as np 
import h5py 

## paths 
abi_path = '/data1/rubinov_lab/brain_genomics/data_Allen'
out_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/predixcan_validation/compare_allen_hcp/allen_expr.hdf5' 

rand_regs_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/assoc/null_pvals_alff'
parc_regs_file = '/data1/rubinov_lab/brain_genomics/data_HCP/hoacer_sn_hth/naming-115.txt'

gene_ids_file = '/data1/rubinov_lab/brain_genomics/analyses_HCP/predixcan_validation/compare_allen_hcp/2022_ensembl_entrez.txt' 
abi_genes_file = '/data1/rubinov_lab/brain_genomics/data_Allen/probes_as_entrez.txt'

## gather regions, and map to parc names 
regions = sorted(os.listdir(rand_regs_path)) ## (no hemi) 
with open(parc_regs_file, 'r') as f: 
    f.readline() 
    lines = [i.split('\t') for i in f.readlines()] 
reg_dict = {j[1]:j[0] for j in lines} ## k: reg abbrev (incl. hemi), v: parc reg name  

## map biomaRt Entrez gene IDs to Ensembl IDs, skip if no mapping  
with open(gene_ids_file, 'r') as f: 
    f.readline() 
    lines = [i.strip().split(',') for i in f.readlines()]
gene_dict = {j[0]:j[1] for j in lines if j[0] != ''} ## k: Entrez, v: Ensembl 

## convert Allen genes from Entrez IDs to Ensembl IDs, and record probe indices of valid genes  
abi_entrez_all = np.loadtxt(abi_genes_file, dtype=str)
abi_probe_genes = [] 
abi_probe_idx = [] 
for p,entrez in enumerate(abi_entrez_all): 
    try: ensembl = gene_dict[entrez]
    except: continue 
    abi_probe_genes.append(ensembl) 
    abi_probe_idx.append(p) 

## gather indices of Allen regions (incl. hemi) per subject 
abi_subjs = [i for i in os.listdir(abi_path) if i.split('_')[0] == 'normalized']
abi_reg_idx = {} ## k: subj, v: {reg-hemi:indices}
for subj_dir in abi_subjs: 
    parcel_path = '{}/{}/hoacer_sn_hth_regions.txt'.format(abi_path, subj_dir)
    abi_parcels = np.loadtxt(parcel_path, delimiter='\t', usecols=[1], dtype=str)  
    subj_reg_idx = {} ## k: reg-hemi, v: indices 
    for reg_hemi, parc_hemi in reg_dict.items():
        subj_reg_idx[reg_hemi] = np.where(abi_parcels==parc_hemi)[0]
    abi_reg_idx[subj_dir] = subj_reg_idx  

## gather Allen regional expression per subject  
dims = (len(abi_subjs), len(abi_probe_idx))
abi_expr = {r:np.zeros(dims) for r in regions} ## k: reg, v: (subjs * probes-as-genes) 
for s,subj_dir in enumerate(abi_subjs): 
    expr_path = '{}/{}/MicroarrayExpression.csv'.format(abi_path, subj_dir)
    expr_all = np.loadtxt(expr_path, delimiter=',') ## (all probes * all regions) 
    expr_all = expr_all[:,1:] ## first column is probe ID 
    
    ## normalize expression matrix 
    #expr_mean = expr_all.mean(axis=1)[:,None] ## probe means 
    #expr_stdv = expr_all.std(axis=1)[:,None] 
    #expr_all = (expr_all - expr_mean) / expr_stdv
    #expr_all = expr_all / expr_mean
    
    ## slice valid probes  
    expr_gene = expr_all[abi_probe_idx] ## (valid probes * all regions) 

    ## then average expression per region -- only consider left hemisphere  
    reg_idx = abi_reg_idx[subj_dir] ## k: reg-hemi, v: indices  
    for reg in regions: 
        
        ## case 1: no hemis 
        if reg in ['hypothalamus', 'substantia-nigra']:
            expr_reg = expr_gene[:,reg_idx[reg]] ## (valid probes * this region)  
            abi_expr[reg][s] = np.mean(expr_reg, axis=1)

        ## case 2: hemis 
        else: 
            ## left hemi average  
            expr_reg_L = expr_gene[:,reg_idx[reg+'-LH']] ## (valid probes * this region left)
            avg_L = np.mean(expr_reg_L, axis=1)
            abi_expr[reg][s] = avg_L

    print(subj_dir)

## write expression and genes to file 
## hdf5 key: subjects, val: ordered list of subjects in expr data 
## hdf5 key: probes-as-genes, val: ordered list of probes by their gene names (applicable to all Allen expr data) 
## hdf5 key: [reg], val: expr matrix (subjects * probes) 
## hdf5 key: [donor#]-[reg]-sites, val: number of sites in left hemisphere 
with h5py.File(out_path, 'w') as f: 
    f['subjects'] = [s.split('_')[-1].encode() for s in abi_subjs] 
    f['genes'] = [g.encode() for g in abi_probe_genes]
    
    for reg in regions: 
        f[reg] = abi_expr[reg]

    for subj_dir in abi_subjs: 
        subj = subj_dir.split('_')[-1]
        for reg in regions: 
            if reg in ['hypothalamus', 'substantia-nigra']:
                nsites = abi_reg_idx[subj_dir][reg].size  
            else:
                nsites = abi_reg_idx[subj_dir][reg+'-LH'].size
            f['{}-{}-sites'.format(subj,reg)] = nsites  
