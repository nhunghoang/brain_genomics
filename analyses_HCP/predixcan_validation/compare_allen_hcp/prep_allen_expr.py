'''
Prep Allen expression for co-expression analysis against HCP. 
Per subject, reduce (probes * regional sites) expression matrix 
into (genes * regions) expression matrix for genes with Ensembl 
IDs and regions with PrediXcan models.  

- Nhung, Feb 2022
'''

import os 
import numpy as np 
import h5py 

## paths 
abi_path = '/data1/rubinov_lab/brain_genomics/data_Allen'
out_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/predixcan_validation/compare_allen_hcp/allen_expr.hdf5' 

rand_regs_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/assoc/null_pvals_alff'
parc_regs_file = '/data1/rubinov_lab/brain_genomics/data_HCP/2021_hoacer_sn_hth/naming-115.txt'

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

## convert Allen genes from Entrez IDs to Ensembl IDs, and record probe indices per valid gene  
abi_entrez_all = np.loadtxt(abi_genes_file, dtype=str)
abi_gene_idx = {} ## k: ensembl, v: indices 
for entrez in np.unique(abi_entrez_all): 
    try: ensembl = gene_dict[entrez]
    except: continue 
    abi_gene_idx[ensembl] = np.where(abi_entrez_all == entrez)[0]
abi_genes = list(abi_gene_idx.keys())

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

## gather Allen regional expression per subject, averaging probes per gene 
abi_expr = {} ## k: subj, v: {reg: expr}  
for subj_dir in abi_subjs: 
    expr_path = '{}/{}/MicroarrayExpression.csv'.format(abi_path, subj_dir)
    expr_all = np.loadtxt(expr_path, delimiter=',') ## (all probes, all regions) 
    
    ## average probe expression per gene 
    dims = (len(abi_genes), expr_all.shape[1])
    expr_gene = np.zeros(dims) ## (genes, all regions) 
    for i,gene in enumerate(abi_genes):
        g_idx = abi_gene_idx[gene] 
        expr_gene[i] = np.mean(expr_all[g_idx], axis=0)

    ## then average expression per region 
    reg_idx = abi_reg_idx[subj_dir] ## k: reg-hemi, v: indices  
    subj_expr = {} ## k: reg, v: expr 
    for reg in regions: 
        
        ## case 1: no hemis 
        if reg in ['hypothalamus', 'substantia-nigra']:
            expr_reg = expr_gene[:,reg_idx[reg]] ## (genes, this region)  
            subj_expr[reg] = np.mean(expr_reg, axis=1)

        ## case 2: hemis 
        else: 
            ## left hemi average  
            expr_reg_L = expr_gene[:,reg_idx[reg+'-LH']] ## (genes, this region left)
            avg_L = np.mean(expr_reg_L, axis=1)

            ## case 2A: no right hemi 
            if reg_idx[reg+'-RH'].size == 0: 
                subj_expr[reg] = avg_L  

            ## case 2B: right hemi 
            else: 
                ## right hemi average 
                expr_reg_R = expr_gene[:,reg_idx[reg+'-RH']] ## (genes, this region right)
                avg_R = np.mean(expr_reg_R, axis=1)
                ## overall average 
                subj_expr[reg] = (avg_L + avg_R)/2 ## (genes, this region) 

    abi_expr[subj_dir] = subj_expr 
    print(subj_dir)

## write expression and genes to file 
## hdf5 key: genes, val: ordered list of genes (applicable to all Allen expr data) 
## hdf5 key: [donor#]-[reg], val: expr array 
## hdf5 key: [donor#]-[reg]-sites, val: total number of sites from both hemispheres 
with h5py.File(out_path, 'w') as f: 
    f['genes'] = [g.encode() for g in abi_genes]
    for subj_dir in abi_subjs: 
        subj = subj_dir.split('_')[-1]
        for reg in regions: 
            f['{}-{}'.format(subj,reg)] = abi_expr[subj_dir][reg] 

            if reg in ['hypothalamus', 'substantia-nigra']:
                f['{}-{}-sites'.format(subj,reg)] = reg_idx[reg].size  
            else:
                nsites = reg_idx[reg+'-LH'].size + reg_idx[reg+'-RH'].size
                f['{}-{}-sites'.format(subj,reg)] = nsites  
             

