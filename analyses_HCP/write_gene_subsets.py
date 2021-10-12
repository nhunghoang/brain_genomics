'''

'''

import h5py 
import numpy as np 
import os 

expr_dir = '/data1/rubinov_lab/brain_genomics/data_HCP/expression/filtered_quality_r0.3_p0.01' 
idx_file = '/data1/rubinov_lab/brain_genomics/analyses_HCP/subj_samp_assoc_order.hdf5'
pval_dir = '/data1/rubinov_lab/brain_genomics/analyses_HCP/single_gene_assoc/genes_by_phen' 

phens = ['regional_homogeneity', 'alff', 'gm_volume']

## sample/subject order 
with h5py.File(idx_file, 'r') as f:
    samp_idx = np.array(f['sample_idx_1142'])

# filter genes 
filter_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/genes_filter'
if not os.path.exists(filter_path): os.mkdir(filter_path)

single_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/genes_single'
if not os.path.exists(single_path): os.mkdir(single_path)

for expr_file in os.listdir(expr_dir):
    if expr_file[-4:] != 'hdf5': continue
    region = expr_file.split('.')[0]
    if region in ['cerebellum', 'cortex', 'spinal-cord']: continue 

    ## read filter subset 
    with h5py.File('{}/{}'.format(expr_dir, expr_file), 'r') as f:
        expr = np.array(f['pred_expr'])[:,samp_idx] ## (genes * subjects)
        genes = np.array(f['genes']) ## (genes) 
        #genes = np.array([g.decode("utf-8") for g in np.array(f['genes'])]) ## (genes) 

    ## write filter subset 
    with h5py.File('{}/{}.hdf5'.format(filter_path, region), 'w') as f:
        f['expr'] = expr
        f['genes'] = genes

    ## read single subset 
    for phen in phens:  
        in_path = '{}/{}/pearsonr_{}.txt'.format(pval_dir, phen, region) 
        out_path = '{}/{}_{}.hdf5'.format(single_path, phen, region) 

        gdata = np.loadtxt(in_path, dtype=str, delimiter='\t', skiprows=1) 
        assert(np.array_equal(np.array([g.decode("utf-8") for g in genes]), gdata[:,0]))
        
        pvals = gdata[:,3].astype(float) 
        mask = np.zeros_like(pvals, dtype=bool) 
        mask[pvals <= 0.05] = True 
        pv_expr = expr[mask] 
        pv_genes = genes[mask]

        print('{:>20s} - {:<22s}: {}'.format(phen, region, mask.sum()))

        ## write single subset 
        with h5py.File(out_path, 'w') as f: 
            f['expr'] = pv_expr 
            f['genes'] = pv_genes 

