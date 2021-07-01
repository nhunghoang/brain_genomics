'''
Compute and save the constants that are needed for the 
CCC Simulated Annealing scripts. These are values that 
do not change during the gene selection/search. 

Constants: 
- ordered lists of pairwise common genes 
- Gi: region i expression matrix of shape (subjects, genes)  
- Gj: region j expression matrix of shape (subjects, genes)  
- GiGj: element-wise multiplication of Gi and Gj AFTER subtracting SUBJECT means from both 
- Gi2: element-wise squaring of mean-centered Gi 
- Gi2: element-wise squaring of mean-centered Gj 
- Fs: co-activity(reg1,reg2) array of shape (subjects, 1)

Nhung, June 2021 
'''

import h5py
import numpy as np 
import os 

## output file 
ccc_file = '/data1/rubinov_lab/brain_genomics/analyses_HCP/multi_gene_assoc/ccc_constants.hdf5'
cgenes_file = '/data1/rubinov_lab/brain_genomics/analyses_HCP/multi_gene_assoc/pairwise_common_genes.txt' 
if os.path.exists(cgenes_file): os.remove(cgenes_file) 

## coactivity 
coact_file = '/data1/rubinov_lab/brain_genomics/data_HCP/phenotypes/coactivity.hdf5'
subjs_file = '/data1/rubinov_lab/brain_genomics/data_HCP/timeseries_order.hdf5' 

with h5py.File(coact_file, 'r') as f:
    coacts = {reg_pair: np.array(f[reg_pair]) for reg_pair in f.keys()}
with h5py.File(subjs_file, 'r') as f: 
    subj_order_all = np.array([str(s) for s in np.array(f['subjects'])])

## expression 
expr_dir = '/data1/rubinov_lab/brain_genomics/data_HCP/expression/filtered_quality_r0.3_p0.01'
samps_file = '/data1/rubinov_lab/brain_genomics/data_HCP/predixcan_v8/sample_ids.txt' 

genes = {}; exprs = {}
for expr_file in os.listdir(expr_dir):
    if expr_file[-5:] != '.hdf5': continue
    with h5py.File('{}/{}'.format(expr_dir, expr_file), 'r') as f:
        reg = expr_file.split('.')[0]
        genes[reg] = np.array([g.decode("utf-8") for g in np.array(f['genes'])])
        exprs[reg] = np.array(f['pred_expr']).T ## (samps * genes)
with open(samps_file, 'r') as f:
    samp_order_all = np.array([s.split('\t')[0] for s in f.readlines()])

## subject/sample IDs
id_file = '/data1/rubinov_lab/brain_genomics/data_HCP/neuro_genetic_ids.txt'
samp_to_subj = {}
with open(id_file, 'r') as f:
    f.readline()
    for line in f.readlines():
        [subj, samp] = line.strip().split('\t')
        samp_to_subj[samp] = subj

## only keep if both data exist
subj_idx = []; samp_idx = []
for i in range(samp_order_all.shape[0]): ## samp is smaller list
    try:
        subj_of_samp = samp_to_subj[samp_order_all[i]]
        idx_of_subj = np.argwhere(subj_order_all == subj_of_samp)[0][0]
        samp_idx.append(i)
        subj_idx.append(idx_of_subj)
    except:
        continue

## reorder data
region_pairs = np.array(list(coacts.keys()))
regions = np.array(list(exprs.keys()))
subj_order = subj_order_all[subj_idx]
samp_order = samp_order_all[samp_idx]
for reg_pair in region_pairs:
    coacts[reg_pair] = coacts[reg_pair][subj_idx]
for reg in regions: 
    exprs[reg] = exprs[reg][samp_idx]

## loop through region pairs, compute CCC constants  
pairwise_genes = {} 
ccc_constants = {} 
for reg_pair in region_pairs:
    [reg1, reg2] = reg_pair.split('_')

    ## gene intersect 
    r1_order = np.argsort(genes[reg1])
    r1_genes = genes[reg1][r1_order]
    r1_exprs = exprs[reg1][:,r1_order] 

    r2_order = np.argsort(genes[reg2])
    r2_genes = genes[reg2][r2_order]
    r2_exprs = exprs[reg2][:,r2_order] 
    
    common_genes = np.intersect1d(r1_genes, r2_genes)
    r1_idx = np.nonzero(np.isin(r1_genes, common_genes))[0]
    r2_idx = np.nonzero(np.isin(r2_genes, common_genes))[0]

    Gi = r1_exprs[:, r1_idx] 
    Gj = r2_exprs[:, r2_idx] 

    print('{:>21s} & {:>21s}: {:3d} common genes'.format(reg1, reg2, common_genes.shape[0]))
    with open(cgenes_file, 'a') as cgf: 
        line = '{}_{}:'.format(reg1, reg2) 
        line += (','.join(common_genes))
        cgf.write(line + '\n')

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

    pairwise_genes[reg_pair] = (Gi,Gj)
    ccc_constants[reg_pair] = (GiGj,Gi2,Gj2)

## save CCC constants: coactivity(F), Gi, Gj, GiGj, Gi2, Gj2 
## also save 2D array of corresponding sample & subject IDs
id_array = np.zeros((2, samp_order.size), dtype=int) 
for s,samp in enumerate(samp_order): 
    id_array[0,s] = samp[3:] ## remove 'MH0' to save in hdf5 format
    id_array[1,s] = subj_order[s] 

with h5py.File(ccc_file, 'w') as f: 
    f['IDs'] = id_array 
    for reg_pair in region_pairs: 
        (Gi, Gj) = pairwise_genes[reg_pair]
        (GiGj, Gi2, Gj2) = ccc_constants[reg_pair] 
        Fs = coacts[reg_pair] 
        
        f[reg_pair + '.Gi'] = Gi 
        f[reg_pair + '.Gj'] = Gj 
        f[reg_pair + '.GiGj'] = GiGj
        f[reg_pair + '.Gi2'] = Gi2 
        f[reg_pair + '.Gj2'] = Gj2 
        f[reg_pair + '.Fs'] = Fs
