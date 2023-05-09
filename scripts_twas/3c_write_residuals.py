'''
Regress confounders (PCs, age, gender) from 
each gene, as well as each phenotype. 

- Nhung (updated April 2023) 
'''

import numpy as np 
import pandas as pd 
import h5py
from sklearn.linear_model import LinearRegression

## paths 
top_path = '/data1/rubinov_lab/brain_genomics/'

hcp_paths = {}

ukb_paths = { \
    'phen': top_path + 'data_UKB/downloads/phenotypes_40k.csv', \
    'expr': top_path + 'scripts_twas/inputs_UKB/grex_JTI', \
    'covs': top_path + 'scripts_twas/inputs_UKB/covariates.csv', \
    'oute': top_path + 'scripts_twas/inputs_UKB/regress_expr_JTI', \
    'outp': top_path + 'scripts_twas/inputs_UKB/regress_phen', \
} 

## func: regress out all covariates 
## X: (n_samples, n_features)
## y: (n_samples) 
def compute_y_residuals(cov_data, y): 
    lr = LinearRegression().fit(cov_data, y) 
    coefs = lr.coef_ ## (n_features) 
    b = lr.intercept_ ## float 
    ccovs = np.matmul(cov_data, coefs) ## (n_samples) 
    y_residuals = y - ccovs - b ## (n_samples) 
    return y_residuals  

## func: format HCP phenotypes 

## func: get UKB phenotypes and covariates 
def get_UKB_data(paths): 

    regions = ['hippocampus', 'amygdala', 'caudate', 'putamen', \
               'cerebellar-hemisphere', 'nucleus-accumbens']

    leftt_ids = ['26562-2.0', '26563-2.0', '26559-2.0', '26560-2.0', '26557-2.0', '26564-2.0']
    right_ids = ['26593-2.0', '26594-2.0', '26590-2.0', '26591-2.0', '26588-2.0', '26595-2.0']

    keys = {} 
    for i in range(len(regions)): 
        keys[leftt_ids[i]] = 'left_vol_{}'.format(REGS[i])
        keys[right_ids[i]] = 'right_vol_{}'.format(REGS[i])
    
    covs_data = pd.read_csv(paths['covs'])
    eids = covs_data['eid'] 
    covs_vals = covs_data.values[:,1:]

    phen_data = pd.read_csv(paths['phen'], index_col='eid')
    phen_data = phen_data.reindex(eids)
    phen_data.rename(columns=keys, inplace=True)

    return phen_data, covs_vals 

## func: apply covariate regression to phenotypes  
def apply_phen_regr(phen_data, cov_mat, out_path): 

    for i in ['left_vol', 'right_vol']: 
        data = {} ## k: reg, v: subj vector  
        for reg in REGS: 
            col = '{}_{}'.format(i, reg) 
            yres = compute_y_residuals(cov_mat, phen_data[col].values)
            data[reg] = yres
            
        out_file = '{}/{}.hdf5'.format(out_path, i)
        with h5py.File(out_file, 'w') as f: 
            for reg, vals in data.items():
                f[reg] = vals

        print(i)

## func: apply covariate regression to expression 
def apply_expr_regr(expr_path, cov_mat, out_path): 
    
    for reg in REGS: 
        expr_file = '{}/{}.hdf5'.format(expr_path, reg) 
        with h5py.File(expr_file, 'r') as f: 
            genes = f['genes'][()] 
            expr = f['pred_expr'][()] 

        expr_regress = np.zeros_like(expr)
        for g in range(genes.size):
            yres = compute_y_residuals(cov_mat, expr[g])
            expr_regress[g] = yres 

        out_file = '{}/{}.hdf5'.format(out_path, reg)
        with h5py.File(out_file, 'w') as f: 
            (n_genes, n_samples) = expr_regress.shape
            n_gene_chunks = np.min([n_genes, 10])

            ds_preds = f.create_dataset('pred_expr', data=expr_regress, \
                                        chunks=(n_gene_chunks, n_samples), \
                                        dtype=np.dtype('float32'), \
                                        scaleoffset=4, \
                                        compression='gzip')

            ds_genes = f.create_dataset('genes', data=genes, \
                                        dtype=h5py.string_dtype())

        print(reg, 'expr')


phen_data, cov_mat = get_UKB_data(ukb_paths)
#apply_phen_regr(phen_data, cov_mat, ukb_paths['outp'])
apply_expr_regr(ukb_paths['expr'], cov_mat, ukb_paths['oute'])

import sys; sys.exit()









dset = sys.argv[1] ## HCP, UKB

## paths 
cov_file = '/data1/rubinov_lab/brain_genomics/scripts_assoc_clean/inputs_{}/covariates.txt'.format(dset)
phen_dir = '/data1/rubinov_lab/brain_genomics/scripts_assoc_clean/inputs_{}/phenotypes'.format(dset)

expr_dir = '/data1/rubinov_lab/brain_genomics/scripts_assoc_clean/inputs_{}/predixcan_grex_v8'.format(dset)
if dset == 'UKB': expr_dir += '/cohort941'
#if dset == 'HCP': expr_dir += '/cohort890'

phen_out = '/data1/rubinov_lab/brain_genomics/scripts_assoc_clean/inputs_{}/phen_regress'.format(dset)
expr_out = '/data1/rubinov_lab/brain_genomics/scripts_assoc_clean/inputs_{}/expr_regress'.format(dset)
#if not os.path.exists(phen_out): os.mkdir(phen_out)
#if not os.path.exists(expr_out): os.mkdir(expr_out)

## region names 
phen_regions = ['amygdala', 'anterior-cingulate', 'caudate', 'cerebellar-hemisphere', 'frontal-pole', \
                'hippocampus', 'hypothalamus', 'nucleus-accumbens', 'putamen', 'substantia-nigra']

############################################################################################################

## read in covariates 
cov_cols = [-4, -3, -2, -1] ## isFemale, age, PC1, PC2 
COV_DATA = np.loadtxt(cov_file, delimiter='\t', skiprows=1, usecols=cov_cols)

## if UKB: get covariates for cohort only 
if dset == 'UKB':
    subj_order1 = '/data1/rubinov_lab/brain_genomics/data_UKB/UKB_subjs_2732.txt' ## covariate order 
    subj_order2 = '/data1/rubinov_lab/brain_genomics/data_UKB/UKB_subjs_941.txt' ## cohort order  

    order1 = np.loadtxt(subj_order1) 
    order2 = np.loadtxt(subj_order2) 
    subj_idx = np.where(np.in1d(order1, order2))[0] ## indices of order2 wrt to order1
    COV_DATA = COV_DATA[subj_idx] 

## function: regress out all covariates 
## X: (n_samples, n_features)
## y: (n_samples) 
def compute_y_residuals(y): 
    lr = LinearRegression().fit(COV_DATA, y) 
    coefs = lr.coef_ ## (n_features) 
    b = lr.intercept_ ## float 
    ccovs = np.matmul(COV_DATA, coefs) ## (n_samples) 
    y_residuals = y - ccovs - b ## (n_samples) 
    return y_residuals  

############################################################################################################

## function: compute and save expression residuals 
def save_expr_residuals():
    for reg in phen_regions: 

        ## read original expression data (and genes) 
        with h5py.File('{}/{}.hdf5'.format(expr_dir, reg), 'r') as f: 
            expr_old = np.array(f['pred_expr']) ## (n_genes, n_samples)
            gene_old = np.array(f['genes']).astype(str) ## (n_genes) 

        ## compute residuals per gene 
        ## remove version number from gene names 
        expr_new = np.zeros_like(expr_old) 
        for g in range(gene_old.size): 
            expr_new[g] = compute_y_residuals(expr_old[g])
            gene_old[g] = gene_old[g].split('.')[0]

        ## save new expression data and genes  
        with h5py.File('{}/{}.hdf5'.format(expr_out, reg), 'w') as f: 
            f['pred_expr'] = expr_new 
            f['genes'] = gene_old.astype(bytes) 

## function: compute and save phenotype residuals 
def save_phen_residuals():
    phens = ['gm_volume', 'myelination', 'alff', 'reho_noGS', 'connmean_noGS']
    for phen in phens:

        ## read original phenotype values 
        phen_old = {} ## k: region, v: (n_samples) 
        with h5py.File('{}/{}.hdf5'.format(phen_dir, phen), 'r') as f: 
            for reg in phen_regions: 
                phen_old[reg] = np.array(f[reg])

        ## compute residuals per regional phenotype  
        phen_new = {} ## k: region, v: (n_samples) 
        for reg in phen_regions: 
            phen_new[reg] = compute_y_residuals(phen_old[reg])

        ## save new phenotype values
        with h5py.File('{}/{}.hdf5'.format(phen_out, phen), 'w') as f: 
            for reg in phen_regions:
                f[reg] = phen_new[reg]

############################################################################################################

save_expr_residuals()
save_phen_residuals() 

