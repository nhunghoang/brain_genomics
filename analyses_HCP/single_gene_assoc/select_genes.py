'''
For a given phenotype, identify the genes that independently and 
significantly correlate with the phenotype across individuals. 

NOTE: the gene models that are being considered have already been 
filtered based on completeness or quality thresholds. 

[07.01.21] Update: Record the r^2 value for each gene (i.e., rho^2)  

- Nhung, updated July 2021 
'''

import matplotlib 
matplotlib.use('Agg')
import matplotlib.pyplot as plt 
import sys 
import os 
import numpy as np  
import h5py 
from scipy.stats import pearsonr, spearmanr 
from statsmodels.stats.multitest import fdrcorrection 

atlas = sys.argv[1] 
phenotype = sys.argv[2] 

## association outputs 
assoc_dir = '/data1/rubinov_lab/brain_genomics/analyses_HCP/single_gene_assoc/genes_by_phen/{}'.format(phenotype)  
if not os.path.exists(assoc_dir): os.mkdir(assoc_dir) 

## phenotype 
phen_file = '/data1/rubinov_lab/brain_genomics/data_HCP/{}/phenotypes/{}.hdf5'.format(atlas, phenotype)
with h5py.File(phen_file, 'r') as f: 
    phens = {reg: np.array(f[reg]) for reg in f.keys()}
regions = np.array(list(phens.keys()))

## expression 
expr_dir = '/data1/rubinov_lab/brain_genomics/data_HCP/expression/filtered_quality_r0.3_p0.01'
samp_file = '/data1/rubinov_lab/brain_genomics/data_HCP/predixcan_v8/sample_ids.txt' 

genes = {}; exprs = {} 
for expr_file in os.listdir(expr_dir): 
    if expr_file[-5:] != '.hdf5': continue 
    with h5py.File('{}/{}'.format(expr_dir, expr_file), 'r') as f: 
        reg = expr_file.split('.')[0]
        genes[reg] = np.array([g.decode("utf-8") for g in np.array(f['genes'])])
        exprs[reg] = np.array(f['pred_expr']) ## (genes * samps)

## subject-sample ordering 
order_file = '/data1/rubinov_lab/brain_genomics/analyses_HCP/subj_samp_assoc_order.hdf5'
with h5py.File(order_file, 'r') as f: 
    samp_order = np.array(f['sample_idx_1142']) 
for reg in regions: 
    exprs[reg] = exprs[reg][:,samp_order]

## single gene associations
for reg in regions: 

    phen_array = phens[reg]
    gene_array = genes[reg]
    expr_matrx = exprs[reg]

    ## BANDAID 
    if phenotype=='myelination' and reg=='frontal-pole': 
        continue 

    ## BANDAID 
    if phenotype in ['fa', 'md']: 
        miss_path = '/data1/rubinov_lab/brain_genomics/data_HCP/hoacer_sn_hth/phenotypes/famd_missing_idx.txt'
        with open(miss_path, 'r') as f: lines = f.readlines()
        miss_idx = [int(i) for i in lines]
        good_idx = np.setdiff1d(np.arange(890), miss_idx)
        expr_matrx = expr_matrx[:,good_idx]

    ## gene loop starts here 
    data = np.zeros((gene_array.shape[0], 3)) ## rho, pval, fdr
    for g,gene in enumerate(gene_array): 
        gene_expr = expr_matrx[g]
        
        ## return nan for genes with no expression variance  
        if np.var(gene_expr) == 0: 
            data[g,:] = np.nan
            continue  

        ## otherwise, compute correlation 
        rho, pval = pearsonr(phen_array, gene_expr) 
        data[g,0] = rho; data[g,1] = pval 

    ## compute FDR-correct pvals 
    pvals = data[:,1]
    pvals = pvals[~np.isnan(pvals)]
    rejected, corrected = fdrcorrection(pvals) 
    c = 0
    for g in range(gene_array.shape[0]): 
        if np.isnan(data[g,2]): continue 
        else: data[g,2] = corrected[c]; c += 1 

    ## save data (gene, r^2, rho, pval, fdr)
    with open('{}/pearsonr_{}.txt'.format(assoc_dir, reg), 'w') as all_gene_corrs:
        header = '\t'.join(['GENE', 'R^2', 'RHO', 'PVAL', 'FDR', '\n'])
        all_gene_corrs.write(header) 
        for g,d in zip(gene_array,data): 
            line = '{}\t{:.5f}\t{:.5f}\t{:.5f}\t{:.5f}\n'.format(g, d[0], d[0]**2, d[1], d[2])
            all_gene_corrs.write(line) 

    #with open('{}/p-selected_{}.txt'.format(assoc_dir, reg), 'w') as p_sel: 
    #    genes_psig = gene_array[data[:,1] <= 0.05]
    #    r2_values = data[:,1][data[:,1] <= 0.05]
    #    for gp, r2 in zip(genes_psig, r2_values): 
    #        gp = gp.split('.')[0]
    #        p_sel.write('{}: {:.3f}\n'.format(gp, r2))
    #        #p_sel.write(gp + '\n')
    #    print('{:>25s} (p): ({:.3f}, {:.3f})'.format(reg, r2_values.min(), r2_values.max()))

    #with open('{}/fdr-selected_{}.txt'.format(assoc_dir, reg), 'w') as f_sel: 
    #    genes_fsig = gene_array[data[:,2] <= 0.05]
    #    r2_values = data[:,2][data[:,2] <= 0.05]
    #    for gf, r2 in zip(genes_fsig, r2_values): 
    #        gf = gf.split('.')[0]
    #        f_sel.write('{}: {:.3f}\n'.format(gf, r2))
    #        #f_sel.write(gf + '\n')
    #    if genes_fsig.size != 0: 
    #        print('{:>25s} (f): ({:.3f}, {:.3f})'.format(reg, r2_values.min(), r2_values.max()))

    ## plot some FDR genes 
    n_plot = 0 
    genes_fsig = gene_array[data[:,2] <= 0.05]
    for gf in genes_fsig: 
        if n_plot > 5: break 
        #if gf in unkept_models[reg]: continue  
        n_plot += 1 
        gf_idx = np.argwhere(genes[reg] == gf)[0][0]
        gf_expr = expr_matrx[gf_idx] 
        rho, pval = pearsonr(gf_expr, phen_array) 
        b, m = np.polynomial.polynomial.polyfit(gf_expr, phen_array, 1) 

        fig, ax = plt.subplots(1,1,figsize=(15,15))
        ax.scatter(gf_expr, phen_array, c='k')
        ax.plot(gf_expr, b+(m*gf_expr), '-', c='y') 
        ax.set_xlabel('expression', fontsize=30)
        ax.set_ylabel(phenotype, fontsize=30) 
        ax.tick_params(labelsize=30)
        
        xleft, xright = ax.get_xlim()
        ybottom, ytop = ax.get_ylim()
        ax.set_aspect(abs((xright-xleft)/(ybottom-ytop)))
    
        title = '{}, {}\nr={:.3f} (p ≤ {:.3f})'.format(gf, reg, rho, pval)
        ax.set_title(title, size=30)
        fname = '{}/{}_{}_{:.3f}.png'.format(assoc_dir, reg, gf, rho) 
        plt.savefig(fname)  
        plt.close('all') 

    ## print summary 
    print('\n- {} -'.format(reg))
    #no_vars = np.count_nonzero(np.isnan(data[:,0]))
    p_sig = np.count_nonzero(data[:,1] <= 0.05)
    f_sig = np.count_nonzero(data[:,2] <= 0.05) 
    print('TOTAL  : {}'.format(gene_array.shape[0]))
    #print('NO VAR : {}'.format(no_vars))
    print('P SIG  : {}'.format(p_sig))
    print('FDR SIG: {}'.format(f_sig))
    #print('P SIG  : {} ({} incomp)'.format(p_sig, p_incomp.shape[0]))
    #print('FDR SIG: {} ({} incomp)'.format(f_sig, f_incomp.shape[0]))
