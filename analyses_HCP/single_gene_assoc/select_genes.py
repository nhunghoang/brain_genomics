'''
For a given phenotype, identify the genes that independently and 
significantly correlate with the phenotype across individuals. 

NOTE: the gene models that are being considered have already been 
filtered based on completeness or quality thresholds. 

- Nhung, updated June 2021 
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

phenotype = sys.argv[1] 

## association outputs 
assoc_dir = '/data1/rubinov_lab/brain_genomics/analyses_HCP/single_gene_assoc/genes_by_reg/{}'.format(phenotype)  
if not os.path.exists(assoc_dir): os.mkdir(assoc_dir) 

summary_file = '{}/summary_table.txt'.format(assoc_dir) 
if os.path.exists(summary_file): os.remove(summary_file) 

## phenotype 
phen_file = '/data1/rubinov_lab/brain_genomics/data_HCP/phenotypes/{}.hdf5'.format(phenotype)
subj_file = '/data1/rubinov_lab/brain_genomics/data_HCP/timeseries_order.hdf5'

with h5py.File(phen_file, 'r') as f: 
    phens = {reg: np.array(f[reg]) for reg in f.keys()}
with h5py.File(subj_file, 'r') as f: 
    subj_order_all = np.array([str(s) for s in np.array(f['subjects'])])

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
with open(samp_file, 'r') as f: 
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
subj_order = subj_order_all[subj_idx]
samp_order = samp_order_all[samp_idx] 
for reg in regions: 
    phens[reg] = phens[reg][subj_idx]
    exprs[reg] = exprs[reg][:,samp_idx]

## read genes per region that HCP doesn't have all SNPs for
'''
unkept_file = '/data1/rubinov_lab/brain_genomics/data_HCP/expression/all_genes_not_kept.txt'
with open(unkept_file, 'r') as f: lines = f.readlines()
unkept_models = {} ## k: region, v: list of genes
for line in lines:
    data = line.split('\t')
    reg = data[0]
    mod = data[1]
    try: unkept_models[reg].append(mod)
    except KeyError: unkept_models[reg] = [mod]
'''

## single gene associations
for reg in regions: 

    phen_array = phens[reg]
    gene_array = genes[reg]
    expr_matrx = exprs[reg]

    n_incomp = 0 
    n_plot = 0 

    ## gene loop starts here 
    data = np.zeros((gene_array.shape[0], 3))
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

    ## save data 
    with open('{}/pearsonr_{}.txt'.format(assoc_dir, reg), 'w') as all_gene_corrs:
        header = '\t'.join(['GENE', 'RHO', 'PVAL', 'FDR', '\n'])
        all_gene_corrs.write(header) 
        for g,d in zip(gene_array,data): 
            line = '{}\t{:.5f}\t{:.5f}\t{:.5f}\n'.format(g, d[0], d[1], d[2])
            all_gene_corrs.write(line) 

    with open('{}/p-selected_{}.txt'.format(assoc_dir, reg), 'w') as p_sel: 
        genes_psig = gene_array[data[:,1] <= 0.05]
        for gp in genes_psig: 
            gp = gp.split('.')[0]
            p_sel.write(gp + '\n')

    with open('{}/fdr-selected_{}.txt'.format(assoc_dir, reg), 'w') as f_sel: 
        genes_fsig = gene_array[data[:,2] <= 0.05]
        for gf in genes_fsig: 
            gf = gf.split('.')[0]
            f_sel.write(gf + '\n')

    ## count number of significant genes with missing SNPs 
    '''
    p_incomp = np.intersect1d(genes_psig, unkept_models[reg])
    f_incomp = np.intersect1d(genes_fsig, unkept_models[reg])

    with open('{}/{}.log'.format(expr_dir, reg), 'r') as df:
        data_header = df.readline()
        lines = df.readlines() 
    data_lines = {} 
    for line in lines: 
        data_lines[line.split('\t')[0]] = line 
    with open('{}/incomp_{}.txt'.format(assoc_dir, reg), 'w') as cf: 
        cf.write(data_header)
        for c in p_incomp:
            cf.write(data_lines[c])
    '''

    ## plot some FDR genes 
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
    no_vars = np.count_nonzero(np.isnan(data[:,0]))
    p_sig = np.count_nonzero(data[:,1] <= 0.05)
    f_sig = np.count_nonzero(data[:,2] <= 0.05) 
    print('TOTAL  : {}'.format(gene_array.shape[0]))
    #print('NO VAR : {}'.format(no_vars))
    print('P SIG  : {}'.format(p_sig))
    print('FDR SIG: {}'.format(f_sig))
    #print('P SIG  : {} ({} incomp)'.format(p_sig, p_incomp.shape[0]))
    #print('FDR SIG: {} ({} incomp)'.format(f_sig, f_incomp.shape[0]))

    ## save summary 
    with open(summary_file, 'a') as sf: 
        info = [reg, str(gene_array.shape[0]), str(no_vars), str(p_sig), str(f_sig)]
        #info = [reg, str(gene_array.shape[0]), str(no_vars), str(p_sig), \
        #        str(p_incomp.shape[0]), str(f_sig), str(f_incomp.shape[0])]
        line = '\t'.join(info)
        sf.write(line + '\n') 
