'''
For a given phenotype, identify the genes that independently and 
significantly correlate with the phenotype across individuals. 

Notes/Updates: 
- The gene models that are being considered have already been 
  filtered based on completeness or quality thresholds. 
- Permutation testing is used to compute p-values, and twin 
  structures are maintained during the shuffling. 
- Input expression and phenotype data should already be in their 
  residual form (i.e. age, gender, PC1, PC2 confounders have been 
regressed out). 
- Associations are based on the training group, rather than the 
  entire group. 

- Nhung, updated Oct 2021 
'''

#import matplotlib 
# matplotlib.use('Agg')
#import matplotlib.pyplot as plt 
import sys 
import os 
import numpy as np  
import h5py 
from scipy.stats import pearsonr, spearmanr 
from statsmodels.stats.multitest import fdrcorrection 

#atlas = sys.argv[1] 
atlas = 'hoacer_sn_hth'
phenotype = sys.argv[1] 

train = True

## paths 
phn_file = '/data1/rubinov_lab/brain_genomics/analyses_HCP/phen_regress/{}.hdf5'.format(phenotype)
expr_dir = '/data1/rubinov_lab/brain_genomics/analyses_HCP/expr_regress'
twn_file = '/data1/rubinov_lab/brain_genomics/analyses_HCP/twin_indices.hdf5'
out_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/new_single/regress_{}'.format(phenotype)

if train:
    phn_file = '/data1/rubinov_lab/brain_genomics/analyses_HCP/phen_train_regress/{}.hdf5'.format(phenotype)
    expr_dir = '/data1/rubinov_lab/brain_genomics/analyses_HCP/expr_train_regress'
    twn_file = '/data1/rubinov_lab/brain_genomics/analyses_HCP/train_twin_indices.hdf5'
    out_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/new_single/regress_train_{}'.format(phenotype)

if not os.path.exists(out_path): os.mkdir(out_path) 

## get list of indices in order of (flattened twin array) and (non_twin array) 
## so it will be like [t1a, t1b, t2a, t2b, t3a, t3b..., nt1] 
with h5py.File(twn_file, 'r') as f: 
    twin_idx = np.array(f['twins_idx'])
    non_twin_idx = np.array(f['non_twins_idx'])
flat_twin = twin_idx.flatten()  
full_idx = np.concatenate((flat_twin, non_twin_idx))

## phenotype (arrange in full_idx order)  
with h5py.File(phn_file, 'r') as f: 
    phens = {reg: np.array(f[reg])[full_idx] for reg in f.keys()}
regions = np.array(list(phens.keys()))

## expression (arrange in full_idx order) 
genes = {}; exprs = {} 
for expr_file in os.listdir(expr_dir): 
    if expr_file[-5:] != '.hdf5': continue 
    with h5py.File('{}/{}'.format(expr_dir, expr_file), 'r') as f: 
        reg = expr_file.split('.')[0]
        genes[reg] = np.array([g.decode("utf-8") for g in np.array(f['genes'])])
        exprs[reg] = np.array(f['pred_expr'])[:,full_idx] ## (genes * samps)

## function: generate indices for one null permutation
def permute_indices(): 
    tidx = twin_idx * 1
    ntidx = non_twin_idx * 1
    np.random.shuffle(tidx)
    np.random.shuffle(ntidx) 
    return np.concatenate((tidx.flatten(), ntidx))

## function: compute p-value for specific gene association 
def compute_pvalue(expr, phen, rho): 
    N = 10000
    null = np.zeros(N) 
    for n in range(N): 
        nidx = permute_indices()
        null[n] = pearsonr(expr[nidx], phen)[0]
    pval = np.mean(np.abs(null) >= np.abs(rho))
    return pval 

## function: save gene-wise results 
## (gene name, rho, pval, fdr)  
def save_results(gene_array, data, filename):
    with open(filename, 'w') as f:
        header = '\t'.join(['GENE', 'RHO', 'PVAL', 'FDR', '\n'])
        f.write(header) 
        for g,d in zip(gene_array,data): 
            line = '{}\t{:.5f}\t{:.5f}\t{:.5f}\n'.format(g, d[0], d[1], d[2])
            f.write(line) 

## function: scatter-plot association (for a specific gene) 
#def plot_gene_result(): 
#    ## plot some FDR genes 
#    n_plot = 0 
#    genes_fsig = gene_array[data[:,2] <= 0.05]
#    for gf in genes_fsig: 
#        if n_plot > 5: break 
#        #if gf in unkept_models[reg]: continue  
#        n_plot += 1 
#        gf_idx = np.argwhere(genes[reg] == gf)[0][0]
#        gf_expr = expr_matrx[gf_idx] 
#        rho, pval = pearsonr(gf_expr, phen_array) 
#        b, m = np.polynomial.polynomial.polyfit(gf_expr, phen_array, 1) 
#
#        fig, ax = plt.subplots(1,1,figsize=(15,15))
#        ax.scatter(gf_expr, phen_array, c='k')
#        ax.plot(gf_expr, b+(m*gf_expr), '-', c='y') 
#        ax.set_xlabel('expression', fontsize=30)
#        ax.set_ylabel(phenotype, fontsize=30) 
#        ax.tick_params(labelsize=30)
#        
#        xleft, xright = ax.get_xlim()
#        ybottom, ytop = ax.get_ylim()
#        ax.set_aspect(abs((xright-xleft)/(ybottom-ytop)))
#    
#        title = '{}, {}\nr={:.3f} (p ≤ {:.3f})'.format(gf, reg, rho, pval)
#        ax.set_title(title, size=30)
#        fname = '{}/{}_{}_{:.3f}.png'.format(assoc_dir, reg, gf, rho) 
#        plt.savefig(fname)  
#        plt.close('all') 

## single gene associations
for reg in regions: 
    print('\n' + reg) 

    phen_array = phens[reg]
    gene_array = genes[reg]
    expr_matrx = exprs[reg]

    ngenes = gene_array.size

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
        pval = compute_pvalue(gene_expr, phen_array, rho)  
        data[g,0] = rho; data[g,1] = pval 

        ## status 
        if (g+1)%10 == 0: 
            perc = ((g+1)/ngenes) * 100 
            print('  {:>.2f}% [{}/{}]'.format(perc, g+1, ngenes))

    ## compute FDR-correct pvals 
    pvals = data[:,1]
    pvals = pvals[~np.isnan(pvals)]
    rejected, corrected = fdrcorrection(pvals) 
    c = 0
    for g in range(gene_array.shape[0]): 
        if np.isnan(data[g,2]): continue 
        else: data[g,2] = corrected[c]; c += 1 

    ## save results 
    filename = '{}/pearsonr_{}.txt'.format(out_path, reg) 
    save_results(gene_array, data, filename)

    ## print summary 
    p_sig = np.sum(data[:,1] <= 0.05)
    f_sig = np.sum(data[:,2] <= 0.05)
    print('> {:>3d} p, {:>3d} f'.format(p_sig, f_sig))
