'''
Concatentate the TWAS and GWAS results into one big table. 

- Nhung, updated Feb 2024
'''

import pandas as pd  
import numpy as np
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests as sm

import sys 
import os

from time import time 
from multiprocessing import Pool

## params 
group = sys.argv[1] ## UKB
phens = sys.argv[2] ## vol_mean

regs = ['dlpfc', 'anterior-cingulate', 'caudate', 'nucleus-accumbens', \
        'putamen', 'amygdala', 'hippocampus', 'cerebellar-hemisphere']

## paths 
#main_path = '/Users/nhunghoang/Desktop/remote_platypus'
main_path = '/data1/rubinov_lab/brain_genomics'
wjti_path = f'{main_path}/models_JTI/weights_by_tissue'

twas_path = f'{main_path}/paper_twas/outputs_{group}/twas_JTI/{phens}' 
gwas_path = f'{main_path}/paper_twas/outputs_{group}/gwas/{phens}'

## out paths 
summ_path = f'{main_path}/paper_twas/outputs_{group}/TWAS_GWAS_summary.csv'
sigs_path = f'{main_path}/paper_twas/outputs_{group}/TWAS_GWAS_sigbars.csv'

summ_rtmp = f'{main_path}/paper_twas/outputs_{group}/tgwa' 
sigs_rtmp = f'{main_path}/paper_twas/outputs_{group}/sigs' 

################################################################################

## get TWAS and GWAS data
def get_data(reg):
    
    print(f'starting on {reg}...')

    ## phen, GWAS_FDR, GWAS_BON, TWAS_FDR, GWAS_BON
    sig_data = {'phenotype': [reg]}

    ## read GWAS
    gwas_file = f'{gwas_path}_{reg}.regenie'
    gwas_data = pd.read_table(gwas_file, usecols=['ID', 'BETA', 'LOG10P'], sep=' ')
    gwas_data = gwas_data.rename(columns={'ID': 'rsid', 'BETA': 'beta', 'LOG10P': '-log10p'})

    pvals = gwas_data['-log10p']
    gwas_data['pval'] = np.power(10, -pvals)

    fkeep, fvals, _, _ = sm(gwas_data['pval'], method='fdr_bh', alpha=0.05)
    bkeep, bvals, _, alpha = sm(gwas_data['pval'], method='bonferroni', alpha=0.05)
    
    gwas_data['FDR'] = fvals
    gwas_data['BON'] = bvals  

    sig_data['GWAS_FDR'] = [gwas_data['pval'][fkeep].max()]
    sig_data['GWAS_BON'] = [alpha]

    ## read TWAS
    twas_file = f'{twas_path}/{reg}.txt'
    twas_data = pd.read_table(twas_file, usecols=['gene', 'pvalue', 'effect'])
    twas_data = twas_data.rename(columns={'pvalue': 'pval', 'effect': 'beta'})

    fkeep, fvals, _, _ = sm(twas_data['pval'], method='fdr_bh', alpha=0.05)
    bkeep, bvals, _, alpha = sm(twas_data['pval'], method='bonferroni', alpha=0.05)
    
    twas_data['FDR'] = fvals
    twas_data['BON'] = bvals  

    sig_data['TWAS_FDR'] = [twas_data['pval'][fkeep].max()]
    sig_data['TWAS_BON'] = [alpha]

    gwas_data = gwas_data[['rsid', 'beta', 'pval', 'FDR', 'BON']]
    twas_data = twas_data[['gene', 'beta', 'pval', 'FDR', 'BON']]
    return gwas_data, twas_data, sig_data 

################################################################################

## set region-specific tables 
def set_tables(reg): 

    rstart = time()

    ## get table of gene-snp pairings 
    mapp = pd.read_csv(f'{wjti_path}/{reg}.csv', usecols=['rsid', 'gene'])

    ## get TWAS/GWAS data 
    gwas, twas, stable = get_data(reg) 

    ## set sig table
    stable = pd.DataFrame(stable)

    rtime = time() - rstart
    (hr, mn, sc) = (rtime//3600, (rtime%3600)//60, (rtime%3600)%60)
    print('{:d} hr {:d} mn {:d} sc to get data ({})'.format(int(hr), int(mn), int(sc), reg))

    ## merge twas table with mapp (this should add rsids to the table) 
    rstart = time()
    mtable = twas.merge(mapp, on='gene', how='inner')

    ## merge gwas table with the new table 
    rtable = mtable.merge(gwas, on='rsid', how='inner', suffixes=['_TWAS', '_GWAS'])

    ## make sure we have all the TWAS genes 
    gene_left = np.setdiff1d(twas['gene'], rtable['gene'])
    twas_left = twas.loc[twas['gene'].isin(gene_left)]
    twas_left = twas_left.rename(columns={'beta': 'beta_TWAS', 'pval': 'pval_TWAS', \
                                          'FDR': 'FDR_TWAS', 'BON': 'BON_TWAS'})

    rtable = rtable.append(twas_left)
    rtable['phenotype'] = reg

    print(rtable.loc[rtable['FDR_TWAS'] <= 0.05]['gene'].drop_duplicates().shape, reg)

    ## merge time
    rtime = time() - rstart
    (hr, mn, sc) = (rtime//3600, (rtime%3600)//60, (rtime%3600)%60)
    print('{:d} hr {:d} mn {:d} sc to merge ({})'.format(int(hr), int(mn), int(sc), reg))

    ## save tables
    rstart = time()

    stable.to_csv(f'{sigs_rtmp}_{reg}.csv', index=False)
    rtable.to_csv(f'{summ_rtmp}_{reg}.csv', index=False)

    rtime = time() - rstart
    (hr, mn, sc) = (rtime//3600, (rtime%3600)//60, (rtime%3600)%60)
    print('{:d} hr {:d} mn {:d} sc to write table ({})'.format(int(hr), int(mn), int(sc), reg))

################################################################################

## compare two gene sets 
def venn(a, b):

    a_not_b = np.setdiff1d(a, b).size
    b_not_a = np.setdiff1d(b, a).size
    a_and_b = np.intersect1d(a, b).size

    vd = '( {} ( {} ) {} )'.format(a_not_b, a_and_b, b_not_a)
    return vd

################################################################################

## main call
#pool = Pool(processes=len(regs))
#pool.map(set_tables, regs)
#
### concat regional TWAS/GWAS summaries
#r0 = regs[0]
#cmd = f'cp {summ_rtmp}_{r0}.csv {summ_path}'
#os.system(cmd)
#for reg in regs[1:]: 
#    cmd = f'tail -n +2 {summ_rtmp}_{reg}.csv >> {summ_path}'
#    os.system(cmd)
#
### concat regional sigs summaries
#cmd = f'cp {sigs_rtmp}_{r0}.csv {sigs_path}'
#os.system(cmd)
#for reg in regs[1:]: 
#    cmd = f'tail -n +2 {sigs_rtmp}_{reg}.csv >> {sigs_path}'
#    os.system(cmd)

df = pd.read_csv(summ_path)

rdfs = df.groupby('phenotype')
for reg, rdf in rdfs:

    ## venn diagram
    tpass = rdf.loc[rdf['BON_TWAS'] <= 0.05]['gene'].drop_duplicates()
    gpass = rdf.loc[rdf['BON_GWAS'] <= 0.05]['gene'].drop_duplicates()
    spass = rdf.loc[rdf['BON_GWAS'] <= 0.05]['rsid'].drop_duplicates()

    print(reg)
    print(venn(tpass, gpass))

    ### stats
    print('{} TWAS genes\n'.format(tpass.shape))
    print('{} GWAS genes\n'.format(gpass.shape))
    print('{} GWAS SNPs\n'.format(spass.shape))

    ## keep the best GWAS for every (TWAS) gene
    pdf = rdf.dropna().sort_values('pval_GWAS', ascending=True).drop_duplicates('gene')

    ## case: multiple genes share the same best GWAS
    ## keep the best TWAS gene for every SNP
    tdf = pdf.sort_values('pval_TWAS', ascending=True).drop_duplicates('rsid')

    print(reg)
    rho, pvl = spearmanr(tdf['pval_TWAS'], tdf['pval_GWAS'])
    print('top: {:.3f} ({})'.format(rho, pvl))

    rho, pvl = spearmanr(rdf.dropna()['pval_TWAS'], rdf.dropna()['pval_GWAS'])
    print('all: {:.3f} ({})\n'.format(rho, pvl))
