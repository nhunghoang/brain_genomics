'''
Count the overlap of gene associations 
between UKB volumes and BioVU phecodes. 


Cleaned up version - treat every gene by region-specificity.  

For all UKB volume genes that pass FDR, 
take that subset from the PrediXVU catalog. 

Then, FDR correct PrediXVU for those genes only. 
Then, find the overlap of genes that pass FDR
for both volume and clinical. 

Save the following information: 

volume, phecode, phename, phecat, phecode_count, phecat_count, gene, gene_tissue, gene_symbol

Just save this info without collapsing genes to lists. 

'''

import numpy as np 
import pandas as pd
from statsmodels.stats.multitest import multipletests as sm

import sys

####################################################################################

## params 
ptype = sys.argv[1] ## FDR, BON

## paths 
top_path = '/data1/rubinov_lab/brain_genomics'

sym_path = f'{top_path}/models_JTI/syms_by_tissue/jti_ens_to_sym.csv'
phe_path = f'{top_path}/paper_twas/aux_files/predixvu_phenames.csv'

pdx_path = f'{top_path}/paper_twas/aux_files/predixvu_assocs.csv.gz'
twa_path = f'{top_path}/paper_twas/outputs_UKB/twas_JTI/summary.csv'

## out paths 
mid_path = f'{top_path}/paper_twas/outputs_UKB/pdx_overlap/{ptype}_mdf.csv'
gxt_path = f'{top_path}/paper_twas/outputs_UKB/pdx_overlap/{ptype}_table_genxtis.csv'
pxv_path = f'{top_path}/paper_twas/outputs_UKB/pdx_overlap/{ptype}_table_phexvol.csv'

####################################################################################

## regions 
mod_names = ['Nucleus_accumbens_basal_ganglia', 'Hippocampus', 'Frontal_Cortex_BA9', \
             'Caudate_basal_ganglia', 'Cerebellar_Hemisphere', 'Amygdala', \
             'Anterior_cingulate_cortex_BA24', 'Putamen_basal_ganglia']

mod_names = [f'Brain_{m}_combo' for m in mod_names]
reg_names = ['nucleus-accumbens', 'hippocampus', 'dlpfc', 'caudate', \
             'cerebellar-hemisphere', 'amygdala', 'anterior-cingulate', 'putamen']

reg_map = {m:r for m,r in zip(mod_names, reg_names)}

####################################################################################

## get ensembl-symbol map
sym_map = pd.read_csv(sym_path, index_col='ens').to_dict()['sym']

## get phecodes and categories 
phe = pd.read_csv(phe_path, index_col='phecode').to_dict()
phe_map = phe['phenotype'] 
cat_map = phe['category']

####################################################################################

## func: FDR correction 
def fcorrect(tdf): 
    tdf['FDR'] = sm(tdf['pval'], method='fdr_bh', alpha=0.05)[1]
    return tdf 

## func: Bonf correction 
def bcorrect(tdf): 
    tdf['BON'] = sm(tdf['pval'], method='bonferroni', alpha=0.05)[1]
    return tdf 

####################################################################################

## func: read TWAS tables
def read_twas():
    pfunc = None
    if ptype == 'FDR': pfunc = fcorrect
    if ptype == 'BON': pfunc = bcorrect 

    ## read TWAS summary
    cols = {'phen': 'volume', 'grex': 'tissue', 'ens': 'ens', ptype: f'{ptype}_UKB'}
    twas = pd.read_csv(twa_path, usecols=cols.keys())
    twas = twas.rename(columns=cols)
    twas['sym'] = twas['ens'].map(sym_map)

    ## get tissue-gene of sig assocs
    twas = twas.loc[twas[f'{ptype}_UKB'] <= 0.05]
    twas_genes = twas[['tissue', 'ens']].drop_duplicates()
    print('read UKB TWAS')

    ## read PrediXVU catalog  
    cols = {'ensembl': 'ens', 'phecode': 'phecode', 'p-value': 'pval', 'tissue': 'tissue'}
    pdx = pd.read_csv(pdx_path, usecols=cols.keys()) 
    pdx = pdx.rename(columns=cols)
    pdx['tissue'] = pdx['tissue'].map(reg_map)
    print('read PDX TWAS')

    ## filter PDX based on UKB TWAS
    pdx1 = pdx.merge(twas_genes, on=['tissue', 'ens'], how='inner')
    print('filtered PDX TWAS')

    ## group PDX by (ensembl, tissue) and apply FDR/BON 
    pdx2 = pdx1.groupby(['phecode', 'tissue']).apply(pfunc)
    pdx2 = pdx2.rename(columns={ptype: f'{ptype}_PDX'})
    pdx3 = pdx2.loc[pdx2[f'{ptype}_PDX'] <= 0.05]
    print('corrected PDX pvals')

    ## merge UKB volume and PDX phecode results 
    mdf = pdx3.merge(twas, on=['tissue', 'ens'], how='inner')

    ## add phenames and phecats
    mdf['phename'] = mdf['phecode'].map(phe_map)
    mdf['phecate'] = mdf['phecode'].map(cat_map)

    cols = ['tissue', 'sym', 'ens', 'volume', f'{ptype}_UKB', \
            'phecode', f'{ptype}_PDX', 'phename', 'phecate']

    mdf = mdf[cols]
    mdf.to_csv(mid_path, index=False)
    return mdf

####################################################################################

mdf = read_twas()
mdf = pd.read_csv(mid_path)

####################################################################################

## add a 'tissue_gene' column to main table 
mdf['tgene'] = mdf['tissue'] + '_' + mdf['sym']

## table: phecat * volume => unique reg_genes
phexvol = mdf.groupby(['volume', 'phecate'])['tgene'].nunique()
phexvol = phexvol.reset_index(name='phecat_count')
phexvol = phexvol.pivot(index='phecate', columns='volume', values='phecat_count')
phexvol = phexvol.fillna(0)
phexvol.to_csv(pxv_path)

####################################################################################

## add a 'volume_phecat' column to main table 
mdf['volcat'] = mdf['volume'] + '_' + mdf['phecate']

## table: gene * tissue => instances in main table
genxtis = mdf.groupby(['tissue', 'sym'])['volcat'].nunique()
genxtis = genxtis.reset_index(name='volcat_count')
genxtis = genxtis.pivot(index='sym', columns='tissue', values='volcat_count')
genxtis = genxtis.fillna(0)
genxtis.to_csv(gxt_path)

####################################################################################
