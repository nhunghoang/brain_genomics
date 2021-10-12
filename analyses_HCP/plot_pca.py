'''
Script for generating PCA plots for the HCP 890 cohort. 
- genotype plot: based on first two PCs, as computed via EIGENSTRAT
- expression plot: based on first two PCs, only using filtered genes (regions-concatenate) 
- show for train and test groups  

Plot description: 
> 2 PCs 
> edge colors (white: red, black: green, Asian: blue, 
               mixed family reporting: yellow, Am Indian/Alaskan Nat: magenta) 
> fill colors (cyan if reported Hisp/Lat, none otherwise)

Nhung, updated Oct 2021
'''

import matplotlib.pyplot as plt
import numpy as np
import os
from sklearn.decomposition import PCA
import sys
import h5py

## paths 
expr_dir = '/data1/rubinov_lab/brain_genomics/data_HCP/expression/filtered_quality_r0.3_p0.01/cohort890' 
dem_file = '/data1/rubinov_lab/brain_genomics/data_HCP/subject_demographics/cohort890_demographics.txt'

spl_file = '/data1/rubinov_lab/brain_genomics/analyses_HCP/DATA_OUTPUT/train_test_assoc_split.hdf5' 
cov_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/DATA_OUTPUT' 
png_file = '/data1/rubinov_lab/brain_genomics/analyses_HCP/DATA_OUTPUT/pca_split.png'

## get train/test group indices 
with h5py.File(spl_file, 'r') as f: 
    train_idx = np.array(f['train_idx_890'])
    testt_idx = np.array(f['test_idx_890'])
    idxs = [train_idx, testt_idx] 

## parse demographics  
with open(dem_file, 'r') as f: 

    ## subjID sampID famID gender age race ethnicity 
    info = np.loadtxt(dem_file, delimiter='\t', skiprows=1, dtype=str)
    all_subjs = info[:,0] 
    all_race = info[:,-2]
    all_ethnicity = info[:,-1]

## define marker colors 
edge_dict = {'White':'r', 
             'Black or African Am.':'g',
             'Asian/Nat. Hawaiian/Othr Pacific Is.':'b',
             'Unknown or Not Reported': 'k',   
             'Am. Indian/Alaskan Nat.': 'm', 
             'More than one': 'y'}

## initialize PC plot
fig, axes = plt.subplots(2,2,figsize=(30,30))

## loop through train/test 
for g,group in enumerate(['train', 'test']):

    ## get indices of the group 
    idx = idxs[g] 
    subjs = all_subjs[idx] 

    ## group race/ethnicity 
    race = all_race[idx]
    ethnicity = all_ethnicity[idx] 

    ## read genotype PCs 
    cov_file = '{}/{}_covariates.txt'.format(cov_path, group)
    geno_pcs = np.loadtxt(cov_file, delimiter='\t', skiprows=1, usecols=[-2,-1])
    geno_pc1 = geno_pcs[:,0]
    geno_pc2 = geno_pcs[:,1] 

    ## compute expression PCs
    all_expr = [] 
    for regf in os.listdir(expr_dir): 
        if regf[-4:] != 'hdf5': continue  
        with h5py.File('{}/{}'.format(expr_dir, regf), 'r') as f: 
            reg_expr = np.array(f['pred_expr']).T[idx]
            all_expr.append(reg_expr)
    cat_expr = np.concatenate(all_expr, axis=1) ## (people * expr) 

    model = PCA(n_components=2)
    model.fit(cat_expr)
    expr_pcs = model.transform(cat_expr)
    expr_pc1 = expr_pcs[:,0]
    expr_pc2 = expr_pcs[:,1]

    ### plot PCs 
    titles = ['PCA on Genotypes ({})'.format(group), 'PCA on Expression ({})'.format(group)]
    edge_colors = [edge_dict[r] for r in race]
    fill_colors = ['c' if e=='Hispanic/Latino' else 'none' for e in ethnicity] 
    for i, pcs in enumerate([geno_pcs, expr_pcs]):
        pc1 = pcs[:,0]
        pc2 = pcs[:,1]
        ax = axes[g,i] 
        ax.scatter(pc1, pc2, s=200, color=fill_colors, edgecolors=edge_colors, linewidth=3)
        
        ax.set_title(titles[i], fontsize=30, fontweight='bold')
        ax.set_xlabel('1st PC', fontsize=30)
        ax.set_ylabel('2nd PC', fontsize=30)
        
        #ax.set_xlim(-20,50)
        #ax.set_ylim(-20,50)
        
        xleft, xright = ax.get_xlim()
        ybottom, ytop = ax.get_ylim()
        ax.set_aspect(abs((xright-xleft)/(ybottom-ytop)))
        ax.tick_params(labelsize=30)
    
#plt.savefig(PLOT_FILE, format='svg', dpi=1200)
plt.savefig(png_file)
plt.close('all')
