'''
Figure 1D: PCA of gr-expression for HCP cohort (& legend) 

- Nhung, updated Sept 2022 
'''

import matplotlib 
matplotlib.use('tkagg')
import matplotlib.pyplot as plt 
import matplotlib.patches as mpatches 
import numpy as np 
import h5py 
from sklearn.decomposition import PCA
import os 

expr_dir = '/Users/nhunghoang/Desktop/local_platypus/expr_noregr' 
dem_file = '/Users/nhunghoang/Desktop/local_platypus/cohort890_demographics.txt'
cov_file = '/Users/nhunghoang/Desktop/local_platypus/covariates.txt'

## parse demographics  
with open(dem_file, 'r') as f: 

    ## subjID sampID famID gender age race ethnicity 
    info = np.loadtxt(dem_file, delimiter='\t', skiprows=1, usecols=[5,6], dtype=str)
    race = info[:,-2]
    ethnicity = info[:,-1]

## compute expression PCs
all_expr = [] 
for regf in os.listdir(expr_dir): 
    with h5py.File('{}/{}'.format(expr_dir, regf), 'r') as f: 
        reg_expr = np.array(f['pred_expr']).T
        all_expr.append(reg_expr)
cat_expr = np.concatenate(all_expr, axis=1)

model = PCA(n_components=2)
model.fit(cat_expr)
expr_pcs = model.transform(cat_expr)

## define marker plots 
edge_dict = {'White':'r', 
             'Black or African Am.':'g',
             'Asian/Nat. Hawaiian/Othr Pacific Is.':'b',
             'Unknown or Not Reported': 'k',   
             'Am. Indian/Alaskan Nat.': 'm', 
             'More than one': 'y'}
edge_colors = [edge_dict[r] for r in race]
fill_colors = np.where(ethnicity=='Hispanic/Latino', 'c', 'none')

## plot PCA
plt.ion()
fig, axes = plt.subplots(2,1, figsize=(3,6))

ax = axes[0]
kws = {'s':30, 'alpha':0.7, 'linewidth':0.7}

mask = np.array(ethnicity) == 'Hispanic/Latino'
a0 = ax.scatter(expr_pcs[:,0][mask], \
                expr_pcs[:,1][mask], \
                facecolor=fill_colors[mask], \
                edgecolor='none', \
                **kws)

a1 = ax.scatter(expr_pcs[:,0], \
                expr_pcs[:,1], \
                facecolor='none', \
                edgecolor=edge_colors, \
                **kws)
               
xleft, xright = ax.get_xlim()
ybottom, ytop = ax.get_ylim()
ax.set_aspect(abs((xright-xleft)/(ybottom-ytop)))

ax.set_xticks([-5, 0, 5, 10, 15])
ax.set_xticklabels(['-5', '', '', '', '15'], fontsize=10)

ax.set_yticks([-5, 0, 5, 10, 15])
ax.set_yticklabels(['-5', '', '', '', '15'], fontsize=10)

ax.set_xlabel('PC 1 (5.63% var)', size=10, weight='bold', loc='center', labelpad=-10)
ax.set_ylabel('PC 2 (2.10% var)', size=10, weight='bold', loc='center', labelpad=-14)

## legend 
ax = axes[1] 
ax.axis('off')

kws = {'radius':0.02, 'clip_on':False, 'facecolor':'none', 'alpha':0.7}
for y, ec in zip(range(95, 35, -10), ['g', 'b', 'r', 'm', 'y', 'c']):
    c = mpatches.Circle((0, 0.01*y), edgecolor=ec, **kws) 
    ax.add_patch(c) 
c.set_facecolor('c') 

pops = ['African', 'Asian', 'European', 'Nat. American', 'Mixed', 'Hispanic/Latino']
for y, txt in zip(range(93, 33, -10), pops): 
    ax.text(0.08, 0.01*y, txt, size=10)

plt.savefig('fig_plots/pca.pdf', format='pdf')
