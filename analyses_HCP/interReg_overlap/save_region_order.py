'''
Create an ordered region list for heatmap visualization based on 
running BCT reorder_mod on the subject-average co-activity matrix. 

- Nhung, Feb 2022 
'''

import os 
import h5py 
import bct 
import numpy as np 

## gather regions 
regs = os.listdir('/data1/rubinov_lab/brain_genomics/analyses_HCP/lofo_assoc/train_pvals_alff')
regs = np.array(sorted(regs)) ## abc order

## gather all co-activity values 
coact_path = '/data1/rubinov_lab/brain_genomics/data_HCP/2021_hoacer_sn_hth/phenotypes/coactivity.hdf5' 
coacts = {} ## k: (reg1,reg2, v: coact value 
with h5py.File(coact_path, 'r') as f: 
    for key in f.keys():
        k = tuple(key.split('_'))
        coacts[k] = np.array(f[key]).mean()

## populate co-activity matrix (diagonal is zero'd out)  
coact_mat = np.zeros((10,10), dtype=float)
for i, reg1 in enumerate(regs):
    for ii, reg2 in enumerate(regs[i+1:]):
        j = i + ii + 1
        coact_mat[i][j] = coacts[(reg1,reg2)] 
        coact_mat[j][i] = coacts[(reg1,reg2)] 

#out_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/interReg_overlap/avg_coact_mat.hdf5'
#with h5py.File(out_path, 'w') as f: 
#    f['avg_coact'] = coact_mat
#import sys; sys.exit()


## reorder matrix  
#c = np.arange(10)
#new_order, new_mat = bct.reorder_mod(coact_mat, c) 
new_order = np.array([10,7,1,6,2,5,4,3,9,8]) - 1 
out_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/interReg_overlap/region_order.txt'
with open(out_path, 'w') as f: 
    for reg in regs[new_order]: 
        f.write(reg + '\n') 


