'''
'''

import matplotlib.pyplot as plt 
import os 

path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/new_multi'

phens = ['alff', 'regional_homogeneity', 'gm_volume']
regions = os.listdir(path + '/results_alff')

fig, ax = plt.subplots(3,1,figsize(30,8))
for a in range(3):
