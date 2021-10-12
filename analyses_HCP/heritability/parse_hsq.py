'''
'''

import sys
import os

tag = sys.argv[1] 
phenotype = sys.argv[2]

#try:
#    version = sys.argv[2] 
#    if version == '7': results_dir += '_v7' 
#except: 
#    pass 

regions = ['amygdala', 'anterior-cingulate', 'caudate', 'cerebellar-hemisphere', 'frontal-pole', 'hippocampus', 'hypothalamus', 'nucleus-accumbens', 'putamen', 'substantia-nigra']
#regions = ['amygdala', 'anterior-cingulate', 'caudate', 'cerebellar-hemisphere', 'frontal-pole', 'hippocampus', 'nucleus-accumbens', 'putamen']
#regions = ['amygdala', 'caudate', 'cerebellar-hemisphere', 'hippocampus', 'hypothalamus', 'nucleus-accumbens', 'putamen', 'substantia-nigra']

'''
print('ALL')
results_dir = '/data1/rubinov_lab/brain_genomics/analyses_HCP/heritability/results_gcta_{}/all'.format(tag)
for region in regions:
    result = '{}/{}_{}.hsq'.format(results_dir, phenotype, region)
    try:
        with open(result, 'r') as f: lines = f.readlines()
    except: 
        print('{:>22}: error'.format(region))
        continue
    H = float(lines[4].split('\t')[1])
    S = float(lines[4].split('\t')[2])
    print('{:>22}: {:>.6f} ± {}'.format(region, H, S))
'''

print('PDX')
results_dir = '/data1/rubinov_lab/brain_genomics/analyses_HCP/heritability/results_gcta_{}/pdx'.format(tag)
for region in regions:
    result = '{}/{}_{}.hsq'.format(results_dir, phenotype, region)
    try:
        with open(result, 'r') as f: lines = f.readlines()
    except: 
        print('{:>22}: error'.format(region))
        continue
    H = float(lines[4].split('\t')[1])
    S = float(lines[4].split('\t')[2])
    print('{:>22}: {:>.3f} ± {:>.3f}'.format(region, H, S))
