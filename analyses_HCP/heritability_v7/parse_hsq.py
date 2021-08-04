'''
'''

import sys
import os

results_dir = '/data1/rubinov_lab/brain_genomics/analyses_HCP/heritability_v7/results_gcta'

phenotype = sys.argv[1]

regions = ['amygdala', 'hippocampus', 'putamenBG', 'aCingulateCortex', 'frontalCortex', 'nAccumbensBG', 'caudateBG', 'cerebellum']

for region in regions:
    result = '{}/{}_{}.hsq'.format(results_dir, phenotype, region)
    with open(result, 'r') as f: lines = f.readlines()
    H = lines[4].split('\t')[1]
    print('{:>22}: {}'.format(region, H))
