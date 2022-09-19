'''
Gather null data into uniform files. 

- Nhung, updated July 2022
'''

import numpy as np 
import h5py 
import sys 
from time import time 

dset = sys.argv[1] 

phenotypes = ['gm_volume', 'alff', 'reho_noGS', 'connmean_noGS', 'myelination']
regions = ['hippocampus', 'amygdala', 'hypothalamus', 'substantia-nigra',\
        'caudate', 'putamen', 'nucleus-accumbens', 'anterior-cingulate',\
        'frontal-pole', 'cerebellar-hemisphere']

## paths 
path_expr = '/data1/rubinov_lab/brain_genomics/analyses_{}/DATA_OUTPUT/r0.3_p0.01/expr_regress'.format(dset) 
path_null = '/data1/rubinov_lab/brain_genomics/analyses_{}/assoc_1M/nulls/pvals'.format(dset)  

nperms = int(1e4) 
run_start = time()

for phen in phenotypes: 
    for reg in regions: 

        loop_start = time() 

        ## get number of regional genes 
        efile = '{}/{}.hdf5'.format(path_expr, reg) 
        with h5py.File(efile, 'r') as f: ngenes = f['genes'].size

        ## init null matrix 
        null_data = np.empty((nperms, ngenes, 3))

        ## load permutations 
        for perm in range(nperms): 
            pfile = '{}_{}/{}_{}.hdf5'.format(path_null, phen, reg, perm) 
            try: 
                with h5py.File(pfile, 'r') as f: 
                    pdata = np.array(f['pearson'])
                null_data[perm] = pdata 
            except: 
                print(pfile)

        ## save collective permutation data 
        nfile = '{}_{}/{}.hdf5'.format(path_null, phen, reg) 
        with h5py.File(nfile, 'w') as f: 
            f['null_pearson'] = null_data 

        secs = time() - loop_start
        hr = int(secs//3600)
        mn = int((secs%3600)//60)
        sc = int((secs%3600)%60)
        print('{:d} hr, {:d} mn, {:d} sc | {} {}\n'.format(hr, mn, sc, phen, reg))

secs = time() - run_start 
hr = int(secs//3600)
mn = int((secs%3600)//60)
sc = int((secs%3600)%60)
print('total runtime: {:d} hr, {:d} mn, {:d} sc'.format(hr, mn, sc))
            
