import numpy as np 
from scipy.stats import pearsonr
import h5py 
import os 

mat_dir = '/data1/rubinov_lab/brain_genomics/data_HCP_neuro' 
mat_files = os.listdir(mat_dir) 

vf = open('valid_timeseries.txt', 'w')  

movt_rrms = np.ones((1206,4), dtype=float) * (100)
for i,mat_file in enumerate(mat_files): 

    bad_scans = 0

    f = h5py.File('{}/{}'.format(mat_dir,mat_file), 'r') 
    rg = f['regressors'] 
    for j in range(4): 
        try: 
            mvt = np.array(f[rg[0][j]]['Movt_RRMS'])
            mvt_max = np.max(mvt)
            movt_rrms[i][j] = mvt_max 
        except: 
            bad_scans += 1 
            continue 

    #if bad_scans > 0: print('{}: {}'.format(mat_file.split('.')[0], bad_scans))
    line = '{}\t{}\n'.format(mat_file.split('.')[0], 4-bad_scans)
    vf.write(line)

    f.close()

vf.close()

with h5py.File('movt_rrms.hdf5', 'w') as f:
    f['max-movt'] = movt_rrms 
