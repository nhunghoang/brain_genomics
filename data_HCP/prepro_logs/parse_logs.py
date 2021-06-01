import numpy as np 
import h5py
import os 

with open('summary.log', 'r') as f: lines = f.readlines()

data = {} 
for line in lines: 
    subj = line.split(' ')[0]
    info = line.split(':')[1].strip() 
    try: 
        data[subj].append(info) 
    except KeyError: 
        data[subj] = [info] 

valid = []
regrs = [] 
qualc = [] 
bdata = [] 
for subj in data.keys(): 
    if 'no error' in data[subj]: valid.append(subj) 
    if 'no regressors' in data[subj]: regrs.append(subj)
    if 'qc_issue' in data[subj]: qualc.append(subj)
    if 'no data' in data[subj]: bdata.append(subj)

print('{} subjs have at least one ts'.format(len(valid)))
print('{} subjs have at least one regr issue'.format(len(regrs)))
print('{} subjs have at least one qc issue'.format(len(qualc)))
print('{} subjs have at least one data issue'.format(len(bdata)))

mat0_dir = '/data1/rubinov_lab/brain_genomics/data_HCP/timeseries'
mat1_dir = '/data1/rubinov_lab/brain_genomics/data_HCP/ts_missed'
'''
new_is_better = 0 
old_is_better = 0 
for subj in valid: 
    n0 = 0 
    n1 = 0 
    with h5py.File('{}/{}.mat'.format(mat0_dir, subj), 'r') as f0: 
        objs = f0.get('timeseries')
        for t in range(4):
            try:
                obj = objs[0][t]
                name = h5py.h5r.get_name(obj, f0.id)
                curr_ts = np.array(f0[name])
                if curr_ts.size > 2: n0 += 1 
            except: 
                pass  
    with h5py.File('{}/{}.mat'.format(mat1_dir, subj), 'r') as f1: 
        objs = f1.get('timeseries')
        for t in range(4):
            try:
                obj = objs[0][t]
                name = h5py.h5r.get_name(obj, f1.id)
                curr_ts = np.array(f1[name])
                if curr_ts.size > 2: n1 += 1 
            except: 
                pass  
    if n1 > n0: 
        new_is_better += 1 
    if n0 >= n1: 
        old_is_better += 1 

print('new is better: {}'.format(new_is_better))
print('old is better: {}'.format(old_is_better))
print('total: {}'.format(len(valid)))
'''

## move valid ts_missed files to timeseries dir 
for subj in valid:
    old_file = '{}/{}.mat'.format(mat0_dir, subj)
    new_file = '{}/{}.mat'.format(mat1_dir, subj)
    os.system('cp {} {}'.format(new_file, old_file)) 

