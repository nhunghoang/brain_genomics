'''
Script for gathering ensembl-symbol mappings 
for genes that pass the PrediXcan threshold, 
such that these mappings refer to PrediXVU filenames. 

- Nhung, Feb 2023 
'''

import numpy as np 
import h5py 
import sqlite3 as sql 
import os 
import subprocess as sp

## paths 
exp_path = '/data1/rubinov_lab/brain_genomics/scripts_assoc_clean/inputs_HCP/expr_regress'
pdx_path = '/data1/rubinov_lab/brain_genomics/models_PrediXcan_v8/models_by_tissue/brain'
pvu_path = '/dors/capra_lab/data/predixVU/byGene'

## out files 
map_path = '/data1/rubinov_lab/brain_genomics/models_PrediXcan_v8/predixvu_ens_sym_map.hdf5' 
dne_path = '/data1/rubinov_lab/brain_genomics/models_PrediXcan_v8/predixvu_ens_sym_dne.txt' 

## parse genes used in our assoc analyses 
assoc_ens = [] 
for rfile in os.listdir(exp_path): 
    with h5py.File(exp_path + '/' + rfile, 'r') as f: 
        rgenes = f['genes'][()].astype(str) 
    assoc_ens.extend(rgenes) 
assoc_ens = np.unique(assoc_ens)

## parse PrediXcan model dbs for ensembl-symbol mapping 
pdx_map = {} ## k: ensembl, v: symbol 
for pfile in os.listdir(pdx_path): 
    con = sql.connect(pdx_path + '/' + pfile)
    cur = con.cursor() 
    res = cur.execute('SELECT gene, genename FROM extra')
    qry = res.fetchall() 
    for (ens, sym) in qry: 
        pdx_map[ens.split('.')[0]] = sym 

## try to find ens-sym mappings for genes (relevant 
## to downstream analyses) that correspond to existing 
## PrediXVU files  
sym_map = {} ## k: ensembl, v: sym 
dne_map = {} ## k: ensembl, v: sym 

for i, ens in enumerate(assoc_ens): 
    if (i%100==0): print('{}/{}'.format(i, assoc_ens.size)) 

    sym = pdx_map[ens] 
    pvu_file = '{}/{}_{}_predixVU.csv.gz'.format(pvu_path, sym, ens)

    ## case 1: file exists, add to map dict   
    if os.path.exists(pvu_file): 
        sym_map[ens] = sym 
        continue 

    ## case 2: file DNE, do a manual search 
    try:
        cmd = 'ls {}/*{}*'.format(pvu_path, ens)
        sp_out = sp.check_output(cmd, shell=True, stderr=sp.DEVNULL)
    except:
        try:
            cmd = 'ls {}/*{}*'.format(pvu_path, sym)
            sp_out = sp.check_output(cmd, shell=True, stderr=sp.DEVNULL)
        except:
            ## case 3: no file found for this gene, move on 
            dne_map[ens] = sym
            print('[{}] DNE: {} {}'.format(len(dne_map), ens, sym))
            continue

    ## case 2 cont.: add to map dict
    fname = sp_out.decode().strip()
    ens1 = 'ENSG' + fname.split('/')[-1].split('ENSG')[1].split('_')[0]
    sym1 = fname.split('/')[-1].split('ENSG')[0][:-1]
    sym_map[ens1] = sym1

print('{}/{}'.format(i+1, assoc_ens.size)) 

## save mapping 
with h5py.File(map_path, 'w') as f: 
    for ens, sym in sym_map.items(): 
        k = ens.split('.')[0]
        f[k] = sym.encode()  

## save DNE mapping 
with open(dne_path, 'w') as f: 
    for ens, sym in dne_map.items(): 
        line = '{}\t{}\n'.format(ens, sym) 
        f.write(line) 
    


