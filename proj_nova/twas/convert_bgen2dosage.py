'''
Convert bgen files to dosage files. 
'''

import numpy as np 
import pandas as pd 
import bgen_reader as bgr

from multiprocessing import Pool
from time import time
import subprocess

top_path = '/data1/rubinov_lab/brain_genomics/proj_nova'
gen_path = f'{top_path}/data/genotypes'

out_path = f'{top_path}/data/dosages'

def compute_dosages(params): 

    ver = params['ver']
    chrom = str(params['chr'])
    snps = params['snps']
    gens = params['gens']

    dos_lines = [] 
    for snp, gen in zip(snps.itertuples(), gens):

        rsid = snp.rsid
        post = str(snp.pos)
        [ref, alt] = snp.allele_ids.split(',')

        dose = np.argmax(gen, axis=1) ## (samples,)

        maf = '9' ## not used
        rsid_info = '\t'.join([chrom, rsid, post, ref, alt, maf])
        dose_vals = '\t'.join(dose.astype(str))
        dose_line = f'{rsid_info}\t{dose_vals}\n'
        dos_lines.append(dose_line)

        #if (i%1000) == 0: print('chr {}.{}: {}'.format(chrom, ver, i))

    dfile = f'{out_path}/tmp/c{chrom}_v{ver}.dosage.txt'  
    with open(dfile, 'w') as f: 
        f.writelines(dos_lines)
    print('done converting chr {}.{}'.format(chrom, ver))

## main 
max_variant = 5000
num_threads = 10
pool = Pool(processes=num_threads)

for chrom in range(7, 0, -1): 
    bgen = bgr.read_bgen(f'{gen_path}/c{chrom}.bgen', verbose=False)
    snps = bgen['variants']
    gens = bgen['genotype']

    nvar = len(snps)
    nblocks = nvar // max_variant
    snps = snps.repartition(npartitions=nblocks)
    snp_blocks = snps.to_delayed()

    num_versions = np.ceil(nblocks / num_threads).astype(int)
    print(f'\nChr {chrom} contains {nvar} variants')
    print(f'-- {nblocks} blocks across {num_versions} multi-threaded runs expected')

    if nblocks != len(snp_blocks): 
        lb = len(snp_blocks)
        print(f'-- NOTE: {nblocks} blocks expected, {lb} blocks created')

    params = [] 
    start = 0
    end = None
    t0 = time()
    for i, block in enumerate(snp_blocks): 
        snp_set = block.compute()
        num_snp = snp_set.shape[0]
        end = start + num_snp

        gen_set = [g.compute()['probs'] for g in gens[start:end]]
        start = end
        print(f'processing chr {chrom} block {i}...')
        
        pp = {'ver': i, 'chr': chrom, 'snps': snp_set, 'gens': gen_set}
        params.append(pp)

        if (len(params) == num_threads) or (i == (nblocks-1)):  
            pool.map(compute_dosages, params)
            params = [] 

            tpass = time() - t0
            hr = int(tpass // 3600)
            mn = int((tpass % 3600) // 60)
            sc = int((tpass % 3600) % 60)

            print(f'> {hr} hr {mn} mn {sc} sc < for this run')
            t0 = time()

    ## concat versions 
    vpath = f'{out_path}/tmp/c{chrom}_v*.dosage.txt'
    cpath = f'{out_path}/c{chrom}.dosage.txt' 
    cmd = f'cat {vpath} >> {cpath}'

    subprocess.run(cmd, shell=True)
    print(chrom)
         
