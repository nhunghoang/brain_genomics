'''
From the PrediXcan v8 brain models, gather all unique variants and sort into chrom-specific vcf files. 
Also write tissue-specific lists of all gene-snp pairs.   

-Nhung, updated May 2021 
'''

import os 
import numpy as np
import sqlite3

snps_per_gene_dir = '/data1/rubinov_lab/brain_genomics/models_PrediXcan_v8/snps_per_gene'
snps_by_chr_dir = '/data1/rubinov_lab/brain_genomics/models_PrediXcan_v8/snps_by_chr'
models_dir = '/data1/rubinov_lab/brain_genomics/models_PrediXcan_v8/models_by_tissue/brain'

chr_set = {'chr'+str(c):[] for c in range(1,23)} # k: chr, v: [[pos,ref,alt,rsid]]

## loop through model files 
for model_file in os.listdir(models_dir): 
    if model_file[-3:] != '.db': continue 
    tissue = model_file[9:-3] ## elastic net naming convention 
    print('- query: {}'.format(tissue.upper()))

    ## sqlite3 query 
    conn = sqlite3.connect(models_dir + '/' + model_file)
    c = conn.cursor()
    c.execute('SELECT varID,rsid,gene FROM weights')
    query = c.fetchall()

    ## open model-specific snps_per_gene file 
    f0 = open('{}/{}.txt'.format(snps_per_gene_dir, tissue), 'w') 

    for q in query: 
        ## save snp info to chrom dict 
        [chrom, pos, ref, alt, build] = q[0].split('_')
        rsid = q[1]  
        chr_set[chrom].append([pos,ref,alt,rsid])
        ## write gene-snp info to file  
        gene = q[2]
        line = '{}\t{}\n'.format(gene, rsid)
        f0.write(line) 

    f0.close()

## loop through chromosomes 
print('unique snps out of total across all tissues:')
uniq_snps = {} # k: chr, v: [[pos,ref,alt,rsid]]
for chrom in chr_set.keys():
    total = len(chr_set[chrom])
    uniqs = np.unique(chr_set[chrom], axis=0)
    uniq_snps[chrom] = uniqs
    print('- save: {}, {} / {} '.format(chrom, uniqs.shape[0], total))     

    vfile = '{}/{}.vcf'.format(snps_by_chr_dir, chrom) 
    with open(vfile, 'w') as f1: 
        f1.write('##fileformat=VCFv4.2\n')
        f1.write('\t'.join(('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', '\n')))
        for [pos,ref,alt,rsid] in uniq_snps[chrom]:
            f1.write('\t'.join((chrom, pos, rsid, ref, alt, '.', '.', '.', '\n')))
    

