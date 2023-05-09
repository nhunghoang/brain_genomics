'''
Convert genotype probabilities to 
dosages as described by MetaXcan. 

Note: Genotype files contain only the SNPs 
needed for the given GREx models. 

- Nhung, May 2023 
'''

from multiprocessing import Pool
import numpy as np 
import sys 

group = sys.argv[1]
model = sys.argv[2]

## paths 
top_path = '/data1/rubinov_lab/brain_genomics'
vcf_path = '{}/scripts_twas/inputs_{}/vcf_{}'.format(top_path, group, model)

## out path 
dos_path = '{}/scripts_twas/inputs_{}/dosage_{}'.format(top_path, group, model)

## func: read vcf file and convert to dosage 
def convert2dosage(data): 

    chrom = data['chrom']
    start = data['start']
    end = data['end']
    version = data['version']

    ## group-specific format differences 
    header = {'HCP': 0, 'UKB': 9}

    ## chr rsid pos a0 a1 MAF=0 dos ... 
    new_lines = [] 

    ## read vcf 
    vcf_file = '{}/c{}.vcf'.format(vcf_path, chrom)
    with open(vcf_file, 'r') as f: 
        for i, line in enumerate(f): 
            if i < header[group]: continue 

            ## versioning 
            if i < start: continue 
            if i == end: break

            ## var description 
            info = line.strip().split('\t')
            rsid = info[2].split(';')[0]
            posn = info[1]
            refa = info[3]
            alta = info[4] 
            maf0 = '9' 
            
            desc = '\t'.join([chrom, rsid, posn, refa, alta, maf0])
            
            ## dosages
            prob = info[9:]
            if group == 'HCP': 
                dosg = [np.array(p.split(','), dtype=float) \
                        .argmax() for p in prob]

            if group == 'UKB':
                dosg = [np.array(p.split(':')[1].split(','), dtype=float) \
                        .argmax() for p in prob]

            dosg = np.array(dosg, dtype=str)
            vals = '\t'.join(dosg)

            ## bring together 
            new_line = desc + '\t' + vals + '\n'
            new_lines.append(new_line)

            if (i%500) == 0: print(i)
    print('done converting chr {}.{}'.format(chrom, version))

    ## write to file 
    dos_file = '{}/c{}_v{}.dosage.txt'.format(dos_path, chrom, version)
    with open(dos_file, 'w') as f: 
        f.writelines(new_lines)

chr_ver = [] 
for i, s in enumerate(np.arange(0, 20000, 1000)):
    for c in np.arange(22, 0, -1):
        chr_ver.append({'chrom':str(c), \
                        'start':s, \
                        'end':s+1000, \
                        'version':str(i) \
                       })

pool = Pool(processes=15)
pool.map(convert2dosage, chr_ver)

