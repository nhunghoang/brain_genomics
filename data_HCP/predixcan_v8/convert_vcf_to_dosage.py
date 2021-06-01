'''
'''

vcf_dir = '/data1/rubinov_lab/brain_genomics/data_HCP/predixcan_v8/vcf_format'
dos_dir = '/data1/rubinov_lab/brain_genomics/data_HCP/predixcan_v8/dosage_format'

## save sample IDs (this is the order PrediXcan will run through) 
header_file = '/data1/rubinov_lab/brain_genomics/data_HCP/predixcan_v8/header_v8.vcf'
sample_file = '/data1/rubinov_lab/brain_genomics/data_HCP/predixcan_v8/sample_ids.txt'
sample_line = ''
with open(header_file, 'r') as hf: 
    line = hf.readlines()[-1] 
    for sid in line.strip().split('\t')[9:]: 
        sample_line += '{}\t{}\n'.format(sid,sid)
with open(sample_file, 'w') as sf: 
    sf.write(sample_line)

## function for computing dosage from genotype probability (GP)
## REF/REF REF/ALT ALT/ALT --> 0 1 2 
def compute_dosage(GP_list): 
    dosages = []
    for indiv in GP_list:
        indiv = indiv.split(',')
        probs = [float(j) for j in indiv]
        idx = probs.index(max(probs))
        dosages.append(str(idx))
    line = '\t'.join(dosages) 
    return line

## convert vcf files to dosage files 
## dosage format: chr rsid pos ref alt MAF dosage-in-sampleID-order 
## dosage: num of alt alleles (between 0 and 2)
MAF = '0' ## PrediXcan v8 does not use MAF but still looks for a number 
for chrom in range(22,0,-1): 

    with open('{}/chr{}.vcf'.format(vcf_dir, chrom), 'r') as vcf_file: 
        lines = vcf_file.readlines()

    dos_file = open('{}/chr{}.dosage.txt'.format(dos_dir, chrom), 'w') 
    for line in lines[5:]: 
        info = line.strip().split('\t') 
        [pos, rsid, ref, alt] = info[1:5] 
        dosages = compute_dosage(info[9:])
        new_info = [str(chrom), rsid, pos, ref, alt, MAF, dosages] # NOTE: PrediXcan expects rsid, then pos 
        new_line = '\t'.join(new_info) 
        dos_file.write(new_line + '\n')
    dos_file.close()  
    print('{}/22'.format(23-chrom))
