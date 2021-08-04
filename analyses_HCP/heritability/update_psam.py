'''

'''

import sys

psamf = sys.argv[1]
FAM_RACE_ETH = '/data1/rubinov_lab/brain_genomics/data_HCP/subject_demographics/sampleID_race_familyID_ethnicity.txt'

# get family info
with open(FAM_RACE_ETH, 'r') as f:
    info = [i.strip().split('\t') for i in f.readlines()[1:]]
family = {ii[0]:ii[2] for ii in info} # k: sample id, v: family id
# remove twin indicator
for sampl in family.keys():
    if family[sampl][-1]=='T': family[sampl] = family[sampl][:-1]

with open(psamf, 'r') as f: lines = f.readlines()[1:]

with open(psamf, 'w') as f:
    f.write('#FID\tIID\tSEX\n')
    for line in lines:
        iid = line.split('\t')[0]
        fid = family[iid]
        new_line = fid + '\t' + line
        f.write(new_line)
