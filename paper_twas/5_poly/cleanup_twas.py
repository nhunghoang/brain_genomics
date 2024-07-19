'''
Clean up TWAS output files. 

- Nhung, Feb 2024
'''

import sys 
import os 
import subprocess

phen = sys.argv[1]
#main_path = '/data1/rubinov_lab/brain_genomics/paper_twas/outputs'
#coho_path = f'{main_path}_{cohort}'
#coho_path = f'/data/rubinov_lab/brain_genomics_project/polygenic/twas/{phen}'
coho_path = f'/data1/rubinov_lab/brain_genomics/paper_twas/outputs_HCP/nonTwin/twas_JTI/perms/{phen}'

tfiles = os.listdir(coho_path)
for t, tpath in enumerate(tfiles): 
    tfile = f'{coho_path}/{tpath}'

    substring = "b"
    cmd = f"sed -n '2p' {tfile} | grep -q {substring}"
    res = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if res.returncode != 0:
        if (t%100) == 0: 
            print('{} / {}'.format(t, len(tfiles)))
        continue

    with open(tfile, 'r') as f: 
        lines = f.readlines()

    new_lines = [lines[0]]
    for line in lines[1:]: 

        if '"' in line: 
            continue 

        if 'pheno' in line: 
            continue 

        if line[0] == 'b': 
            new_line = line.replace('b', '').replace("'", "")
            new_lines.append(new_line)

    with open(tfile, 'w') as f: 
        f.writelines(new_lines)

    if (t%100) == 0: 
        print('{} / {}'.format(t, len(tfiles)))

        
    

