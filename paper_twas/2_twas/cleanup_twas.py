'''
Clean up TWAS output files. 

- Nhung, Feb 2024
'''

import sys 
import os 

cohort = sys.argv[1]
main_path = '/data1/rubinov_lab/brain_genomics/paper_twas/outputs'
coho_path = f'{main_path}_{cohort}'

tfiles = os.listdir(coho_path)
for t, tpath in enumerate(tfiles): 
    tfile = f'{coho_path}/{tpath}'
    with open(tfile, 'r') as f: 
        lines = f.readlines()

    if lines[1][0] != 'b': continue 

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

        
    

