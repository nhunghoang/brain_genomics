'''
Gather HCP subject demographics. 

- Nhung, updated Feb 2024
'''

import pandas as pd 

main_path = '/data1/rubinov_lab/brain_genomics'

idds_path = f'{main_path}/data_HCP/HCP_dbGaP_all/GenotypeFiles/origgeno/HCPyoungadult_ImputationSubjectSampleMapping_MEGA.xlsx' 
dem1_path = f'{main_path}/data_HCP/subject_demographics/COMPLETE_DEMOGRAPHICS.csv'
dem2_path = f'{main_path}/data_HCP/subject_demographics/HCP_unrestricted_4_23_2019_21_59_10.csv'

outs_path = f'{main_path}/paper_twas/inputs_HCP/demographics.csv'

## sample id, subject id, age, twin GT, family id, mom id, dad id, race, ethnicity, gender 

## subject/sample mapping 
id_map = pd.read_excel(idds_path, index_col='SUBJECT_ID').to_dict()['SAMPLE_ID']

## load demographics 
c0 = ['Subject', 'Age_in_Yrs', 'ZygosityGT', 'Family_ID', 'Mother_ID', 'Father_ID', 'Race', 'Ethnicity'] 
d0 = pd.read_csv(dem1_path, usecols=c0)

c1 = ['Subject', 'Gender']
d1 = pd.read_csv(dem2_path, usecols=c1)

## clean up table 
df = d0.merge(d1, on='Subject', how='inner')

cdict = {'Subject': 'subject', 'Age_in_Yrs': 'age', 'Family_ID': 'family_id', \
         'Mother_ID': 'mother_id', 'Father_ID': 'father_id', 'Race': 'race', \
         'Ethnicity': 'ethnicity', 'Gender': 'gender', 'ZygosityGT': 'zygosity_GT'}
df = df.rename(columns=cdict)

## add info 
df['sample'] = df['subject'].map(id_map)
df['twin_id'] = df.apply(lambda x: x['family_id'] if x['zygosity_GT'] == 'MZ' else '', axis=1)

cols = ['subject', 'sample', 'age', 'zygosity_GT', 'twin_id', 'family_id', 'mother_id', 'father_id', 'race', 'ethnicity', 'gender']
df = df[cols] 

## save 
df.to_csv(outs_path, index=False)

