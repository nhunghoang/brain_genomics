#!/bin/bash

## Use JTI to impute GREx

## params (temp) 
reg_model="JTI_Brain_Putamen_basal_ganglia"
phen_name="putamen"

## path
top_path="/data1/rubinov_lab/brain_genomics/proj_nova" 
run_path="${top_path}/software/MetaXcan/software/Predict.py"

mod_path="${top_path}/models/JTI/weights/${reg_model}.db"
gen_path="${top_path}/data/dosages"
sam_path="${top_path}/data/cohort.txt"

out_path="${top_path}/outputs/grex/jti_${phen_name}"

## grex
python -u ${run_path} \
        --model_db_path ${mod_path} \
        --text_genotypes ${gen_path}/c*.dosage.txt \
        --text_sample_ids ${sam_path} \
        --prediction_output ${out_path}.hdf5 HDF5 \
        --prediction_summary_output ${out_path}.log 

