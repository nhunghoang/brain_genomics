#!/bin/bash

## Script for running PrediXcan v8 models (only the 10 brain 
## tissues of interest) that pass predictionn quality thresholds. 
##
## Usage Notes: 
## ACCRE modules: GCC/8.2.0 OpenMPI/3.1.4 Python/3.7.2
## python venv: /data1/rubinov_lab/brain_genomics/models_PrediXcan_v8/venv_predixcan
## venv packages: [...]/models_PrediXcan_v8/packages_predixcan.txt
##
## Variables to modify: DATASET, PRED_R2, PRED_PV, 0/1-based SNP position

## Nhung, Tim, Jee Hyun; updated Jan 2022  

## genotype data
DATASET=HCP ## or UKB 
dosages="/data1/rubinob_lab/brain_genomics/data_${DATASET}/predixcan_v8/dosage_format/zipped/chr*.dosage.txt.gz" 
samples="/data1/rubinov_lab/brain_genomics/data_${DATASET}/predixcan_v8/sample_ids.txt"
outputs="/data1/rubinov_lab/brain_genomics/models_PrediXcan_v8/expression_${DATASET}"

## PrediXcan models 
new_names=(hippocampus frontal-pole amygdala anterior-cingulate caudate cerebellar-hemisphere hypothalamus nucleus-accumbens putamen substantia-nigra)
pdx_names=(Hippocampus Frontal_Cortex_BA9 Amygdala Anterior_cingulate_cortex_BA24 Caudate_basal_ganglia Cerebellar_Hemisphere Hypothalamus Nucleus_accumbens_basal_ganglia Putamen_basal_ganglia Substantia_nigra)

pdx_script="/data1/rubinov_lab/brain_genomics/models_PrediXcan_v8/MetaXcan-master/software/Predict.py"
pdx_models="/data1/rubinov_lab/brain_genomics/models_PrediXcan_v8/models_by_tissue/brain"

PRED_R2="0.30"
PRED_PV="0.01"

for idx in ${!new_names[*]}
do 

    ## query gene models that pass prediction quality check   
    prompt="SELECT gene FROM extra \
            WHERE \"pred.perf.R2\" >= ${PRED_R2} \
            AND \"pred.perf.pval\" <= ${PRED_PV};"
    query=$(sqlite3 ${pdx_models}/en_Brain_${pdx_names[$idx]}.db "${prompt}")

    ## call PrediXcan for only these models 
    echo -e "- ${new_names[$idx]} -"
    python -u ${pdx_script} \
    --model_db_path ${pdx_models}/en_Brain_${pdx_names[$idx]}.db \
    --only_entries ${query} \
    --text_genotypes ${dosages} \
    --text_sample_ids ${samples} \
    --prediction_output ${outputs}/${new_names[$idx]}.hdf5 HDF5 \
    --prediction_summary_output ${outputs}/${new_names[$idx]}.log \
    --zero_based_positions ## this argument needs to be checked for each dataset  
done

