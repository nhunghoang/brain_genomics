#!/bin/bash

## Script for running any imputation models (only the brain 
## tissues of interest) that pass prediction quality thresholds. 
##
## Usage Notes: 
## ACCRE modules: GCC/8.2.0 OpenMPI/3.1.4 Python/3.7.2
## python venv: /data1/rubinov_lab/brain_genomics/models_PrediXcan_v8/venv_predixcan
## venv packages: [...]/models_PrediXcan_v8/packages_predixcan.txt
##
## Variables to modify: GROUP, MODEL, PRED_R2, PRED_PV, 0/1-based SNP position

## updated by Nhung (April 2023)  

module purge 
ml GCC/8.2.0 OpenMPI/3.1.4 Python/3.7.2
source /data1/rubinov_lab/brain_genomics/models_PrediXcan_v8/venv_predixcan/bin/activate

GROUP=HCP
MODEL=JTI  

#PRED_R2="0.30"
#PRED_PV="0.01"

top_dir="/data1/rubinov_lab/brain_genomics"

## input genotype dosages 
genotyp="${top_dir}/paper_twas/inputs_${GROUP}/dosage_${MODEL}/c*.dosage.txt"

## samples 
samples="${top_dir}/paper_twas/inputs_${GROUP}/cohort.txt"

## GREx scripts
grex_script="${top_dir}/models_PrediXcan_v8/MetaXcan/software/Predict.py"
models_path="${top_dir}/models_${MODEL}/models_by_tissue"

## output GREx path
outgrex="${top_dir}/paper_twas/inputs_${GROUP}/grex_${MODEL}"

## brain model names 
short_names=(hippocampus amygdala caudate nucleus-accumbens putamen \
            cerebellar-hemisphere anterior-cingulate dlpfc)
model_names=(Hippocampus Amygdala Caudate_basal_ganglia \
            Nucleus_accumbens_basal_ganglia Putamen_basal_ganglia \
            Cerebellar_Hemisphere Anterior_cingulate_cortex_BA24 Frontal_Cortex_BA9)

short_names=(dlpfc_psychencode)
model_names=(DLPFC_PsychEncode)

## GREx inference
for idx in ${!short_names[*]}
do 
    ## query gene models that pass prediction quality check   
    #prompt="SELECT gene FROM extra \
    #        WHERE \"pred.perf.R2\" >= ${PRED_R2} \
    #        AND \"pred.perf.pval\" <= ${PRED_PV};"
    #query=$(sqlite3 ${models_path}/${MODEL}_Brain_${model_names[$idx]}.db "${prompt}")

    ## call PrediXcan/JTI/etc for only these models 
    echo -e "- ${short_names[$idx]} -"
    python -u ${grex_script} \
    --model_db_path ${models_path}/${MODEL}_Brain_${model_names[$idx]}.db \
    --text_genotypes ${genotyp} \
    --text_sample_ids ${samples} \
    --prediction_output ${outgrex}/${short_names[$idx]}.hdf5 HDF5 \
    --prediction_summary_output ${outgrex}/logs/${short_names[$idx]}.log & #\
    #--only_entries ${query} &
    #--zero_based_positions ## this argument needs to be checked for each dataset  
done

