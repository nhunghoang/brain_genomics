#!/bin/bash

## Script for running PrediXcan v8 models (all 13 available brain tissues)   
## NOTE: this script needs to be ran with modules (GCC/8.2.0 OpenMPI/3.1.4 Python/3.7.2) and 
##       with venv '/data1/rubinov_lab/brain_genomics/models_PrediXcan_v8/venv_predixcan' 
##       which contains packages listed in '[...]/models_PrediXcan_v8/packages_predixcan.txt' 

## Tim, Jee Hyun, Nhung; modified May 2021  

prefix=(hippocampus frontal-pole amygdala anterior-cingulate caudate cerebellar-hemisphere cerebellum cortex hypothalamus nucleus-accumbens putamen spinal-cord substantia-nigra)
models=(Hippocampus Frontal_Cortex_BA9 Amygdala Anterior_cingulate_cortex_BA24 Caudate_basal_ganglia Cerebellar_Hemisphere Cerebellum Cortex Hypothalamus Nucleus_accumbens_basal_ganglia Putamen_basal_ganglia Spinal_cord_cervical_c1 Substantia_nigra)

predixcan8='/data1/rubinov_lab/brain_genomics/models_PrediXcan_v8/MetaXcan-master/software/Predict.py'
models_dir='/data1/rubinov_lab/brain_genomics/models_PrediXcan_v8/models_by_tissue/brain'
dosg_files='/data1/rubinov_lab/brain_genomics/data_HCP/predixcan_v8/dosage_format/zipped/chr*.dosage.txt.gz'
sample_ids='/data1/rubinov_lab/brain_genomics/data_HCP/predixcan_v8/sample_ids.txt' 
output_dir='/data1/rubinov_lab/brain_genomics/models_PrediXcan_v8/expression_HCP'

for idx in ${!prefix[*]}
do 
    echo -e "- ${prefix[$idx]} -"
    python -u ${predixcan8} \
    --model_db_path ${models_dir}/en_Brain_${models[$idx]}.db \
    --text_genotypes ${dosg_files} \
    --text_sample_ids ${sample_ids} \
    --prediction_output ${output_dir}/${prefix[$idx]}.hdf5 HDF5 \
    --prediction_summary_output ${output_dir}/${prefix[$idx]}.log \
    --zero_based_positions 
done

