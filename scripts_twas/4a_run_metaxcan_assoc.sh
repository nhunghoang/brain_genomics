#!/bin/bash 
## Script for running PrediXcan Association 

## Nhung, April 2023 

module purge
ml GCC/8.2.0 OpenMPI/3.1.4 Python/3.7.2
source /data1/rubinov_lab/brain_genomics/models_PrediXcan_v8/venv_predixcan/bin/activate

GROUP=UKB
MODEL=JTI
top_dir="/data1/rubinov_lab/brain_genomics"

## UKB data 
grex_path="${top_dir}/scripts_twas/outputs_${GROUP}/grex_${MODEL}"
phen_file="${top_dir}/scripts_twas/inputs_${GROUP}/volume_phenotypes.csv"
covs_file="${top_dir}/scripts_twas/inputs_${GROUP}/covariates.csv"

## output 
outs_dir="${top_dir}/scripts_twas/outputs_${GROUP}/twas_metaxcan_${MODEL}"

## brain regions
brain_names=(hippocampus amygdala caudate nucleus-accumbens putamen \
             cerebellar-hemisphere) # frontal-pole anterior-cingulate)

## covariates 
covs=(eid age isMale) 
for i in {1..40}
do
    covs+=(PC${i})
done

twas_script="${top_dir}/models_PrediXcan_v8/MetaXcan/software/PrediXcanAssociation.py"

for idx in ${!brain_names[*]}
do 
    echo -e "- ${brain_names[$idx]} -"

    python -u ${twas_script} \
        --hdf5_expression_file ${grex_path}/${brain_names[$idx]}.hdf5 \
        --input_phenos_file ${phen_file} \
        --input_phenos_column left_vol_${brain_names[$idx]} \
        --covariates_file ${covs_file} \
        --covariates ${covs[*]} \
        --output ${outs_dir}/left_vol_${brain_names[$idx]}.txt &

    python -u ${twas_script} \
        --hdf5_expression_file ${grex_path}/${brain_names[$idx]}.hdf5 \
        --input_phenos_file ${phen_file} \
        --input_phenos_column right_vol_${brain_names[$idx]} \
        --covariates_file ${covs_file} \
        --covariates ${covs[*]} \
        --output ${outs_dir}/right_vol_${brain_names[$idx]}.txt &
done
         
