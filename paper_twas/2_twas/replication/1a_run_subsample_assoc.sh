#!/bin/bash 

## Script for running PrediXcan Association 
## on a random subsample of HCP that matches 
## the cohort size of non-twins

## Specific to HCP vol mean 

## Nhung, Feb 2024 

module purge
ml GCC/8.2.0 OpenMPI/3.1.4 Python/3.7.2
source /data1/rubinov_lab/brain_genomics/models_PrediXcan_v8/venv_predixcan/bin/activate

top_dir="/data1/rubinov_lab/brain_genomics"
twas_script="${top_dir}/models_PrediXcan_v8/MetaXcan/software/nh_PrediXcanAssociation.py"

## inputs 
grex_path="${top_dir}/paper_twas/inputs_HCP/allEuro/grex_JTI"
phen_file="${top_dir}/paper_twas/inputs_HCP/allEuro/phenotypes/vol_mean.csv"
covs_file="${top_dir}/paper_twas/inputs_HCP/allEuro/covariates.csv"

## output 
outs_dir="${top_dir}/paper_twas/outputs_HCP/allEuro/subsampling"

## brain regions
brain_names=(hippocampus amygdala caudate nucleus-accumbens putamen \
             cerebellar-hemisphere anterior-cingulate dlpfc)

## covariates 
covs=(age isMale) 
for i in {1..40}
do
    covs+=(PC${i})
done

## TWAS
for itr in {0..99}
do
    for idx in ${!brain_names[*]}
    do 
        echo -e "[${itr}] ${brain_names[$idx]}"
        if [ $((itr % 3)) -eq 0 ] && [ "${brain_names[$idx]}" = 'dlpfc' ]; then
            nice -n 18 python -u ${twas_script} \
                --hdf5_expression_file ${grex_path}/${brain_names[$idx]}.hdf5 \
                --input_phenos_file ${phen_file} \
                --input_phenos_column ${brain_names[$idx]} \
                --covariates_file ${covs_file} \
                --covariates ${covs[*]} \
                --sample_iter HCP_${itr} \
                --output ${outs_dir}/nonTwin_${itr}_${brain_names[$idx]}.txt 
                #--verbosity 100 \
        else 
            nice -n 18 python -u ${twas_script} \
                --hdf5_expression_file ${grex_path}/${brain_names[$idx]}.hdf5 \
                --input_phenos_file ${phen_file} \
                --input_phenos_column ${brain_names[$idx]} \
                --covariates_file ${covs_file} \
                --covariates ${covs[*]} \
                --sample_iter HCP_${itr} \
                --verbosity 100 \
                --output ${outs_dir}/nonTwin_${itr}_${brain_names[$idx]}.txt & 
        fi
    done 
done
