#!/bin/bash 

## Script for running PrediXcan Association 
## on permutations of HCP that match the 
## cohort size of non-twins 

## Nhung, July 2024 

module purge
ml GCC/8.2.0 OpenMPI/3.1.4 Python/3.7.2
source /data1/rubinov_lab/brain_genomics/models_PrediXcan_v8/venv_predixcan/bin/activate

top_dir="/data1/rubinov_lab/brain_genomics"
twas_script="${top_dir}/models_PrediXcan_v8/MetaXcan/software/PrediXcanAssociation.py"

## inputs 
grex_path="${top_dir}/paper_twas/inputs_HCP/nonTwin/grex_JTI"
phen_path="${top_dir}/paper_twas/inputs_HCP/nonTwin/phenotypes"
covs_file="${top_dir}/paper_twas/inputs_HCP/nonTwin/covariates.csv"

## output 
outs_dir="${top_dir}/paper_twas/outputs_HCP/nonTwin/twas_JTI/perms"

## brain regions
brain_names=(hippocampus amygdala caudate nucleus-accumbens putamen \
             cerebellar-hemisphere anterior-cingulate dlpfc)

## phenotypes 
phens=(vol_mean alff_mean reho_noGS_mean connmean_noGS_mean)

## covariates 
covs=(age isMale) 
for i in {1..40}
do
    covs+=(PC${i})
done

## TWAS
acc=1
for itr in {2000..2999}
do
    for pidx in ${!phens[*]}
    do
        for ridx in ${!brain_names[*]}
        do 
            echo -e "[${itr}] ${phens[$pidx]} ${brain_names[$ridx]}"
            gfile="${grex_path}/${brain_names[$ridx]}.hdf5"
            pfile="${phen_path}/perm_${phens[$pidx]}.csv"
            pcols="${brain_names[$ridx]}_${itr}"
            ofile="${outs_dir}/${phens[$pidx]}/${brain_names[$ridx]}_${itr}.txt"

            if (( acc % 35 == 0)); then
                nice -n 18 python -u ${twas_script} \
                    --hdf5_expression_file ${gfile} \
                    --input_phenos_file ${pfile} \
                    --input_phenos_column ${pcols} \
                    --covariates_file ${covs_file} \
                    --covariates ${covs[*]} \
                    --output ${ofile}
            else 
                nice -n 18 python -u ${twas_script} \
                    --hdf5_expression_file ${gfile} \
                    --input_phenos_file ${pfile} \
                    --input_phenos_column ${pcols} \
                    --covariates_file ${covs_file} \
                    --covariates ${covs[*]} \
                    --verbosity 100 \
                    --output ${ofile} &
            fi
            (( acc++ ))
        done
    done 
done
