#!/bin/bash 
## Script for running PrediXcan Association 
## for regional or inter-regional associations. 

## Nhung, April 2023 

module purge
ml GCC/8.2.0 OpenMPI/3.1.4 Python/3.7.2
source /data1/rubinov_lab/brain_genomics/models_PrediXcan_v8/venv_predixcan/bin/activate

######################
GROUP=UKB/replications/50_0 #UKB/nonBrit
MODEL=JTI
WHICH=same ## same or cross 
PHENS=vol_mean
######################

top_dir="/data1/rubinov_lab/brain_genomics"
twas_script="${top_dir}/models_PrediXcan_v8/MetaXcan/software/PrediXcanAssociation.py"

## inputs 
grex_path="${top_dir}/paper_twas/inputs_${GROUP}/grex_${MODEL}"
phen_file="${top_dir}/paper_twas/inputs_${GROUP}/phenotypes/${PHENS}.csv"
covs_file="${top_dir}/paper_twas/inputs_${GROUP}/covariates.csv"

## output 
outs_dir="${top_dir}/paper_twas/outputs_${GROUP}/twas_${MODEL}"

## brain regions
brain_names=(hippocampus amygdala caudate nucleus-accumbens putamen \
             dlpfc cerebellar-hemisphere anterior-cingulate)

## covariates 
covs=(age isMale) 
for i in {1..40}
do
    covs+=(PC${i})
done

## TWAS
if [[ $WHICH = same ]]
then

    ## (same reg for grex and phen) 
    for idx in ${!brain_names[*]}
    do 
        echo -e "- ${brain_names[$idx]} -"
        python -u ${twas_script} \
            --hdf5_expression_file ${grex_path}/${brain_names[$idx]}.hdf5 \
            --input_phenos_file ${phen_file} \
            --input_phenos_column ${brain_names[$idx]} \
            --covariates_file ${covs_file} \
            --covariates ${covs[*]} \
            --output ${outs_dir}/${PHENS}/${brain_names[$idx]}.txt & 
    done

elif [[ $WHICH = cross ]]
then
    
    ## (cross reg for grex and phen)
    for idx1 in ${!brain_names[*]}
    do
        for idx2 in ${!brain_names[*]}
        do
            if [[ $idx1 -eq $idx2 ]]; then
                continue
            fi

            b1=${brain_names[$idx1]}
            b2=${brain_names[$idx2]}

            if [ ${b1} != "dlpfc" ] && [ ${b1} != "dlpfc_psychencode" ]; then
                continue 
            fi

            echo -e "- ${b1} expr, ${b2} phen -"
            python -u ${twas_script} \
                --hdf5_expression_file ${grex_path}/${b1}.hdf5 \
                --input_phenos_file ${phen_file} \
                --input_phenos_column ${b2} \
                --covariates_file ${covs_file} \
                --covariates ${covs[*]} \
                --output ${outs_dir}/cross_regs/${PHENS}/grex_${b1}_phen_${b2}.txt &
        done
    done

fi

