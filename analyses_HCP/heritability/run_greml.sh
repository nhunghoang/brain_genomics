#!/usr/bin/env bash

#REGS=(amygdala anterior-cingulate caudate cerebellar-hemisphere frontal-pole hippocampus nucleus-accumbens putamen) 
#REGS=(amygdala caudate cerebellar-hemisphere hippocampus hypothalamus nucleus-accumbens putamen substantia-nigra) 
PHNS=(connectivity_mean connectivity_variance timeseries_variance alff falff)
REGS=(amygdala anterior-cingulate caudate cerebellar-hemisphere frontal-pole hippocampus hypothalamus nucleus-accumbens putamen substantia-nigra) 
#PHNS=(connectivity_mean connectivity_variance participation_coefficient regional_homogeneity alff falff)

############################################################################################

## SM ## 
#PHN_DIR='/data1/rubinov_lab/brain_genomics/analyses_HCP/heritability/phen_files_SM'
#
#PLK_PDX='/data1/rubinov_lab/brain_genomics/analyses_HCP/heritability/plink_files_full/pdx'
#RES_PDX='/data1/rubinov_lab/brain_genomics/analyses_HCP/heritability/results_gcta_SM/pdx'
#
#PLK_ALL='/data1/rubinov_lab/brain_genomics/analyses_HCP/heritability/plink_files_full/all'
#RES_ALL='/data1/rubinov_lab/brain_genomics/analyses_HCP/heritability/results_gcta_SM/all'

## fSM ## 
PHN_DIR='/data1/rubinov_lab/brain_genomics/analyses_HCP/heritability/phen_files_fSM'

PLK_PDX='/data1/rubinov_lab/brain_genomics/analyses_HCP/heritability/plink_files_full/pdx'
RES_PDX='/data1/rubinov_lab/brain_genomics/analyses_HCP/heritability/results_gcta_fSM/pdx'

PLK_ALL='/data1/rubinov_lab/brain_genomics/analyses_HCP/heritability/plink_files_full/all'
RES_ALL='/data1/rubinov_lab/brain_genomics/analyses_HCP/heritability/results_gcta_fSM/all'

## HOA ## 
#PHN_DIR='/data1/rubinov_lab/brain_genomics/analyses_HCP/heritability/phen_files_HOA'
#
#PLK_PDX='/data1/rubinov_lab/brain_genomics/analyses_HCP/heritability/plink_files_full/pdx'
#RES_PDX='/data1/rubinov_lab/brain_genomics/analyses_HCP/heritability/results_gcta_HOA/pdx'
#
#PLK_ALL='/data1/rubinov_lab/brain_genomics/analyses_HCP/heritability/plink_files_full/all'
#RES_ALL='/data1/rubinov_lab/brain_genomics/analyses_HCP/heritability/results_gcta_HOA/all'

## fHOA ## 
#PHN_DIR='/data1/rubinov_lab/brain_genomics/analyses_HCP/heritability/phen_files_fHOA'
#
#PLK_PDX='/data1/rubinov_lab/brain_genomics/analyses_HCP/heritability_v7/subset_plink_files'
#RES_PDX='/data1/rubinov_lab/brain_genomics/analyses_HCP/heritability/results_gcta_fHOA/pdx' 
#
#PLK_ALL='/data1/rubinov_lab/brain_genomics/analyses_HCP/heritability_v7/plink_files'
#RES_ALL='/data1/rubinov_lab/brain_genomics/analyses_HCP/heritability/results_gcta_fHOA/all' 

## fHOA2 ## 
#PHN_DIR='/data1/rubinov_lab/brain_genomics/analyses_HCP/heritability/phen_files_fHOA2'
#
#PLK_PDX='/data1/rubinov_lab/brain_genomics/analyses_HCP/heritability/plink_files_full/pdx'
#RES_PDX='/data1/rubinov_lab/brain_genomics/analyses_HCP/heritability/results_gcta_fHOA2/pdx'
#
#PLK_ALL='/data1/rubinov_lab/brain_genomics/analyses_HCP/heritability/plink_files_full/all'
#RES_ALL='/data1/rubinov_lab/brain_genomics/analyses_HCP/heritability/results_gcta_fHOA2/all'

## fMPM ##
#PHN_DIR='/data1/rubinov_lab/brain_genomics/analyses_HCP/heritability/phen_files_fMPM'
#
#PLK_PDX='/data1/rubinov_lab/brain_genomics/analyses_HCP/heritability/plink_files_full/pdx'
#RES_PDX='/data1/rubinov_lab/brain_genomics/analyses_HCP/heritability/results_gcta_fMPM/pdx'
#
#PLK_ALL='/data1/rubinov_lab/brain_genomics/analyses_HCP/heritability/plink_files_full/all'
#RES_ALL='/data1/rubinov_lab/brain_genomics/analyses_HCP/heritability/results_gcta_fMPM/all'

## regional phenotype loops   
for reg in ${!REGS[*]}
do
    for phn in ${!PHNS[*]}
    do    
        region=${REGS[$reg]}
        phenotype=${PHNS[${phn}]}

        ## gcta greml (all) 
        ./gcta64 \
            --grm ${PLK_ALL}/chr_all/chr_all \
            --pheno ${PHN_DIR}/${phenotype}_${region}.phen \
            --reml \
            --out ${RES_ALL}/${phenotype}_${region}

        ## gcta greml (pdx) 
        ./gcta64 \
            --grm ${PLK_PDX}/chr_all/chr_all \
            --pheno ${PHN_DIR}/${phenotype}_${region}.phen \
            --reml \
            --out ${RES_PDX}/${phenotype}_${region}
    done
done

