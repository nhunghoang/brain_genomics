#!/usr/bin/env bash

VCF_DIR='/data1/rubinov_lab/brain_genomics_accre/data_MarchiniPrediXcan_v7/vcf_format'
PHN_DIR='/data1/rubinov_lab/brain_genomics/analyses_HCP/heritability_v7/phen_files'

REGS=(amygdala aCingulateCortex caudateBG cerebellum frontalCortex hippocampus nAccumbensBG putamenBG)
PHNS=(falff)

#mkdir subset_plink_files
#mkdir subset_plink_files/chr{1..22}
#mkdir subset_plink_files/chr_all
#
#mkdir subset_results_gcta
#
#for chr in {22..1}
#do
#    ## vcf to plink
#    ./plink2 --vcf ${VCF_DIR}/chr${chr}_subset.vcf dosage=GP --out subset_plink_files/chr${chr}/chr${chr}
#    python update_psam.py subset_plink_files/chr${chr}/chr${chr}.psam
#
#    ## GCTA GRM
#    ./gcta64 \
#        --pfile subset_plink_files/chr${chr}/chr${chr} \
#        --make-grm \
#        --out subset_plink_files/chr${chr}/chr${chr}
#done
#
### merge chr files
#./gcta64 --mgrm subset_grm_chrs.txt --make-grm --out subset_plink_files/chr_all/chr_all

for reg in ${!REGS[*]}
do
    for phn in ${!PHNS[*]}
    do
        region=${REGS[$reg]}
        phenotype=${PHNS[${phn}]}

        ## GCTA GREML -- all chr
        ./gcta64 \
            --grm subset_plink_files/chr_all/chr_all \
            --pheno ${PHN_DIR}/${phenotype}_${region}.phen \
            --reml \
            --out subset_results_gcta/${phenotype}_${region}
    done
done
