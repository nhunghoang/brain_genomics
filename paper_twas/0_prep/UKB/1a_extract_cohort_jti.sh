#!/bin/bash

## Extract SNPs used in the JTI models for the  UKB cohort 
## of interest. Then convert the BGEN files to VCF format. 

## Nhung, Feb 2024 

coho="UKB/nonBrit" # UKB/nonBrit

path_top="/data1/rubinov_lab/brain_genomics"
path_jti="${path_top}/models_JTI/rsids_by_tissue/all_rsids.txt"
path_inn="${path_top}/data_UKB/downloads/full_genotypes_hg37"

path_sample="${path_top}/data_UKB/downloads/archives/ukb22828_all.sample"
path_cohort="${path_top}/paper_twas/inputs_${coho}/cohort.txt"

path_bgn="${path_top}/paper_twas/inputs_${coho}/bgen_JTI"
path_vcf="${path_top}/paper_twas/inputs_${coho}/vcf_JTI"

for CHR in {22..1} 
do 
    if [[ $CHR == 16 || $CHR == 10 || $CHR == 4 || $CHR == 1 ]]; then
        plink2 \
            --bgen ${path_inn}/ukb22828_c${CHR}_b0_v3.bgen ref-first \
            --sample ${path_sample} \
            --keep ${path_cohort} \
            --extract ${path_jti} \
            --export bgen-1.2 \
            --snps-only \
            --out ${path_bgn}/c${CHR}  
    else
        plink2 \
            --bgen ${path_inn}/ukb22828_c${CHR}_b0_v3.bgen ref-first \
            --sample ${path_sample} \
            --keep ${path_cohort} \
            --extract ${path_jti} \
            --export bgen-1.2 \
            --snps-only \
            --out ${path_bgn}/c${CHR} &
    fi
done

for CHR in {22..1}
do
    bgenix \
        -g ${path_bgn}/c${CHR}.bgen \
        -index 

    bgenix \
        -g ${path_bgn}/c${CHR}.bgen \
        -vcf \
        > ${path_vcf}/c${CHR}.vcf & 
done
