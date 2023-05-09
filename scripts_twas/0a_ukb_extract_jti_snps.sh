#!/bin/bash

## Extract UKB SNPs used in the JTI models. 
## Then convert the BGEN files to VCF format. 

## Nhung, April 2023 

path_top="/data1/rubinov_lab/brain_genomics"
path_jti="${path_top}/models_JTI/rsids_by_tissue/all_rsids.txt"
path_inn="${path_top}/data_UKB/downloads/genotypes_hg37"

path_bgn="${path_top}/scripts_twas/inputs_UKB/bgen_JTI"
path_vcf="${path_top}/scripts_twas/inputs_UKB/vcf_JTI"

for CHR in {22..1} 
do 
    plink2 \
        --bgen ${path_inn}/c${CHR}_filtered.bgen ref-first \
        --sample ${path_inn}/cohort_filtered.sample \
        --extract ${path_jti} \
        --export bgen-1.2 \
        --snps-only \
        --out ${path_bgn}/c${CHR}

    bgenix \
        -g ${path_bgn}/c${CHR}.bgen \
        -vcf \
        > ${path_vcf}/c${CHR}.vcf &  
done 
