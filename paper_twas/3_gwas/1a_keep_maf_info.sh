#!/bin/bash

## Extract UKB SNPs with filters MAF > 0.01 and INFO SCORE > 0.2. 

## Nhung, Oct 2023 

path_top="/data1/rubinov_lab/brain_genomics"
path_jti="${path_top}/data_UKB/downloads/imputations_maf"
path_inn="${path_top}/data_UKB/downloads/genotypes_hg37"

path_out="${path_top}/analyses_UKB/qual_ctrl/chrs_maf_info"

for aCHR in {1..11}
do 
    plink2 \
        --bgen ${path_inn}/c${aCHR}_filtered.bgen ref-first \
        --sample ${path_inn}/cohort_filtered.sample \
        --extract ${path_jti}/rsids_keep/rsids_c${aCHR}.txt \
        --export bgen-1.2 \
        --snps-only \
        --out ${path_out}/c${aCHR} &

    bCHR=$((23-${aCHR}))
    plink2 \
        --bgen ${path_inn}/c${bCHR}_filtered.bgen ref-first \
        --sample ${path_inn}/cohort_filtered.sample \
        --extract ${path_jti}/rsids_keep/rsids_c${bCHR}.txt \
        --export bgen-1.2 \
        --snps-only \
        --out ${path_out}/c${bCHR}
done 

for CHR in {22..1}
do 
    bgenix \
        -g ${path_out}/c${CHR}.bgen \
        -index & 
done
