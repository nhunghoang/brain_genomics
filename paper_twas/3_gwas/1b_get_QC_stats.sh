#! /bin/bash

## At this stage, we already have variants with MAF > 0.01 
## and INFO SCORE > 0.2 based on all 500k UKB subjects. 

## Now, further filter variants based on: 
## EDIT (hard-call 0.1 | MAC 5), (missing 0.05 | MAF 0.01 | HWE 1*10^-5)

## Oct 2023 

for CHR in {22..1}
do 
    ## save new chr files with hard call genotypes  
    plink2 --bgen chr_maf_info/c${CHR}.bgen ref-first \
           --sample chr_maf_info/c${CHR}.sample \
           --snps-only \
           --hard-call-threshold 0.1 \
           --make-pgen \
           --out chr_plink_init/hardcall_c${CHR} 

    ## apply filters for missingness, then MAF 
    plink2 --pfile chr_plink_init/hardcall_c${CHR} \
           --geno 0.05 \
           --maf 0.01 \
           --make-pgen \
           --out chr_plink_init/c${CHR} 

    ## remove duplicate SNPs (by id, keep one), then apply HWE filter 
    ## export as bgen
    plink2 --pfile chr_plink_init/c${CHR} \
           --rm-dup exclude-mismatch \
           --hwe 0.00001 \
           --export bgen-1.2 'bits=8' \
           --out chr_bgen_final/c${CHR}

    ## index bgen file 
    bgenix -g chr_bgen_final/c${CHR}.bgen -index -clobber &

    ## sanity check 
    #plink2 --bgen chr_bgen_final/c${CHR}.bgen ref-first \
    #       --sample chr_bgen_final/c${CHR}.sample \
    #       --freq \
    #       --out c${CHR} &
done 
