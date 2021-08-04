#!/usr/bin/env bash

VCF_DIR='/data1/rubinov_lab/brain_genomics/data_HCP/Marchini_phg000989/snps_by_chrom_hg37'
PHN_DIR='/data1/rubinov_lab/brain_genomics/analyses_HCP/heritability_v7/phen_files'

REGS=(amygdala aCingulateCortex caudateBG cerebellum frontalCortex hippocampus nAccumbensBG putamenBG)
PHNS=(alff falff connvar)

#mkdir plink_files
#mkdir plink_files/chr{1..22}
#mkdir plink_files/chr_all
#
#mkdir results_gcta
#
#for chr in {22..1}
#do
#    ## vcf to plink
#    ./plink2 --vcf ${VCF_DIR}/chr${chr}.filtered.sampid.vcf dosage=GP --out plink_files/chr${chr}/chr${chr}
#    python update_psam.py plink_files/chr${chr}/chr${chr}.psam
#
#    ## GCTA GRM
#    ./gcta64 \
#        --pfile plink_files/chr${chr}/chr${chr} \
#        --make-grm \
#        --out plink_files/chr${chr}/chr${chr}
#done
#
### merge chr files
#./gcta64 --mgrm grm_chrs.txt --make-grm --out plink_files/chr_all/chr_all

for reg in ${!REGS[*]}
do
    for phn in ${!PHNS[*]}
    do
        region=${REGS[$reg]}
        phenotype=${PHNS[${phn}]}

        ## GCTA GREML -- all chr
        ./gcta64 \
            --grm plink_files/chr_all/chr_all \
            --pheno ${PHN_DIR}/${phenotype}_${region}.phen \
            --reml \
            --out results_gcta/${phenotype}_${region}
    done
done
