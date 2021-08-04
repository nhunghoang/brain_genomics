#!/usr/bin/env bash

VCF_ALL='/data1/rubinov_lab/brain_genomics/data_HCP/Marchini_phg000989/snps_by_chrom_hg38'
VCF_PDX='/data1/rubinov_lab/brain_genomics/data_HCP/predixcan_v8/vcf_format' 

PLK_ALL='/data1/rubinov_lab/brain_genomics/analyses_HCP/heritability/plink_files_full/all'
PLK_PDX='/data1/rubinov_lab/brain_genomics/analyses_HCP/heritability/plink_files_full/pdx'

GRM_ALL='/data1/rubinov_lab/brain_genomics/analyses_HCP/heritability/grm_chrs_all_full.txt'
GRM_PDX='/data1/rubinov_lab/brain_genomics/analyses_HCP/heritability/grm_chrs_pdx_full.txt'

#for chr in {9..1}
#do
#    ## vcf to plink (all)
#    ./plink2 --vcf ${VCF_ALL}/chr${chr}.vcf dosage=GP --out ${PLK_ALL}/chr${chr}/chr${chr}
#    python update_psam.py ${PLK_ALL}/chr${chr}/chr${chr}.psam
#
#    ## vcf to plink (pdx)
#    ./plink2 --vcf ${VCF_PDX}/chr${chr}.vcf dosage=GP --out ${PLK_PDX}/chr${chr}/chr${chr}
#    python update_psam.py ${PLK_PDX}/chr${chr}/chr${chr}.psam
#
#    ## gcta grm (all)
#    ./gcta64 \
#        --pfile ${PLK_ALL}/chr${chr}/chr${chr} \
#        --make-grm \
#        --out ${PLK_ALL}/chr${chr}/chr${chr}
#
#    ## gcta grm (pdx)
#    ./gcta64 \
#        --pfile ${PLK_PDX}/chr${chr}/chr${chr} \
#        --make-grm \
#        --out ${PLK_PDX}/chr${chr}/chr${chr}
#done

## merge chr files
./gcta64 --mgrm ${GRM_ALL} --make-grm --out ${PLK_ALL}/chr_all/chr_all
./gcta64 --mgrm ${GRM_PDX} --make-grm --out ${PLK_PDX}/chr_all/chr_all
