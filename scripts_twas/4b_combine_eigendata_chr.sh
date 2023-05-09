#!/bin/bash 

dset=UKB ## HCP or UKB 
path=/data1/rubinov_lab/brain_genomics/scripts_assoc_clean/inputs_${dset}/eigendata

for chr in {1..22} 
do
    cat ${path}/chr${chr}_${dset}_cohort.geno >> ${path}/${dset}_cohort.geno 
    cat ${path}/chr${chr}_${dset}_cohort.snp >> ${path}/${dset}_cohort.snp 
done 

mkdir ${path}/chr_files
mv ${path}/chr* ${path}/chr_files 
