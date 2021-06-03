#!/bin/bash

## Script for saving the subset of HCP snps that are needed for PrediXcan v8 

## Also saves the subset of PrediXcan-GTEx snps that are not available in HCP  
## NOTE: this part is commented out, as we don't need this information to filter 
##       out gene models anymore

## Nhung, May 2021 

ml GCCcore/.6.4.0 BEDTools/2.27.1

pdx_dir="/data1/rubinov_lab/brain_genomics/models_PrediXcan_v8/snps_by_chr" 
hcp_dir="/data1/rubinov_lab/brain_genomics/data_HCP/Marchini_phg000989/snps_by_chrom_hg38" 

header="/data1/rubinov_lab/brain_genomics/data_HCP/predixcan_v8/header_v8.vcf" 
hcp_in_pdx_dir="/data1/rubinov_lab/brain_genomics/data_HCP/predixcan_v8/vcf_format"
#pdx_notin_hcp_file="/data1/rubinov_lab/brain_genomics/models_PrediXcan_v8/pdx_notin_hcp.txt" 

for chr in {22..1}
do 
    ## create vcf file with header 
    cp ${header} ${hcp_in_pdx_dir}/chr${chr}.vcf
    echo "- CHR ${chr} -" 
    
    ## use bedtools to get the snp intersection between HCP and PrediXcan-GTEx (PDX) 
    ## awk: keep snp IF at same position AND (ref & alt match OR flipped but match) 
    ## cut: PDX-chr, PDX-pos, PDX-rsid, PDX-ref, PDX-alt, HCP-qual, HCP-filter, HCP-info, HCP-format, HCP-vals 
    ##      >> save this HCP snp subset in vcf
    bedtools intersect -wa -wb \
        -a ${pdx_dir}/chr${chr}.vcf \
        -b ${hcp_dir}/chr${chr}.vcf | \
        awk '{ if ($2==$10) { if (($4==$12 && $5==$13) || ($4==$13 && $5==$12)) {print} } }' | \
        cut -f 1-5,15-18,19- >> ${hcp_in_pdx_dir}/chr${chr}.vcf  
    echo "hcp intersect pdx" 

    ## use bedtools to get all the snps in PDX that aren't in HCP  
    ## flag -v: report those in A that are not in B 
    ## cut: chr, pos, rsid 
    #bedtools intersect \
    #    -a ${pdx_dir}/chr${chr}.vcf \
    #    -b ${hcp_dir}/chr${chr}.vcf \
    #    -v | \
    #    cut -f 1,2,3 >> ${pdx_notin_hcp_file} 
    #echo "pdx intersect-v hcp" 

    ## another way of getting PDX snps that aren't in HCP 
    ## yields the same set of snps 
    #bedtools subtract \
    #    -a ${pdx_dir}/chr${chr}.vcf \
    #    -b ${hcp_dir}/chr${chr}.vcf | \
    #    cut -f 1,2,3 >> ${pdx_notin_hcp_file} 
    #echo "pdx subtract hcp" 
done 
