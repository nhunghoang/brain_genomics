Methods for the 2022 paper. 
- Nhung (updated Sept 19, 2022) 

******************************************************************************
Genotype preprocessing: 
HCP: HCP (hg37) already comes in VCF format 
     --> LiftOver to hg38 
     --> convert SNPs to dosages (only for the subset of SNPs available for PrediXcan, or convert all SNPs?) 
     --> run PrediXcan 
UKB: 
******************************************************************************
0a) Impute gene expression (via PrediXcan-GTEX, filtered by prediction quality) 
>> /data1/rubinov_lab/brain_genomics/scripts_assoc/0a_run_PrediXcan_v8.sh

0b) Format phenotype data 
>> /data1/rubinov_lab/brain_genomics/scripts_assoc/0b_format_phenotypes.py  

1) Compute genotype PCs using EIGENSTRAT 
>> /data1/rubinov_lab/brain_genomics/scripts_assoc/1a_write_eigendata.py
>> /data1/rubinov_lab/brain_genomics/scripts_assoc/1b_combine_eigendata_chr.sh 
>> /data1/rubinov_lab/brain_genomics/analyses_HCP/EIG/EIGENSTRAT/run_eigenstrat.perl

2) Regress confounders (age, gender, PC1, PC2) from every gene and phenotype 
>> /data1/rubinov_lab/brain_genomics/scripts_assoc/2a_write_covariates.py
>> /data1/rubinov_lab/brain_genomics/scripts_assoc/2b_write_residuals.py

3) Compute gene-phenotype associations 
>> /data1/rubinov_lab/brain_genomics/scripts_assoc/3a_hcp_write_twin_perms.py 
>> /data1/rubinov_lab/brain_genomics/scripts_assoc/3a_ukb_write_perms.py 
>> /data1/rubinov_lab/brain_genomics/scripts_assoc/3b_run_associations.py 

4) Compute permutation gene-phenotype associations 
>> /data1/rubinov_lab/brain_genomics/scripts_assoc/4a_null_associations.py
>> /data1/rubinov_lab/brain_genomics/scripts_assoc/4b_gather_nulls.py 

5) Construct similarity matrices based on association results  
>> /data1/rubinov_lab/brain_genomics/scripts_assoc/5_regphen_overlaps.py

6) Compute gene-set associations, including nulls 
>> /data1/rubinov_lab/brain_genomics/scripts_assoc/6_compute_set_associations.py

3a) Run PrediXVU enrichment analysis 

3b) Run GO/HPO enrichment analysis  


******************************************************************************

Scripts for generating figures: 

