Methods for the 2022 paper. 

Genotype preprocessing: 
* HCP: HCP (hg37) already comes in VCF format --> LiftOver to hg38 --> convert SNPs to dosages (only for the subset of SNPs available for PrediXcan, or convert all SNPs?) --> run PrediXcan
* UKB: ???

1) SNPs to GREx (via PrediXcan-GTEx, filtered by prediction quality) 
>> /data1/rubinov_lab/brain_genomics/scripts_data_agnostic/run_PrediXcan_v8.sh 

2a) Compute genotype PCs using EIGENSTRAT
>> /data1/rubinov_lab/brain_genomics/scripts_data_agnostic/[HCP/UKB]_write_eigendata.py  
>> /data1/rubinov_lab/brain_genomics/scripts_data_agnostic/combine_eigendata_chr.sh **change dset in script 
>> /data1/rubinov_lab/brain_genomics/analyses_HCP/EIG/EIGENSTRAT/run_eigenstrat.perl

2b) Regress confounders (age, gender, PC1, PC2) from each gene and phenotype 
>> /data1/rubinov_lab/brain_genomics/scripts_data_agnostic/[HCP/UKB]_write_covariates.py
   Note: separate scripts for HCP/UKB, but output covariate files should be in the same format
>> /data1/rubinov_lab/brain_genomics/scripts_data_agnostic/write_residuals.py [HCP/UKB]

3) Genes that individually correlate with regional phenotypes  
>> /data1/rubinov_lab/brain_genomics/scripts_data_agnostic/[HCP/UKB]_write_assoc_perms.py
>> /data1/rubinov_lab/brain_genomics/scripts_data_agnostic/select_indep_genes.py [HCP/UKB]

not amended yet:
>> /data1/rubinov_lab/brain_genomics/scripts_data_agnostic/null_jobs_select_indep_genes.py 

4-1) Find gene sets that best predicts interregional phenotypes (e.g., CCC)  

4-2) PrediXVU, gene ontology, gene set enrichment analysis 
