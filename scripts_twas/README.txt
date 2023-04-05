last update: March 10, 2023 (Nhung) 

******************************************
Analyses for the 2022 paper. 

Note: all paths have this prefix: 
      '/data1/rubinov_lab/brain_genomics/'
******************************************

 0) Lift the HCP genotypes from hg37 to hg38.
    (input) data_HCP/Marchini_phg000989/snps_by_chrom_hg37
    >> [script??] 
    (output) data_HCP/Marchini_phg000989/snps_by_chrom_hg38
 
 1) Generate dosage files for the cohort genotypes.   
    >> scripts_twas/1_generate_dosage_files.py
    (output) scripts_twas/inputs_HCP/predixcan_samples.txt
    (output) scripts_twas/inputs_HCP/dosages_hg38
 
 2) Infer GReX in all ten brain regions.  
    >> /data1/rubinov_lab/brain_genomics/scripts_twas/2_run_PrediXcan_v8.sh
    (output) scripts_twas/inputs_HCP/predixcan_grex_v8
 
 3) Format phenotype data. 
    >> 3_format_phenotypes.py 
    (output) scripts_twas/inputs_HCP/phenotypes 
 
 4) Compute genotype PCs using EIGENSTRAT. 
    >> 4a_write_eigen_input.py 
    >> 4b_combine_eigendata_chr.sh
    (output) scripts_twas/inputs_HCP/eigen_input
    >> analyses_HCP/EIG/EIGENSTRAT/run_eigenstrat.perl
    (output) scripts_twas/inputs_HCP/eigen_output 
 
 5) Regress confounders (age, gender, PC1, PC2) from every gene and phenotype.  
    >> scripts_twas/5a_write_covariates.py 
    (output) scripts_twas/inputs_HCP/covariates.txt  
    >> scripts_twas/5b_write_residuals.py 
    (output) scripts_twas/inputs_HCP/expr_regress 
    (output) scripts_twas/inputs_HCP/phen_regress 
 
 6) Compute gene-phenotype associations. 
    >> scripts_twas/6a_hcp_write_perms.py 
    (output) scripts_twas/inputs_HCP/permutations_100k.hdf5
    >> scripts_twas/6b_run_associations.py 
    (output) scripts_twas/outputs_HCP/assoc_1M
 
 7) Compute permutation gene-phenotype associations. 
    >> scripts_twas/7a_run_null_associations.py
    (output) scripts_twas/outputs_HCP/assoc_1M/nulls
    >> scripts_twas/7b_gather_nulls.py 
 
 8) Construct similarity matrices based on association results. 
    >> scripts_twas/8_regphen_similarities.py
    (output) scripts_twas/outputs_HCP/regphen_similarities/mats_inter[reg/phen].hdf5 
 
 9) Compute [gene set, phenotype] associations using gene rankings. 
    >> scripts_twas/9_multigene_association.py 
    (output) scripts_twas/outputs_HCP/multigene_ranked_sets 

10) Compute [gene set, phenotype] linear regressions. 
    >> scripts_twas/10_multigene_regression.py 
    (output) scripts_twas/outputs_HCP/multigene_linear_sets  

11) Run PrediXVU enrichment analysis. 
    >> scripts_twas/11a_predixvu_filenames.py
    (output) models_PrediXcan_v8/predixvu_ens_sym_map.hdf5
    (output) models_PrediXcan_v8/predixvu_ens_sym_dne.txt
    >> scripts_twas/11b_predixvu_assocs.py 
    (output) models_PrediXcan_v8/predixvu_assocs.csv
    >> scripts_twas/11c_predixvu_enrichment.py 
    (output) scripts_twas/outputs_HCP/predixvu/cloud_results.hdf5
    (output) scripts_twas/outputs_HCP/predixvu/set_enrichments.hdf5


12) Run GO/HPO enrichment analysis. 

13) Run UKB replication analysis (Mika's scripts). 


### Figure Generation Scripts ### 

