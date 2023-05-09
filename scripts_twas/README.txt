last update: May 8, 2023 (Nhung) 

******************************************************
Analyses for the 2022 paper. 

Note: almost all paths have this prefix: 
      '/data1/rubinov_lab/brain_genomics/scripts_twas'

Note: 'x' refers to input data;
      'f' refers to the script; 
      'y' refers to output data; 
      'cohort' refers to UKB or HCP 
******************************************************

GREx INFERENCE: 

 0) Extract SNPs from the original genotype files 
    for the JTI GREx models of interest. 

    HCP: (x) ../data_HCP/Marchini_phg000989/snps_by_chrom_hg37
         (f) 0a_hcp_extract_jti_snps.py

    UKB: (x) ../data_UKB/downloads/genotypes_hg37
         (f) 0a_ukb_extract_jti_snps.sh

    (y) inputs_{cohort}/vcf_JTI/c*.vcf  

 1) Convert genotype probabilities to dosages for JTI. 

    (x) inputs_{cohort}/vcf_JTI/c*.vcf  
    (f) 1_convert_to_dosage.py 
    (y) inputs_{cohort}/dosage_JTI/c*.dosage.txt 


 2) Infer GREx in all brain regions of interest. 

    (x) inputs_{cohort}/dosage_JTI/c*dosage.txt 
        inputs_{cohort}/{genotype_samples.txt OR cohort_filtered.txt}
    (f) 2_infer_grex.sh
    (y) outputs_{cohort}/grex_JTI/{region}.hdf5  

NEUROIMAGING TWAS:

 3) Apply any necessary preprocessing. 

    HCP: Compute eigenvectors as PCs   
         (x) ../data_HCP/Marchini_phg000989/snps_by_chrom_hg37/*vcf
         (f) 3a_hcp_write_eigen_input.py
         (xy) inputs_HCP/eigen_input/HCP_cohort.*
         (f) ../analyses_HCP/EIG/EIGENSTRAT/run_eigenstrat.perl 
         (y) inputs_HCP/eigen_output/HCP_cohort.pca

    Any: Format covariates of interest 
         (x) (locations of demographic data vary) 
         (f) 3b_format_covariates.py
         (y) inputs_{cohort}/covariates.csv 

    Any: Regress covariates from GREx and phenotypes 
         (x) inputs_{cohort}/phenotypes.csv 
             inputs_{cohort}/grex_JTI/*
             inputs_{cohort}/covariates.csv 
         (f) 3c_write_residuals.py
         (y) inputs_{cohort}/grex_JTI_regress/*
             inputs_{cohort}/phen_regress/*

 4) Run MetaXcan TWAS or permutation-based TWAS 

    (x) inputs_{cohort}/phenotypes.csv
        inputs_{cohort}/covariates.csv 
        outputs_{cohort}/grex_JTI/* 

    MetaXcan:
    (f) 4a_run_metaxcan_assoc.sh 
    (y) outputs_{cohort}/twas_metaxcan_JTI/*

    Permutation-based: 
    (f) 4b_run_1M_perms_assoc.py
    (y) outputs_{cohort}/twas_1M_perms_JTI/*

GENE-SET ANALYSES: 

 5) Generate TWAS permutations. 

    (x) inputs_{cohort}/set_permutations.hdf5 
        inputs_{cohort}/grex_JTI_regress/*
        inputs_{cohort}/phen_regress/*
    (f) 5a_run_twas_perms.py
    (y) outputs_{cohort}/permutation_twas_JTI

 6) Construct similarity matrices based on TWAS results. 
    (x) outputs_{cohort}/twas_1M_perms_JTI
        outputs_{cohort}/permutation_twas_JTI 
    (f) 6_count_regphen_similarities.py
    (y) outputs_{cohort}/regphen_similarities_JTI/* 
 
 7) Compute [gene set, phenotype] associations. 

    Based on gene rankings in TWAS:
    (f) 7a_multigene_association.py 
    (y) outputs_{cohort}/multigene_ranked_sets_JTI 

    Based on linear regression: 
    (f) 7b_multigene_regression.py
    (y) outputs_{cohort}/multigene_linear_sets_JTI 

 8) Run PrediXVU enrichment analysis. 
    >> scripts_twas/11a_predixvu_filenames.py
    (output) models_PrediXcan_v8/predixvu_ens_sym_map.hdf5
    (output) models_PrediXcan_v8/predixvu_ens_sym_dne.txt
    >> scripts_twas/11b_predixvu_assocs.py 
    (output) models_PrediXcan_v8/predixvu_assocs.csv
    >> scripts_twas/11c_predixvu_enrichment.py 
    (output) scripts_twas/outputs_HCP/predixvu/cloud_results.hdf5
    (output) scripts_twas/outputs_HCP/predixvu/set_enrichments.hdf5


 9) Run GO/HPO enrichment analysis. 

