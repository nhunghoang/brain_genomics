UKB: 

    See 'brain_genomics/data_UKB/downloads/README.txt' for 
    details on getting covariates, phenotypes, and genotypes
    for individuals with White British ancestry. 

HCP: 

    See 'brain_genomics/data_HCP/ancestry_clusters/' for 
    details on applying PCA on the genotypes and then 
    filtering for European ancestry. 

################################################################

gr-Expression inference:

    Extract SNPs that are in the JTI models: 
        0_ukb_extract_jti_snps.sh
        0_hcp_extract_jti_snps.py

    Convert genotype probabilities to dosage format: 
        1_convert_to_dosage.py

    Infer gr-expression: 
        2_infer_grex.sh
