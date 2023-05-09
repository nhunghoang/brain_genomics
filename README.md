### Summary of software

This repository contains analysis software that accompanies Hoang N, Sardaripour N, Ramey G, Chen Y, Liao E, Park JH, Benton ML, Capra JA, Rubinov M. Genetically regulated gene expression underlies the individuality of human brain organization and disease risk. doi.org/10.31219/osf.io/xefru

### Main scripts

The project is organized into the following directories and subdirectories:

- `analyses_UKB` and `analyses_HCP`: Main data analyses.
- `scripts_assoc`: Main association analyses.
- `scripts_data_agnostic`: Processing of phenotypes and eigendata.
- `scripts_neuro_prepro`: Processing of neuroimaging data.
- `models_PrediXcan_v8`: Processing of PrediXcan v8 models.
- `data_Allen`: Processing Allen Institute gene-expression data.
- `data_HCP`: Processing Human Connectome Project data.

### Main analysis

0. Preprocessing
    - Impute gene expression
    - Format phenotype data
1. Compute genotype PCs using EIGENSTRAT
2. Regress confounders from every gene and phenotype
3. Compute gene-phenotype associations
4. Compute permutation gene-phenotype associations
5. Construct similarity matrices based on association results
6. Compute gene-set associations, including nulls
7. Run PrediXVU enrichment analysis
