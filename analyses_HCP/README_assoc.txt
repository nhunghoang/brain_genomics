Notes from Nhung (updated 10.10.2021) 

Running log of all the steps and scripts for the single/multi gene association analyses. 
nb: Assume all paths start with '/data1/rubinov_lab/brain_genomics/'
nb: All of this is being applied to the HCP cohort 

Gene expression: 
    + data_HCP/expression/filtered_quality_r0.3_p0.01/cohort890/
Phenotypes: 
    + data_HCP/hoacer_sn_hth/phenotypes/ 

[0] Store auxilary information 
[A] Generate a train/test group split (maintaining family groups)  
    > analyses_HCP/write_assoc_split.py 
    + analyses_HCP/DATA_OUTPUT/train_test_assoc_split.hdf5 
[B] Store twin/non-twin indices for whole/training group 
    > analyses/write_twin_indices.py  
    + analyses_HCP/DATA_OUTPUT/[train_]twin_indices.hdf5 
[C] Generate data files needed for EIGENSTRAT (all SNPs, whole cohort) 
    > analyses_HCP/write_eigendata.py 
    + analyses_HCP/DATA_OUTPUT/eigendata/hcp_cohort.[snp/geno/ind]
[D] Generate another .ind file that differentiates train/test populations 
    > analyses_HCP/write_new_eigen_ind.py
    + analysis_HCP/DATA_OUTPUT/eigendata/hcp_split.ind
    + analysis_HCP/DATA_OUTPUT/eigendata/[train/test]_pop.txt

** Apply [1] on train/test groups separately ** 
[1] Regress out confounders (age, gender, PC1, PC2) from expression and phenotypes 
[A] Compute PCs using EIGENSTRAT  
    > analyses_HCP/EIG/EIGENSTRAT/run_eigenstrat.perl 
    + analyses_HCP/DATA_OUTPUT/eigen_results/[train/test]_cohort.pca 
[B] Assemble covariate data 
    > analyses_HCP/write_covariates.py 
    + analyses_HCP/DATA_OUTPUT/[train/test]_covariates.txt 
[B1] (optional) Create PCA plots and cov-vs-phen plots 
    > analyses_HCP/plot_pca.py 
    + analyses_HCP/DATA_OUTPUT/pca_[train/test].png  
    > analyses_HCP/ #TODO 
[C] Calculate expression and phenotype residuals 
> analyses_HCP/write_residuals.py 

** Apply [2] on train group ** 
[2] Compute single-gene associations  
> analyses_HCP/new_single/single_select_genes.py 
[B] Plot p/fdr distributions and store [p <= 0.05] significant genes
> [...] 

[3] Compute multi-gene associations 
[A] Run gene set search on the observed data 
> new_multi/pool_run_phen_sa.py 
[B] Run gene set search on permutations 
> new_multi/pool_null_run_phen_sa.py

