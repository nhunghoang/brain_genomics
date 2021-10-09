Notes from Nhung (updated 10.09.2021) 

Running log of all the steps and scripts for the single/multi gene association analyses. 
nb: Assume all paths start with '/data1/rubinov_lab/brain_genomics/'
nb: All of this is being ran for the HCP cohort 

[0] Store auxilary information 
[A] Generate a train/test group split (maintaining family groups)  
> analyses_HCP/write_assoc_split.py 
[B] Store twin/non-twin indices for whole/training group 
> analyses/write_twin_indices.py  

[1] Regress out confounders (age, gender, PC1, PC2) from expression and phenotypes 
[A] Generated .geno, .snp, .ind data files (for EIGENSTRAT) - all SNPs, just not PDX SNPs
> data_HCP/write_eigendata.py 
[B] Computed PCs using EIGENSTRAT  
> analyses_HCP/single_gene_assoc/EIG/EIGENSTRAT/run_eigenstrat.perl 
[C] Assembled covariate data 
> data_HCP/write_covariates.py 
[D] Calculated expression and phenotype residuals 
> analyses_HCP/write_residuals.py 

[2] Compute single-gene associations (over the training group) 
> analyses_HCP/new_single/single_select_genes.py 
[B] Plot p/fdr distributions and store [p <= 0.05] significant genes
> [...] 

[3] Compute multi-gene associations 
[A] Run gene set search on the observed data 
> new_multi/pool_run_phen_sa.py 
[B] Run gene set search on permutations 
> new_multi/pool_null_run_phen_sa.py

