Identify brain-specific clinical traits from PrediXVU that associate 
with selected genes (genes that associate with neural phenotypes). 

Two different p-values, based on: 
1) [null] genes selected in the permutations 
2) [top_null] top N genes (by p-value) selected in the permutations, 
   where N is the observed number of genes

Scripts: 
> gather_predixvu.py --> output dir: symbol_ensembl_region  
- Create phenotype files containing regionally significant genes. 

> search_predixvu.py --> output dir: phecodes_genes, need_run  
- Query clinical traits that associate with selected genes.  
- Track PrediXVU gene filenames that aren't found (need_run). 

> search_predixvu2.py --> output dir: phecodes_genes 
- This script is for genes whose expected filenames weren't found. 
- If another file for these genes are found, the corresponding 
  results files are updated accordingly. 

> compute_pvalue.py 
- Count the number of times that observed PrediXVU results 
  appear in the permutations (script can run for either 
  permutation type). 

Top N genes permutations: 
> sort_pvalues.py --> output dir: sorted_pvals 
- For each phenotype, sort all regional genes by 
  FDR and then save the table.  
