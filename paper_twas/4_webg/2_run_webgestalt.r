#!/usr/bin/env Rscript

## don't forget to module load R

library(WebGestaltR)

regs <- list('putamen', 'caudate', 'cerebellar-hemisphere', 'nucleus-accumbens',  
             'hippocampus', 'amygdala', 'dlpfc', 'anterior-cingulate')

group <- 'UKB'
ptype <- 'BON'

#phens <- 'vol_mean' 
#ontol <- 'gwas_catalog' 

phens <- 'vol_mean_interreg' 
ontol <- 'pdx_nom0.001'

topPath <- '/data1/rubinov_lab/brain_genomics/'
gmtPath <- paste(topPath, 'paper_twas/aux_files/', ontol, '.gmt', sep='')
desPath <- paste(topPath, 'paper_twas/aux_files/', ontol, '.des', sep='')

jtiPath <- paste(topPath, 'models_JTI/genes_by_tissue/', sep='')

innPath <- paste(topPath, 'paper_twas/inputs_', group, '/enrich_sets/', ptype, '_', phens, '_', sep='')
outPath <- paste(topPath, 'paper_twas/outputs_', group, '/enrich_', ontol, '/', ptype, '_', phens, sep='')

for (reg in regs) {

    interreg <- grepl('interreg', phens)
    if (interreg) {
        refsPath <- paste(jtiPath, 'genes_8regs.txt', sep='')
    } else {
        refsPath <- paste(jtiPath, reg, '.txt', sep='')
    }
    refsList <- readLines(refsPath)

    genePath <- paste(innPath, reg, '.txt', sep='')
    geneList <- readLines(genePath)

    enrichResult <- WebGestaltR(
        enrichMethod='ORA', 
        organism='hsapiens',
        enrichDatabaseFile=gmtPath, 
        enrichDatabaseType='ensembl_gene_id',
        enrichDatabaseDescriptionFile=desPath,
        interestGene=geneList,
        referenceGene=refsList,
        interestGeneType='ensembl_gene_id',
        referenceGeneType='ensembl_gene_id',
        sigMethod='top',
        topThr=6000,
        isOutput=TRUE,
        outputDirectory=outPath, projectName=reg)

}
