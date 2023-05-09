'''
'''

import numpy as np 
import h5py 
import os 

## paths 
path_main = '/data1/rubinov_lab/brain_genomics/analyses_HCP/enrichment_1M'
path_genes = '/data1/rubinov_lab/brain_genomics/analyses_HCP/enrichment_1M/gene_sets'
path_syms = '/data1/rubinov_lab/brain_genomics/analyses_HCP/enrichment_1M/pdx_symbol_map.hdf5' 
path_out = '/data1/rubinov_lab/brain_genomics/analyses_HCP/enrichment_1M/counts' ## _[GO/HPO]
if not os.path.exists(path_out + '_GO'): 
    os.mkdir(path_out + '_GO')
    os.mkdir(path_out + '_HPO')

## names 
ont_files = {'GO': 'GO_Biological_Process_2021.txt', 'HPO': 'Human_Phenotype_Ontology.txt'}
phenotypes = ['gm_volume', 'myelination', 'alff', 'reho_noGS', 'connmean_noGS']
regions = ['caudate', 'nucleus-accumbens', 'hippocampus', 'amygdala', 'anterior-cingulate', \
    'putamen', 'hypothalamus', 'substantia-nigra', 'frontal-pole', 'cerebellar-hemisphere']

## DNE genes (no symbol) 
dne_genes = ['ENSG00000173209', 'ENSG00000179979', 'ENSG00000180279', \
    'ENSG00000213077', 'ENSG00000214189', 'ENSG00000214776', 'ENSG00000226067', \
    'ENSG00000230606', 'ENSG00000233276', 'ENSG00000269713', 'ENSG00000272710', \
    'ENSG00000273036', 'ENSG00000273136', 'ENSG00000276128', 'ENSG00000276399', \
    'ENSG00000277863', 'ENSG00000278662', 'ENSG00000283033', 'ENSG00000185684'] 
symbols = ['AHSA2P', 'CRIPAK', 'LINC01869', 'OSTCL', 'ZNF788P', '', \
    'LINC00623', '', 'GPX1', 'NBPF9', '', '', 'NBPF26', 'GXYLT1P5', \
    'FLJ36000', '', 'GOLGA6L10', '', 'EP400P1']  

## parse gene names (gene-to-symbol mapping) 
with h5py.File(path_syms, 'r') as f:
    info = np.array(f['mapping']).astype(str) ## array([sym-biomart, sym-pdx, gene]) 
    sym_map = {i[2]:i[0] for i in info} ## k: gene, v: sym-biomart 
    alt_map = {i[0]:i[1] for i in info} ## k: sym-biomart, v: sym-pdx 
    for gene, sym in zip(dne_genes, symbols): 
        sym_map[gene] = sym 
        alt_map[sym] = 'NONE' 

## parse ontologies (gene-to-terms mapping) 
ont_maps = {}  
for ont in ['GO', 'HPO']: 
    mapping = {} ## k: gene symbol, v: terms    
    mpath = '{}/{}'.format(path_main, ont_files[ont])
    with open(mpath, 'r') as f: lines = f.readlines() 
    for line in lines: 
        info = line.strip().split('\t')
        term = info[0] 
        genes = info[2:] ## two tabs between term and genes 
        for gene in genes: 
            try: mapping[gene].append(term) 
            except KeyError: mapping[gene] = [term] 
    ont_maps[ont] = mapping 

## gather gene sets of interest and convert to symbols  
phen_obsv = {} ## k: phen, v: gene symbols 
phen_null = {} ## k: phen, v: [gene symbols] 
nperms = int(1e4) 
for phen in phenotypes: 

    pfile = '{}/{}.hdf5'.format(path_genes, phen) 
    with h5py.File(pfile, 'r') as f: 

        ## observed genes 
        obsv_set = np.array(f['obsv-genes']).astype(str)
        phen_obsv[phen] = [sym_map[g] for g in obsv_set] 

        ## null genes 
        null_sets = [] 
        for i in range(nperms): 
            null_set = np.array(f['null-genes-{}'.format(i)]).astype(str)  
            null_sets.append([sym_map[g] for g in null_set])
        phen_null[phen] = null_sets 

## ontology loop 
for ont in ['GO', 'HPO']: 
    ont_map = ont_maps[ont] ## k: gene, v: terms  
    ont_keys = list(ont_map.keys())
    for phen in phenotypes:
        dne = 0 
        b_miss = 0 
        b_nomap = 0 
        p_miss = 0 
        p_nomap = 0 

        ## count observed genes per term 
        oterms = [] 
        for gene in phen_obsv[phen]:   
            try_alt = False ## start with biomart 

            if gene[:4] in ['NONE', '']: try_alt = True; b_miss += 1 ## missing biomart 
            elif gene not in ont_keys: try_alt = True; b_nomap += 1 ## biomart no map  
            else: gterms = ont_map[gene] ## valid biomart  

            if try_alt: ## try predixvu  
                gene = alt_map[gene]  
                if gene[:4] in ['NONE', '']: dne += 1; p_miss += 1; continue ## missing predixvu 
                elif gene not in ont_keys: dne += 1; p_nomap += 1; continue ## predixvu no map 
                else: gterms = ont_map[gene] ## valid predixvu  

            oterms.extend(gterms) 
        obsv_terms, obsv_counts = np.unique(oterms, return_counts=True) 
        oterm_counts = {t:c for t,c in zip(obsv_terms, obsv_counts)} ## k: term, v: count  

        ## count null genes per term 
        nterm_counts = {} ## k: term, v: array of len(nperms) containing term gene counts  
        for i, null_genes in enumerate(phen_null[phen]): 
            for gene in null_genes: 
                try_alt = False 
                gterms = [] 
                
                if gene[:4] in ['NONE', '']: try_alt = True 
                elif gene not in ont_keys: try_alt = True 
                else: gterms = ont_map[gene] 
    
                if try_alt: 
                    gene = alt_map[gene] 
                    if gene[:4] in ['NONE', '']: try_alt = True 
                    elif gene not in ont_keys: try_alt = True 
                    else: gterms = ont_map[gene] 

                for term in gterms:
                    try: 
                        nterm_counts[term][i] += 1 
                    except KeyError: 
                        nterm_counts[term] = np.zeros(nperms) 
                        nterm_counts[term][i] = 1 

        ## write term counts to file 
        cfile = '{}_{}/{}.hdf5'.format(path_out, ont, phen)         
        with h5py.File(cfile, 'w') as f: 
            for term, ocount in oterm_counts.items(): 
                f['{}-obsv'.format(term)] = ocount 
                f['{}-null'.format(term)] = nterm_counts[term] 

        print('{} - {}'.format(phen, ont))
        print('{} bmiss, {} bnomap, {} pmiss, {} pnomap'.format(b_miss, b_nomap, p_miss, p_nomap))
        print('{} total genes of interest, {} mapped\n'.format(len(phen_obsv[phen]), len(phen_obsv[phen])-dne))

        #print('{} obsv terms ({}-{}) [{} dne]'.format(len(oterm_counts), phen, ont, np.unique(dne).size))
        #total = len(phen_obsv[phen])
        #found = np.intersect1d(phen_obsv[phen], ont_keys).size
        #diff = np.setdiff1d(phen_obsv[phen], ont_keys).size
        #print('{} found, {} missing, {} total ({} - {})'.format(found, diff, total, phen, ont))
        #print('{} total, {} dne ({} - {})'.format(total, dne, phen, ont))
            
