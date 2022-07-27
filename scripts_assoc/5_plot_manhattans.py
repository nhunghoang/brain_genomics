'''
'''

import numpy as np 
import matplotlib.pyplot as plt 
import os 

## paths 
assoc_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/assoc_10k/pvals' ## _[phen]/[reg].txt 
genes_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/assoc/paper_figures/genes_parsed.bed' 
out_path = '/data1/rubinov_lab/brain_genomics/analyses_HCP/assoc_10k/manh_plots'
if not os.path.exists(out_path): os.mkdir(out_path) 

## regional phenotypes 
phens = ['gm_volume', 'alff', 'reho_noGS', 'connmean_noGS', 'myelination']
regs = ['hippocampus', 'amygdala', 'hypothalamus', 'substantia-nigra',\
        'caudate', 'putamen', 'nucleus-accumbens', 'anterior-cingulate',\
        'frontal-pole', 'cerebellar-hemisphere']

## parse gene ordering 
data = np.loadtxt(genes_path, skiprows=1, usecols=[0,3,4], dtype=str)  
ranks = {d[1]: (int(d[2]), int(d[0])) for \
        d in data if d[0][0] not in ['K', 'X']} ## k: gene, v: (rank, chr)   

## function: manhattan plot 
def manhattan(data_by_chr, out_file): 
    colors = ['#bbbbbb', '#8a8a8a'] 
    fig, ax = plt.subplots(1,1, figsize=(6,6)) 
    start_idx = 0 

    chr_ticks = [] 
    lbar = -np.log10(0.005)
    for ch, info in data_by_chr.items(): 
        corrs = info[0][:,0]
        logps = info[0][:,1]
        mcolors = info[1]
        X = range(start_idx, start_idx + logps.size)  
        sizes = (1e4 * (corrs**2)).astype(int) + 5
        color = [colors[ch%2] if logps[i] < lbar else mcolors[i] for i in range(logps.size)] 
        ax.scatter(X, logps, s=35, color=color, marker='.', alpha=0.40, linewidth=1) 
        #ax.scatter(X, logps, s=40, color=color, marker='.') 
        chr_ticks.append(start_idx + (logps.size/2) - 1)
        start_idx += logps.size 

    #ps = [0.1, 0.15, 0.2, 0.25] 
    ps = [bar] 
    lp = [-1*np.log10(p) for p in ps] 
    for p, l in zip(ps, lp):
        #print('{}: {:.2f}'.format(p, l))
        ax.hlines(l, xmin=0, xmax=start_idx, linestyles='dashed', linewidth=0.5, color='k')
    ax.set_ylim([0,4.25])

    ax.set_xticks(chr_ticks) 
    ax.set_xticklabels([str(c) for c in range(1,23)], size=6)

    #ax.set_yticks(lp)
    #ax.set_yticklabels(['p = {:.2f}'.format(p) for p in ps], size=8)

    ax.set_xlabel('chromosome', fontsize=10)
    ax.set_ylabel('$-log_{10}(p)$', fontsize=10) 

    #xleft, xright = ax.get_xlim()
    #ybottom, ytop = ax.get_ylim()
    #ax.set_aspect(abs((xright-xleft)/((ybottom-ytop)*3)))
    ax.set_box_aspect(1/3) ## height = aspect * width

    plt.margins(x=0.0, y=0.01)
    plt.tight_layout()
    plt.savefig(out_file, dpi=400) 
    plt.close('all') 

## regional marker colors 
rcolors = ['#FE3A3E', '#F1A043', '#6FAE61', '#2F6EBA', '#5F2F8D', \
           '#EC7F82', '#F2C164', '#9FCF63', '#77B3D6', '#B176F7']
rcolors = ['#FE3A3E', '#EC7F82','#F1A043', '#F2C164','#6FAE61', \
           '#9FCF63','#2F6EBA', '#77B3D6','#5F2F8D', '#B176F7']
rcolors.reverse()

## parse association data 
bar = 0.005 ## uncorrected 
pcol = 2 ## 2 for uncorrected, 3 for FDR-corrected 
for phen in phens:
    print('\n' + phen)
    adata_by_chr = {c:[] for c in range(1,23)} ## k: chrom, v: genes * ([rank, rho, logp], color)
    for r, reg in enumerate(regs): 

        afile = '{}_{}/{}.txt'.format(assoc_path, phen, reg)
        adata = np.loadtxt(afile, skiprows=1, usecols=[0,1,pcol], dtype=str)  
        agenes = adata[:,0] 
        acorrs = adata[:,1].astype(float) 

        valid_idx = ~np.isnan(acorrs)
        acorrs = np.abs(acorrs[valid_idx]) 
        agenes = agenes[valid_idx] 

        apvals = adata[:,2].astype(float)[valid_idx] 
        apvals = np.where(apvals <= 0.0001, 0.0001, apvals)
        alogps = -1 * np.log10(apvals) 

        print((apvals <= bar).sum(), reg)

        for g,gene in enumerate(agenes):
            corr = acorrs[g] 
            logp = alogps[g]
            try: 
                (rank, chrm) = ranks[gene] 
                adata_by_chr[chrm].append(([rank, corr, logp], rcolors[r])) 
            except: 
                #print('{}: {:.2f}'.format(gene, apvals[g])) 
                continue 

    ## sort genes by rank 
    data_sorted = {} ## k: chrom, v: ([rank-sorted rhos, pvals], colors) 
    for chrm, info in adata_by_chr.items(): 
        info1 = np.array([i[0] for i in info])
        info2 = np.array([i[1] for i in info])
        idx = np.argsort(info1[:,0])
        data_sorted[chrm] = (info1[:,1:][idx], info2[idx])
        
    ofile = '{}/{}.png'.format(out_path, phen) 
    manhattan(data_sorted, ofile)  

