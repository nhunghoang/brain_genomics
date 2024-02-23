'''
Script for applying K-Means to HCP PCA results, 
in order to define the Euro cluster. 

- Nhung, updated Jan 2024
'''

import numpy as np 
import pandas as pd
from sklearn.cluster import KMeans

## paths 
top_path = '/data1/rubinov_lab/brain_genomics/paper_twas/inputs_HCP/'
pca_path = top_path + 'eigen_output/jti2.pca.evec'
ids_path = top_path + 'demographics.csv'
out_path = top_path + 'ancestry2.csv'

## parse PCA results 
with open(pca_path, 'r') as f: 
    _ = f.readline() 
    lines = [line.split() for line in f.readlines()] 
    lines = np.array(lines)

pca = pd.DataFrame({'sample': lines[:,0], 'self_reported': lines[:,-1]})
pca['pc1'] = lines[:,1].astype(float)
pca['pc2'] = lines[:,2].astype(float)

## run K-Means multiple times to check for stability 
km_clusters = []
for i in range(50):

    ## fit K-Means model
    model = KMeans(n_clusters=3) 
    model.fit(pca[['pc1', 'pc2']].values)

    ## largest cluster is Euro, followed by African, then Asian
    clusters, csize = np.unique(model.labels_, return_counts=True) 
    clusters = clusters[np.argsort(csize)]
    cluster_labels = {c: l for c,l in zip(clusters, ['Asian', 'African', 'White'])}

    ## set cluster labels
    pca['cluster_num'] = model.labels_
    pca['cluster'] = pca['cluster_num'].map(cluster_labels)
    km_clusters.append(pca['cluster'].values)

assert(np.all([np.array_equal(km_clusters[0], km_clusters[i]) for i in range(50)]))

## add subject IDs 
id_map = pd.read_csv(ids_path, index_col='sample').to_dict()['subject']
pca['subject'] = pca['sample'].map(id_map)

## save ancestry results  
cols = ['sample', 'subject', 'self_reported', 'cluster']
pca.to_csv(out_path, columns=cols, index=False)

#################################################################

## plot PCA with marker colors for self-reported ancestry
## and marker shapes for cluster ancestry  
def plot_pca():

    import matplotlib 
    matplotlib.use('tkagg') 
    import matplotlib.pyplot as plt 

    ## set up 2D/3D scatter plot 
    plt.ion()
    fig = plt.figure()
    ax = fig.add_subplot() #(projection='3d')
    ax = fig.subplots(1,1) 
     
    ## population marker params 
    color_dict = {'White':'r', 'Black':'g', 'Asian/Nat.':'b', \
                  'Unknown': 'k', 'Am.': 'm', 'More': 'y'}
    shape_dict = {'White': 'o', 'African': '*', 'Asian': 's'}

    ancestry['self_short'] = ancestry['self_reported'].apply(lambda x: x.split(' ')[0])
    ancestry['mcolor'] = ancestry['self_short'].map(color_dict)
    ancestry['mshape'] = ancestry['cluster'].map(shape_dict)

    ## link pca data 
    pc1_map = pca.set_index('sample_id').to_dict()['pc1']
    pc2_map = pca.set_index('sample_id').to_dict()['pc2']

    ancestry['pc1'] = ancestry['sample_id'].map(pc1_map) 
    ancestry['pc2'] = ancestry['sample_id'].map(pc2_map) 

    ## plot 
    sc = ax.scatter(ancestry['pc1'], \
                    ancestry['pc2'], \
                    c='none', \
                    s=40, \
                    edgecolor=ancestry['mcolor'], \
                    marker=ancestry['mshape'], \
                    alpha=0.3)
                    
#################################################################

