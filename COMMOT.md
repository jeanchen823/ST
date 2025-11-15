# COMMOT



```bash
https://commot.readthedocs.io/en/latest/index.html

conda create -n commot1  python=3.10.14  -c https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge

conda activate commot1

pip install commot

conda install--channel conda-forge pygraphviz

pip install rpy2==3.4.2pip install anndata2ri==1.0.6

conda install spyder-kernels=2.5  -c https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge
conda install spyder-notebook  -c https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge

spyder

```



```bash
import os
import gc
import ot
import pickle
import anndata
import scanpy as sc
import pandas as pd
import numpy as np
from scipy import sparse
from scipy.stats import spearmanr, pearsonr
from scipy.spatial import distance_matrix
import matplotlib.pyplot as plt
import commot as ct


adata = sc.datasets.visium_sge(sample_id='V1_Mouse_Brain_Sagittal_Posterior')
adata.var_names_make_unique()
adata.raw = adata

sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)
adata_dis500 = adata.copy()

sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]

sc.tl.pca(adata, svd_solver='arpack')

sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)

sc.tl.leiden(adata, resolution=0.4)

sc.pl.spatial(adata, color='leiden')

df_cellchat = ct.pp.ligand_receptor_database(
    species='mouse', 
    signaling_type='Secreted Signaling', 
    database='CellChat'
)
print(df_cellchat.shape)

df_cellchat_filtered = ct.pp.filter_lr_database(
    df_cellchat, 
    adata_dis500, 
    min_cell_pct=0.05
)
print(df_cellchat_filtered.shape)
print(df_cellchat_filtered.head())

# Running takes about 30 minutes
ct.tl.spatial_communication(
    adata_dis500,
    database_name='cellchat',
    df_ligrec=df_cellchat_filtered,
    dis_thr=500,
    heteromeric=True,
    pathway_sum=True
)
adata_dis500.write("last_adata.commot.h5ad")

ct.tl.communication_direction(
    adata_dis500, 
    database_name='cellchat', 
    pathway_name='PSAP', 
    k=5
)

ct.pl.plot_cell_communication(
    adata_dis500,
    database_name='cellchat',
    pathway_name='PSAP',
    plot_method='grid',
    background_legend=True,
    scale=0.00003,
    ndsize=8,
    grid_density=0.4,
    summary='sender',
    background='image',
    clustering='leiden',
    cmap='Alphabet',
    normalize_v=True,
    normalize_v_quantile=0.995
)

adata_dis500.obs['leiden'] = adata.obs['leiden']

ct.tl.cluster_communication(
    adata_dis500, 
    database_name='cellchat', 
    pathway_name='PSAP', 
    clustering='leiden',
    n_permutations=100
)

# Install pygraphviz
conda install --channel conda-forge pygraphviz

ct.pl.plot_cluster_communication_network(
    adata_dis500,
    uns_names=['commot_cluster-leiden-cellchat-PSAP'],
    nx_node_pos=None,
    nx_bg_pos=False,
    p_value_cutoff=5e-2,
    filename='PSAP_cluster.pdf',
    nx_node_cmap='Light24'
)

ct.tl.cluster_position(adata_dis500, clustering='leiden')

ct.pl.plot_cluster_communication_network(
    adata_dis500,
    uns_names=['commot_cluster-leiden-cellchat-PSAP'],
    clustering='leiden',
    nx_node_pos='cluster',
    nx_pos_idx=np.array([0, 1]),
    nx_bg_pos=True,
    nx_bg_ndsize=0.25,
    p_value_cutoff=5e-2,
    filename='PSAP_cluster_spatial.pdf',
    nx_node_cmap='Light24'
)

adata_dis500 = sc.read_h5ad("./9.27_adata.commot.h5ad")
adata = sc.datasets.visium_sge(sample_id='V1_Mouse_Brain_Sagittal_Posterior')
adata_dis500.layers['counts'] = adata.X

df_deg, df_yhat = ct.tl.communication_deg_detection(
    adata_dis500,
    database_name='cellchat',
    pathway_name='PSAP',
    summary='receiver'
)

df_deg, df_yhat = ct.tl.communication_deg_detection(
    adata_dis500,
    database_name='cellchat',
    pathway_name='PSAP',
    summary='receiver'
)

# For R (BioC_mirror and tradeSeq installation)
# options(BioC_mirror="https://mirrors.westlake.edu.cn/bioconductor")
# BiocManager::install("tradeSeq")
# library(tradeSeq)

df_deg, df_yhat = ct.tl.communication_deg_detection(
    adata_dis500,
    database_name='cellchat',
    pathway_name='PSAP',
    summary='receiver'
)

# Install additional R packages
# BiocManager::install("tradeSeq")
# library(tradeSeq)
# BiocManager::install("clusterExperiment")
# library(clusterExperiment)
# install.packages("howmany_0.3-1.zip", repos=None, type="win.binary")
# library(howmany)

df_deg, df_yhat = ct.tl.communication_deg_detection(
    adata_dis500,
    database_name='cellchat',
    pathway_name='PSAP',
    summary='receiver'
)

# Save DEG results to pickle file
import pickle

deg_result = {"df_deg": df_deg, "df_yhat": df_yhat}

with open('./deg_PSAP.pkl', 'wb') as handle:
    pickle.dump(deg_result, handle, protocol=pickle.HIGHEST_PROTOCOL)

# Load DEG results from pickle file
with open("./deg_PSAP.pkl", 'rb') as file:
    deg_result = pickle.load(file)

df_deg_clus, df_yhat_clus = ct.tl.communication_deg_clustering(
    df_deg, 
    df_yhat, 
    deg_clustering_res=0.4
)

top_de_genes_PSAP = ct.pl.plot_communication_dependent_genes(
    df_deg_clus,
    df_yhat_clus,
    top_ngene_per_cluster=5,
    filename='./heatmap_deg_PSAP.pdf',
    font_scale=1.2,
    return_genes=True
)

X_sc = adata_dis500.obsm['spatial']
fig, ax = plt.subplots(1, 3, figsize=(15, 4))

colors = adata_dis500.obsm['commot-cellchat-sum-receiver']['r-PSAP'].values
idx = np.argsort(colors)
ax[0].scatter(X_sc[idx, 0], X_sc[idx, 1], c=colors[idx], cmap='coolwarm', s=10)

colors = adata_dis500[:, 'Ctxn1'].X.toarray().flatten()
idx = np.argsort(colors)
ax[1].scatter(X_sc[idx, 0], X_sc[idx, 1], c=colors[idx], cmap='coolwarm', s=10)

colors = adata_dis500[:, 'Gpr37'].X.toarray().flatten()
idx = np.argsort(colors)
ax[2].scatter(X_sc[idx, 0], X_sc[idx, 1], c=colors[idx], cmap='coolwarm', s=10)

ax[0].set_title('Amount of received signal')
ax[1].set_title('An example negative DE gene (Ctxn1)')
ax[2].set_title('An example positive DE gene (Gpr37)')

df_impact_PSAP = ct.tl.communication_impact(
    adata_dis500,
    database_name='cellchat',
    pathway_name='PSAP',
    tree_combined=True,
    method='treebased_score',
    tree_ntrees=100,
    tree_repeat=100,
    tree_method='rf',
    ds_genes=top_de_genes_PSAP,
    bg_genes=500,
    normalize=True
)

ct.pl.plot_communication_impact(
    df_impact_PSAP,
    summary='receiver',
    top_ngene=30,
    top_ncomm=5,
    colormap='coolwarm',
    font_scale=1.2,
    linewidth=0,
    show_gene_names=True,
    show_comm_names=True,
    cluster_knn=2,
    filename='heatmap_impact_PSAP.pdf'
)

df_impact_PSAP

```
