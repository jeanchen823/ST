# MERFISH



```bash
import scanpy as sc
import pandas as pd
import seaborn as sns
import squidpy as sq
import os
import sys
import anndata
import numpy as np
from scipy.sparse import csr_matrix
import matplotlib.pylab as plt

sc.settings.verbosity = 1  # verbosity errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
sc.settings.set_figure_params(
    dpi=80,
    frameon=False,
    figsize=(3, 3),
    facecolor='white'
)

print(os.getcwd())


adata = sq.read.vizgen(
    path="MERFISH\\",   counts_file="datasets_mouse_brain_map_BrainReceptorShowcase_Slice1_Replicate1_cell_by_gene_S1R1.csv",    meta_file="datasets_mouse_brain_map_BrainReceptorShowcase_Slice1_Replicate1_cell_metadata_S1R1.csv",
    library_id="spatial"
)
adata


adata = anndata.read_text(
    "MERFISH\\datasets_mouse_brain_map_BrainReceptorShowcase_Slice1_Replicate1_cell_by_gene_S1R1.csv",
    delimiter=",",
    first_column_names=True
)
cell_meta = pd.read_csv(
    "MERFISH\\datasets_mouse_brain_map_BrainReceptorShowcase_Slice1_Replicate1_cell_metadata_S1R1.csv",
    index_col=0
)

gene_s = ["Blank" in v for v in adata.var_names]
from collections import Counter
Counter(gene_s)
Counter(cell_meta["fov"])

blank_genes = np.array(["Blank" in v for v in adata.var_names])
adata.obsm["blank_genes"] = pd.DataFrame(
    adata[:, blank_genes].X.copy(),
    columns=adata.var_names[blank_genes],
    index=adata.obs_names
)

bool_arr = np.array([True, False, True, False])
inverted_bool_arr = ~bool_arr
inverted_bool_arr
adata = adata[:, ~blank_genes].copy()

# Convert original matrix to sparse matrix
adata.X = csr_matrix(adata.X)
cell_meta.set_index(cell_meta.index.astype("str"), inplace=True)

adata.obs = pd.merge(
    adata.obs,
    cell_meta,
    how="left",
    left_index=True,
    right_index=True
)
adata.obsm['spatial'] = adata.obs[["center_x", "center_y"]].values
adata.obs.drop(columns=["center_x", "center_y"], inplace=True)
adata

sc.pl.embedding(adata, basis="spatial", projection="2d")
sc.pp.calculate_qc_metrics(adata, percent_top=(50, 100, 200, 300), inplace=True)
adata
adata.obsm["blank_genes"].to_numpy().sum() / adata.var["total_counts"].sum() * 100
adata.obs.describe()

fig, axs = plt.subplots(1, 3, figsize=(15, 4))
axs[0].set_title("Total transcripts per cell")
axs[0].set_xlim(0, 2000)
axs[0].set_ylim(0, 2000)
sns.histplot(adata.obs["total_counts"], kde=False, ax=axs[0])
axs[1].set_title("Unique transcripts per cell")
axs[1].set_xlim(0, 450)
axs[1].set_ylim(0, 10000)
sns.histplot(adata.obs["n_genes_by_counts"], kde=False, ax=axs[1])
axs[2].set_title("Volume of segmented cells")
axs[2].set_xlim(5, 10000)
axs[1].set_ylim(0, 1000)
sns.histplot(adata.obs["volume"], kde=False, ax=axs[2])

sc.pp.filter_cells(adata, min_counts=10)  # Filter cells with fewer than 10 transcripts
sc.pp.filter_genes(adata, min_cells=10)

adata.layers["counts"] = adata.X.copy()
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=4000)
sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)
sc.pp.pca(adata)
sc.pp.neighbors(adata, n_pcs=20)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5, key_added='leiden_05')
sc.tl.leiden(adata, resolution=1)


sc.pl.umap(adata, wspace=0.4, color=["total_counts", "n_genes_by_counts", 'leiden_05', 'leiden'])

fig, ax = plt.subplots(figsize=(8, 4), dpi=200)
sc.pl.umap(adata, color="leiden", use_raw=False, legend_loc='on data', legend_fontsize=8, ax=ax)

fig, ax = plt.subplots(figsize=(8, 4), dpi=150)
sq.pl.spatial_scatter(
    adata,
    shape=None,
    color=['leiden'],
    size=0.1,
    wspace=0.4,
    dpi=100,
    ax=ax
)
fig, ax = plt.subplots(figsize=(8, 4), dpi=200)
sc.pl.umap(adata,color="leiden",use_raw=False,
		legend_loc='on data',legend_fontsize=8,ax=ax)
		

sc.tl.rank_genes_groups(adata, groupby='leiden', use_raw=False, method="wilcoxon")
sc.pl.rank_genes_groups_matrixplot(adata, groupby='leiden', values_to_plot='logfoldchanges', cmap='bwr')
adata.uns.keys()

sc.tl.filter_rank_genes_groups(
    adata,
    groupby="leiden",
    min_in_group_fraction=0.8,
    max_out_group_fraction=0.2
)
adata.uns.keys()

del adata.uns['rank_genes_groups_filtered']

sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby="leiden",
    standard_scale="var",
    n_genes=5,
    cmap='bwr'
)

df_markers = sc.get.rank_genes_groups_df(adata, group=adata.obs['leiden'].unique())
df_markers = df_markers.loc[~df_markers.names.isna()]
df_markers['abs_score'] = df_markers.scores.abs()
df_markers.sort_values('abs_score', ascending=False, inplace=True)

top0_markers = df_markers.groupby('group').head(0).sort_values('group')['names'].unique().tolist()
top2_markers = df_markers.groupby('group').head(2).sort_values('group')['names'].unique().tolist()
top9_markers = df_markers.groupby('group').head(9).sort_values('group')['names'].unique().tolist()
df_markers

df9 = adata.obsm
sc.pl.umap(
    adata,
    color=['n_counts', 'volume'] + top2_markers[1:14],
    use_raw=False,
    color_map='jet',
    ncols=5
)

sq.pl.spatial_scatter(
    adata,
    spatial_key="spatial",
    shape=None,
    use_raw=False,
    color=['n_counts', 'volume'] + top2_markers[1:14],
    alpha=0.8,
    marker='o',
    size=0.5,
    na_color='white',
    dpi=100,
    colorbar=True,
    cmap='jet',
    ncols=5
)

for cluster in adata.obs.leiden.unique().tolist()[1:10]:
    sq.pl.spatial_scatter(
        adata,
        spatial_key="spatial",
        shape=None,
        color='leiden',
        library_id='spatial',
        groups=[cluster],
        title=cluster,
        alpha=0.8,
        marker='o',
        size=0.5,
        na_color='white',
        dpi=100,
        colorbar=True
    )

fig, ax = plt.subplots(figsize=(8, 4), dpi=150)
sq.pl.spatial_scatter(
    adata,
    spatial_key="spatial",
    shape=None,
    color='leiden',
    groups=['0', '2', '4', '7', '10', '11'],
    alpha=0.8,
    marker='o',
    size=0.5,
    na_color='white',
    ax=ax,
    dpi=200,
    colorbar=True
)

sq.gr.spatial_neighbors(adata, coord_type="generic", delaunay=True)
adata
sq.gr.centrality_scores(adata, cluster_key="leiden")


import copy

sq.gr.centrality_scores(adata, "leiden")

# sc.set_figure_params(figsize=(20, 8))

# copy centrality data to new DataFrame
df_central = copy.deepcopy(adata.uns["leiden_centrality_scores"])

ser_counts = adata.obs["leiden"].value_counts()
ser_counts.name = "cell counts"
meta_leiden = pd.DataFrame(ser_counts)
meta_leiden

df_central.index = meta_leiden.index.tolist()  # sort clusters based on centrality scores

#################################################

# closeness centrality - measure of how close the group is to other nodes.
ser_closeness = df_central["closeness_centrality"].sort_values(ascending=False)

# degree centrality - fraction of non-group members connected to group members.
# [Networkx](https://networkx.org/documentation/stable/reference/algorithms/generated/networkx.algorithms.centrality.degree_centrality.html#networkx.algorithms.centrality.degree_centrality)
# The degree centrality for a node v is the fraction of nodes it is connected to.
ser_degree = df_central["degree_centrality"].sort_values(ascending=False)

# clustering coefficient - measure of the degree to which nodes cluster together.
ser_cluster = df_central["average_clustering"].sort_values(ascending=False)

sq.pl.centrality_scores(adata, cluster_key="leiden", figsize=(16, 6))

inst_clusters = ser_closeness.index.tolist()[:5]
print(inst_clusters)

sq.pl.spatial_scatter(
    adata, groups=inst_clusters, color="leiden", size=15, 
    img=False, figsize=(10, 10), shape=None
)


inst_clusters = ser_closeness.index.tolist()[-5:]
print(inst_clusters)

sq.pl.spatial_scatter(
    adata, groups=inst_clusters, color="leiden", size=15, 
    img=False, figsize=(10, 10), shape=None
)

inst_clusters = ser_degree.index.tolist()[:5]
print(inst_clusters)

sq.pl.spatial_scatter(
    adata, groups=inst_clusters, color="leiden", size=15, 
    img=False, figsize=(10, 10), shape=None
)

inst_clusters = ser_degree.index.tolist()[-5:]
print(inst_clusters)

sq.pl.spatial_scatter(
    adata, groups=inst_clusters, color="leiden", size=15, 
    img=False, figsize=(10, 10), shape=None
)

inst_clusters = ser_cluster.index.tolist()[:5]
print(inst_clusters)

sq.pl.spatial_scatter(
    adata, groups=inst_clusters, color="leiden", size=15, 
    img=False, figsize=(10, 10), shape=None
)

inst_clusters = ser_cluster.index.tolist()[-5:]
print(inst_clusters)

sq.pl.spatial_scatter(
    adata, groups=inst_clusters, color="leiden", size=15, 
    img=False, figsize=(10, 10), shape=None
)

# 为了加快运行速度，进行抽样
adata_subsample = sc.pp.subsample(adata, fraction=0.3, copy=True)
# adata_subsample = sc.pp.subsample(adata, n_obs=2000, copy=True)

## 开始计算
sq.gr.co_occurrence(
    adata_subsample,
    cluster_key="leiden",
)

sq.pl.co_occurrence(
    adata_subsample,
    cluster_key="leiden",
    clusters="12",
    figsize=(15, 10),
)

sq.pl.spatial_scatter(
    adata_subsample,
    color="leiden",
    shape=None,
    size=2,
    figsize=(15, 10),
)

sq.gr.nhood_enrichment(adata, cluster_key="leiden")


fig, ax = plt.subplots(1, 2, figsize=(18, 9))

sq.pl.nhood_enrichment(
    adata,
    cluster_key="leiden",
    figsize=(9, 9),
    title="Neighborhood enrichment adata",
    ax=ax[0],
)

sq.pl.spatial_scatter(adata_subsample, color="leiden", shape=None, size=2, ax=ax[1])


mode = "L"
sq.gr.ripley(adata, cluster_key="leiden", mode=mode)

fig, ax = plt.subplots(1, 2, figsize=(15, 7))

sq.pl.ripley(adata, cluster_key="leiden", mode=mode, ax=ax[0])

sq.pl.spatial_scatter(
    adata_subsample,
    color="leiden",
    groups=["0", "1", "3"],
    shape=None,
    size=2,
    ax=ax[1],
)

# SPARK
# Paper: https://www.nature.com/articles/s41592-019-0701-7#Abs1
# Code: https://github.com/xzhoulab/SPARK

# SpatialDE
# Paper: https://www.nature.com/articles/nmeth.4636
# Code: https://github.com/Teichlab/SpatialDE

# trendsceek
# Paper: https://www.nature.com/articles/nmeth.4634
# Code: https://github.com/edsgard/trendsceek

# HMRF
# Paper: https://www.nature.com/articles/nbt.4260
# Code: https://bitbucket.org/qzhudfci/smfishhmrf-py



# 因为用了抽样的数据  所以重新进行邻域构建
sq.gr.spatial_neighbors(adata_subsample, coord_type="generic", delaunay=True)

sq.gr.spatial_autocorr(
    adata_subsample,
    mode="moran",
    n_perms=1sq.pl.spatial_scatter(
    adata_subsample, color=bot_autocorr, size=20, cmap="Reds", 
    img=False, figsize=(5, 5), shape=None
)

adata_subsample.write_h5ad('adata_subsample_squidpy.h5ad')
adata.write_h5ad('adata_squidpy.h5ad')

del adata.uns['rank_genes_groups_filtered']
adata.write_h5ad('adata_squidpy.h5ad')

del adata_subsample.uns['rank_genes_groups_filtered']
adata_subsample.write_h5ad('adata_subsample_squidpy.h5ad')



sq.gr.ligrec(
    adata_subsample, use_raw=False, 
    threshold=0, n_perms=100, 
    cluster_key="leiden", 
    clusters=["0", "1", "3"],
)

df1=adata_subsample.uns['leiden_ligrec']


sq.pl.ligrec(
    adata_subsample,
    cluster_key="leiden",
    source_groups="0",
    target_groups=["1", "3"],
    means_range=(0, np.inf),
    alpha=1e-4,
    swap_axes=True,
)


sq.pl.ligrec(
    adata_subsample,
    cluster_key="leiden",
    source_groups="0",
    target_groups=["1", "3"],
    means_range=(0, np.inf),
    alpha=1e-4,
    swap_axes=True,
)

00,
    n_jobs=1,
)

adata_subsample.uns["moranI"].head(10)


num_view = 12

top_autocorr = (
    adata.uns["moranI"]["I"].sort_values(ascending=False).head(num_view).index.tolist()
)

bot_autocorr = (
    adata.uns["moranI"]["I"].sort_values(ascending=True).head(num_view).index.tolist()
)

sq.pl.spatial_scatter(
    adata_subsample, color=top_autocorr, size=20, cmap="Reds", 
    img=False, figsize=(5, 5), shape=None
)



sq.pl.spatial_scatter(
    adata_subsample, color=bot_autocorr, size=20, cmap="Reds", 
    img=False, figsize=(5, 5), shape=None
)

adata_subsample.write_h5ad('adata_subsample_squidpy.h5ad')
adata.write_h5ad('adata_squidpy.h5ad')

del adata.uns['rank_genes_groups_filtered']
adata.write_h5ad('adata_squidpy.h5ad')

del adata_subsample.uns['rank_genes_groups_filtered']
adata_subsample.write_h5ad('adata_subsample_squidpy.h5ad')



sq.gr.ligrec(
    adata_subsample, use_raw=False, 
    threshold=0, n_perms=100, 
    cluster_key="leiden", 
    clusters=["0", "1", "3"],
)

df1=adata_subsample.uns['leiden_ligrec']


sq.pl.ligrec(
    adata_subsample,
    cluster_key="leiden",
    source_groups="0",
    target_groups=["1", "3"],
    means_range=(0, np.inf),
    alpha=1e-4,
    swap_axes=True,
)


sq.pl.ligrec(
    adata_subsample,
    cluster_key="leiden",
    source_groups="0",
    target_groups=["1", "3"],
    means_range=(0, np.inf),
    alpha=1e-4,
    swap_axes=True,
)


```

