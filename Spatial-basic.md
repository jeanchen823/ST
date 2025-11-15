# Spatial-basic

```bash
https://scanpy.readthedocs.io/en/stable/tutorials/spatial/basic-analysis.html

import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sc.logging.print_versions()
sc.set_figure_params(facecolor="white", figsize=(8, 8))
sc.settings.verbosity = 3


https://support.10xgenomics.com/spatial-gene-expression/datasets/1.0.0/V1_Human_Lymph_Node?

# 读取数据
adata = sc.datasets.visium_sge(sample_id="V1_Human_Lymph_Node")
adata.var_names_make_unique()

# 计算线粒体基因组并保存为注释
adata.var["mt"] = adata.var_names.str.startswith("MT-")
adata.var["mt"]

# 计算质控指标，qc的var选择var.mt
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

# 可视化质控数据
fig, axs = plt.subplots(1, 4, figsize=(15, 4))
sns.distplot(adata.obs["total_counts"], kde=False, ax=axs[0])
sns.distplot(adata.obs["total_counts"][adata.obs["total_counts"] < 10000], kde=False, bins=40, ax=axs[1])
sns.distplot(adata.obs["n_genes_by_counts"], kde=False, bins=60, ax=axs[2])
sns.distplot(adata.obs["n_genes_by_counts"][adata.obs["n_genes_by_counts"] < 4000], kde=False, bins=60, ax=axs[3])

# 过滤低质量细胞

sc.pp.filter_cells(adata, min_counts=5000)
sc.pp.filter_cells(adata, max_counts=35000)
adata = adata[adata.obs["pct_counts_mt"] < 20]
print(f"#cells after MT filter: {adata.n_obs}")
sc.pp.filter_genes(adata, min_cells=10)


# 归一化处理
sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)

# PCA降维
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata, key_added="clusters")

# 可视化UMAP
adata.obsm['spatial'].shape
adata.obsm['spatial']

plt.rcParams["figure.figsize"] = (8, 8)
sc.pl.spatial(adata, img_key="hires", color=["total_counts", "n_genes_by_counts"])

sc.pl.spatial(adata, img_key="hires", color="clusters", size=1.5)
sc.pl.spatial(adata, img_key="hires", color="clusters", groups=["0", "4"], alpha=0.5, size=1.3)


sc.tl.rank_genes_groups(adata, "clusters", method="t-test")
sc.pl.rank_genes_groups_heatmap(adata, groups="4", n_genes=10, groupby="clusters")
sc.pl.rank_genes_groups_heatmap(adata, groups="2", n_genes=10, groupby="clusters")
sc.pl.rank_genes_groups_heatmap(adata, groups="0", n_genes=10, groupby="clusters")


sc.pl.spatial(adata, img_key="hires", color=["clusters", "CR2"])


os.getcwd()
os.chdir('F:/spe.lesson\\k8\\')
adata.write("8_25.h5ad",compression='gzip')


counts = pd.DataFrame(adata.X.todense(), columns=adata.var_names, index=adata.obs_names)
coord = pd.DataFrame(adata.obsm['spatial'], columns=['x_coord', 'y_coord'], index=adata.obs_names)
coord


pip install spatialde

https://github.com/Teichlab/SpatialDE/issues/37

```



```bash
https://github.com/drieslab/spatial-datasets/tree/master/data/2019_merFISH_PNAS

coordinates = pd.read_excel("./pnas.1912459116.sd15.xlsx", index_col=0)
counts = sc.read_csv("./pnas.1912459116.sd12.csv").transpose()

adata_merfish = counts[coordinates.index, :]
adata_merfish.obsm["spatial"] = coordinates.to_numpy()

adata_merfish
pip  install openpyxl

coordinates = pd.read_excel("pnas.1912459116.sd15.xlsx", index_col=0)
counts = sc.read_csv("pnas.1912459116.sd12.csv").transpose()
adata_merfish = counts[coordinates.index, :]
adata_merfish.obsm["spatial"] = coordinates.to_numpy()
adata_merfish


sc.pp.normalize_per_cell(adata_merfish, counts_per_cell_after=1e6)
sc.pp.log1p(adata_merfish)
sc.pp.pca(adata_merfish, n_comps=15)
sc.pp.neighbors(adata_merfish)
sc.tl.umap(adata_merfish)
sc.tl.leiden(adata_merfish, key_added="clusters", resolution=0.5)


sc.pl.umap(adata_merfish, color="clusters")
sc.pl.embedding(adata_merfish, basis="spatial", color="clusters")


```



```bash
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
sc.logging.print_versions()
sc.set_figure_params(facecolor="white", figsize=(8, 8))
sc.settings.verbosity = 3

# 获取当前工作目录

print(os.getcwd())

print(os.path )

adata = sc.read_visium('GBM1_spaceranger_out\\')
adata.var_names_make_unique()

adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
adata

fig, axs = plt.subplots(1, 4, figsize=(15, 4))
sns.histplot(adata.obs["total_counts"], kde=False, ax=axs[0])
sns.histplot(
    adata.obs["total_counts"][adata.obs["total_counts"] < 10000],
    kde=False,bins=40,ax=axs[1],)
sns.histplot(adata.obs["n_genes_by_counts"], kde=False, bins=60, ax=axs[2])
sns.histplot(
    adata.obs["n_genes_by_counts"][adata.obs["n_genes_by_counts"] < 4000],
    kde=False,bins=60,ax=axs[3],)

sc.pp.filter_cells(adata, min_counts=5000)
sc.pp.filter_cells(adata, max_counts=35000)
adata = adata[adata.obs["pct_counts_mt"] < 20].copy()
print(f"#cells after MT filter: {adata.n_obs}")
sc.pp.filter_genes(adata, min_cells=10)

sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)

sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata, key_added="clusters" )

plt.rcParams["figure.figsize"] = (4, 4)
sc.pl.umap(adata, color=["total_counts", "n_genes_by_counts", "clusters"], wspace=0.4)
sc.pl.umap(adata, color=["total_counts" ], wspace=0.4)
sc.pl.umap(adata, color=[ "n_genes_by_counts"], wspace=0.4)
sc.pl.umap(adata, color=[ "clusters"], wspace=0.4)

plt.rcParams["figure.figsize"] = (8, 8)
sc.pl.spatial(adata, img_key="hires", color=["total_counts", "n_genes_by_counts"])


sc.pl.spatial(adata, img_key="hires", color="clusters", size=1.5)

sc.pl.spatial( adata,img_key="hires",color="clusters",
    groups=["1", "3"],crop_coord=[7000, 10000, 0, 6000],
    alpha=0.5,size=1.3,)


sc.tl.rank_genes_groups(adata, "clusters", method="t-test")
sc.pl.rank_genes_groups_heatmap(adata, groups="2", n_genes=10, groupby="clusters")

sc.pl.spatial(adata, img_key="hires", color=["clusters" ])

sc.pl.spatial(adata, img_key="hires", color=["COL1A2", "SYPL1"], alpha=0.7)


adata.write("GSE194329.data.anndata.h5ad")
```



```bash
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sc.logging.print_versions()
sc.set_figure_params(facecolor="white", figsize=(8, 8))
sc.settings.verbosity = 3

adata = sc.datasets.visium_sge(sample_id="V1_Human_Lymph_Node")
adata.var_names_make_unique()

adata


# 将线粒体基因组保存为注释 var.mt
adata.var["mt"] = adata.var_names.str.startswith("MT-")
adata.var["mt"]

#计算指标, qc的var选择 var.mt
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

adata

fig, axs = plt.subplots(1, 4, figsize=(15, 4))
sns.distplot(adata.obs["total_counts"], kde=False, ax=axs[0])
sns.distplot(adata.obs["total_counts"][adata.obs["total_counts"] < 10000], kde=False, bins=40, ax=axs[1])
sns.distplot(adata.obs["n_genes_by_counts"], kde=False, bins=60, ax=axs[2])
sns.distplot(adata.obs["n_genes_by_counts"][adata.obs["n_genes_by_counts"] < 4000], kde=False, bins=60, ax=axs[3])


sc.pp.filter_cells(adata, min_counts=5000)
sc.pp.filter_cells(adata, max_counts=35000)
adata = adata[adata.obs["pct_counts_mt"] < 20]
print(f"#cells after MT filter: {adata.n_obs}")
sc.pp.filter_genes(adata, min_cells=10)


sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)

sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata, key_added="clusters")

plt.rcParams["figure.figsize"] = (4, 4)
sc.pl.umap(adata, color=["total_counts", "n_genes_by_counts", "clusters"], wspace=0.4)

adata.obsm['spatial'].shape
adata.obsm['spatial']

plt.rcParams["figure.figsize"] = (8, 8)
sc.pl.spatial(adata, img_key="hires", color=["total_counts", "n_genes_by_counts"])
sc.pl.spatial(adata, img_key="hires", color="clusters", size=1.5)
sc.pl.spatial(adata, img_key="hires", color="clusters", groups=["0", "4"], alpha=0.5, size=1.3)


sc.tl.rank_genes_groups(adata, "clusters", method="t-test")
sc.pl.rank_genes_groups_heatmap(adata, groups="4", n_genes=10, groupby="clusters")
sc.pl.rank_genes_groups_heatmap(adata, groups="2", n_genes=10, groupby="clusters")
sc.pl.rank_genes_groups_heatmap(adata, groups="0", n_genes=10, groupby="clusters")

sc.pl.spatial(adata, img_key="hires", color=["clusters", "CR2"])

import os
os.getcwd()
os.chdir('F:/spe.lesson\\k8\\')
adata.write("8_25.h5ad",compression='gzip')


counts = pd.DataFrame(adata.X.todense(), columns=adata.var_names, index=adata.obs_names)
coord = pd.DataFrame(adata.obsm['spatial'], columns=['x_coord', 'y_coord'], index=adata.obs_names)

counts

coord

import SpatialDE
results = SpatialDE.run(coord, counts)
results = SpatialDE.run(np.array(coord), counts)




import urllib.request





coordinates = pd.read_excel("pnas.1912459116.sd15.xlsx", index_col=0)
counts = sc.read_csv("pnas.1912459116.sd12.csv").transpose()
adata_merfish = counts[coordinates.index, :]
adata_merfish.obsm["spatial"] = coordinates.to_numpy()
adata_merfish

sc.pp.normalize_per_cell(adata_merfish, counts_per_cell_after=1e6)
sc.pp.log1p(adata_merfish)
sc.pp.pca(adata_merfish, n_comps=15)
sc.pp.neighbors(adata_merfish)
sc.tl.umap(adata_merfish)
sc.tl.leiden(adata_merfish, key_added="clusters", resolution=0.5)


sc.pl.umap(adata_merfish, color="clusters")
sc.pl.embedding(adata_merfish, basis="spatial", color="clusters")



##########################################
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
sc.logging.print_versions()
sc.set_figure_params(facecolor="white", figsize=(8, 8))
sc.settings.verbosity = 3

# 获取当前工作目录

print(os.getcwd())


print(os.path )

adata = sc.read_visium('GBM1_spaceranger_out\\')
adata.var_names_make_unique()

adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
adata

fig, axs = plt.subplots(1, 4, figsize=(15, 4))
sns.histplot(adata.obs["total_counts"], kde=False, ax=axs[0])
sns.histplot(
    adata.obs["total_counts"][adata.obs["total_counts"] < 10000],
    kde=False,bins=40,ax=axs[1],)
sns.histplot(adata.obs["n_genes_by_counts"], kde=False, bins=60, ax=axs[2])
sns.histplot(
    adata.obs["n_genes_by_counts"][adata.obs["n_genes_by_counts"] < 4000],
    kde=False,bins=60,ax=axs[3],)

sc.pp.filter_cells(adata, min_counts=5000)
sc.pp.filter_cells(adata, max_counts=35000)
adata = adata[adata.obs["pct_counts_mt"] < 20].copy()
print(f"#cells after MT filter: {adata.n_obs}")
sc.pp.filter_genes(adata, min_cells=10)

sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)

sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata, key_added="clusters" )

plt.rcParams["figure.figsize"] = (4, 4)
sc.pl.umap(adata, color=["total_counts", "n_genes_by_counts", "clusters"], wspace=0.4)
sc.pl.umap(adata, color=["total_counts" ], wspace=0.4)
sc.pl.umap(adata, color=[ "n_genes_by_counts"], wspace=0.4)
sc.pl.umap(adata, color=[ "clusters"], wspace=0.4)

plt.rcParams["figure.figsize"] = (8, 8)
sc.pl.spatial(adata, img_key="hires", color=["total_counts", "n_genes_by_counts"])


sc.pl.spatial(adata, img_key="hires", color="clusters", size=1.5)

sc.pl.spatial( adata,img_key="hires",color="clusters",
    groups=["1", "3"],crop_coord=[7000, 10000, 0, 6000],
    alpha=0.5,size=1.3,)


sc.tl.rank_genes_groups(adata, "clusters", method="t-test")
sc.pl.rank_genes_groups_heatmap(adata, groups="2", n_genes=10, groupby="clusters")

sc.pl.spatial(adata, img_key="hires", color=["clusters" ])

sc.pl.spatial(adata, img_key="hires", color=["COL1A2", "SYPL1"], alpha=0.7)

adata.write("GSE194329.data.anndata.h5ad")
```

