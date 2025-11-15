# Xenium

```bash
https://www.cell.com/cell/fulltext/S0092-8674(24)00233-2

2024年 Cell：Cellular architecture of evolving neuroinflammatory lesions and multiple sclerosis pathology

2021年 Nature：Molecular architecture of the developing mouse brain

NC：High resolution mapping of the tumor microenvironment using integrated single-cell, spatial and in situ analysis

conda create -n squidpy_xenium  python=3.10.14  -c https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/
conda-forgeconda activate squidpy_xenium

conda install spyder-kernels=2.5  -c https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge
conda install spyder-notebook  -c https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge

pip install squidpy -i https://mirrors.aliyun.com/pypi/simple/


pip install  spatialdata_io -i https://mirrors.aliyun.com/pypi/simple/
pip install  spatialdata_plot  -i https://mirrors.aliyun.com/pypi/simple/
pip install napari_spatialdata  -i https://mirrors.aliyun.com/pypi/simple/

https://www.nature.com/articles/s41467-023-43458-x

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243168



import spatialdata as sd
from spatialdata_io import xenium

import matplotlib.pyplot as plt
import seaborn as sns

import scanpy as sc
import squidpy as sq
import os
plt.rcParams['pdf.fonttype'] = 42


print(os.getcwd())

#定义存储路径
xenium_path = "GSM7780153_Xenium_FFPE_Human_Breast_Cancer_Rep1_outs\\outs\\"

sdata = xenium(xenium_path)
sdata

SpatialData object
├── Images
│     ├── 'morphology_focus': DataTree[cyx] (1, 25778, 35416), (1, 12889, 17708), (1, 6444, 8854), (1, 3222, 4427), (1, 1611, 2213)
│     └── 'morphology_mip': DataTree[cyx] (1, 25778, 35416), (1, 12889, 17708), (1, 6444, 8854), (1, 3222, 4427), (1, 1611, 2213)
├── Labels
│     ├── 'cell_labels': DataTree[yx] (25778, 35416), (12889, 17708), (6444, 8854), (3222, 4427), (1611, 2213)
│     └── 'nucleus_labels': DataTree[yx] (25778, 35416), (12889, 17708), (6444, 8854), (3222, 4427), (1611, 2213)
├── Points
│     └── 'transcripts': DataFrame with shape: (<Delayed>, 8) (3D points)
├── Shapes
│     ├── 'cell_boundaries': GeoDataFrame shape: (167780, 1) (2D shapes)
│     ├── 'cell_circles': GeoDataFrame shape: (167780, 2) (2D shapes)
│     └── 'nucleus_boundaries': GeoDataFrame shape: (167780, 1) (2D shapes)
└── Tables
      └── 'table': AnnData (167780, 313)
with coordinate systems:
    ▸ 'global', with elements:
        morphology_focus (Images), morphology_mip (Images), cell_labels (Labels), nucleus_labels (Labels), transcripts (Points), cell_boundaries (Shapes), cell_circles (Shapes), nucleus_boundaries (Shapes)
        

adata = sdata.tables["table"]
adata

adata.obs

adata.obsm["spatial"]

sc.pp.calculate_qc_metrics(adata, percent_top=(10, 20, 50, 150), inplace=True)


cprobes = (
    adata.obs["control_probe_counts"].sum() / adata.obs["total_counts"].sum() * 100)
cwords = (
    adata.obs["control_codeword_counts"].sum() / adata.obs["total_counts"].sum() * 100)
print(f"Negative DNA probe count % : {cprobes}")
print(f"Negative decoding count % : {cwords}")



fig, axs = plt.subplots(1, 4, figsize=(15, 4))
axs[0].set_title("Total transcripts per cell")
sns.histplot(
    adata.obs["total_counts"],
    kde=False,ax=axs[0],)

axs[1].set_title("Unique transcripts per cell")
sns.histplot(
    adata.obs["n_genes_by_counts"],
    kde=False,ax=axs[1],)

axs[2].set_title("Area of segmented cells")
sns.histplot(
    adata.obs["cell_area"],
    kde=False,ax=axs[2],)

axs[3].set_title("Nucleus ratio")
sns.histplot(
    adata.obs["nucleus_area"] / adata.obs["cell_area"],
    kde=False,ax=axs[3],)


sc.pp.filter_cells(adata, min_counts=10)
sc.pp.filter_genes(adata, min_cells=5)

adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata)

sc.pl.umap(adata,
    color=["total_counts",
        "n_genes_by_counts",
        "leiden",],
    wspace=0.4,)


sq.pl.spatial_scatter(adata,
    library_id="spatial",
    shape=None,
    color=["leiden",],
    wspace=0.4,)


sq.gr.spatial_neighbors(adata, coord_type="generic", delaunay=True)

sq.gr.centrality_scores(adata, cluster_key="leiden")

sq.pl.centrality_scores(adata, cluster_key="leiden", figsize=(16, 5))

sdata.tables["subsample"] = sc.pp.subsample(adata, fraction=0.5, copy=True)

adata_subsample = sdata.tables["subsample"]

sq.gr.co_occurrence(
    adata_subsample,
    cluster_key="leiden",)
sq.pl.co_occurrence(
    adata_subsample,
    cluster_key="leiden",
    clusters="12",
    figsize=(10, 10),)
sq.pl.spatial_scatter(
    adata_subsample,
    color="leiden",
    shape=None,
    size=2,)


sq.gr.nhood_enrichment(adata, cluster_key="leiden")
fig, ax = plt.subplots(1, 2, figsize=(13, 7))
sq.pl.nhood_enrichment(
    adata,
    cluster_key="leiden",
    figsize=(8, 8),
    title="Neighborhood enrichment adata",
    ax=ax[0],
)
sq.pl.spatial_scatter(adata_subsample, color="leiden", shape=None, size=2, ax=ax[1])


sq.gr.spatial_neighbors(adata_subsample, coord_type="generic", delaunay=True)
sq.gr.spatial_autocorr(
    adata_subsample,
    mode="moran",
    n_perms=100,
    n_jobs=1,)
adata_subsample.uns["moranI"].head(10)
df=adata_subsample.uns["moranI"]
sq.pl.spatial_scatter(
    adata_subsample,
    library_id="spatial",
    color= ["KRT7", "EPCAM"],
    shape=None,
    size=2,
    img=False,)

import spatialdata_plot

gene_name = ["KRT7", "EPCAM"]
for name in gene_name:
    sdata.pl.render_images("morphology_focus").pl.render_shapes(
        "cell_circles",
        color=name,
        table_name="table",
        use_raw=False,
    ).pl.show(
        title=f"{name} expression over Morphology image",
        coordinate_systems="global",
        figsize=(10, 5),)

from napari_spatialdata import Interactive

Interactive(sdata)
```


