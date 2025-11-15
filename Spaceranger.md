# Spaceranger

```python
spaceranger count --id=8_7_k2 \
 --transcriptome=/mnt/f/spe.lesson/k2.spaceranger/ref/hg38/refdata-gex-GRCh38-2020-A \
 --probe-set=Visium_FFPE_Human_Prostate_Acinar_Cell_Carcinoma_probe_set.csv \
 --fastqs=Visium_FFPE_Human_Prostate_Acinar_Cell_Carcinoma_fastqs/ \
 --sample=Prostate \
 --image=Visium_FFPE_Human_Prostate_Acinar_Cell_Carcinoma_image.tif \
 --localcores=16 \
 --slide=V11J26-002 \
  --area=B1     --create-bam  true \
  --localmem=128



spaceranger count --id=huagek2 \
 --transcriptome=/huage/spatial/k2/ref/refdata-gex-GRCh38-2020-A \
 --probe-set=Visium_FFPE_Human_Prostate_Acinar_Cell_Carcinoma_probe_set.csv \
 --fastqs=Visium_FFPE_Human_Prostate_Acinar_Cell_Carcinoma_fastqs/ \
 --sample=Prostate \
 --image=Visium_FFPE_Human_Prostate_Acinar_Cell_Carcinoma_image.tif \
 --localcores=16 \
 --slide=V11J26-002 \
  --area=B1     --create-bam  true \
  --localmem=128

```



```python
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 18:58:37 2024

@author: admin
"""

import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

sc.logging.print_versions()
sc.set_figure_params(facecolor="white", figsize=(8, 8))
sc.settings.verbosity = 3

print(os.getcwd())

os.chdir('GSE194329_RAW\\')

print(os.getcwd())

adata = sc.read_visium('GBM1_spaceranger_out\\')

adata 


adata.var_names_make_unique()

adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)


adata

adata1=adata

adata1=adata.copy()

fig, axs = plt.subplots(1, 4, figsize=(15, 4))
sns.histplot(adata.obs["total_counts"], kde=False, ax=axs[0])
sns.histplot(adata.obs["total_counts"][adata.obs["total_counts"] < 10000],
    kde=False, bins=40,ax=axs[1],)
sns.histplot(adata.obs["n_genes_by_counts"], kde=False, bins=60, ax=axs[2])
sns.histplot(adata.obs["n_genes_by_counts"][adata.obs["n_genes_by_counts"] < 4000],
    kde=False, bins=60, ax=axs[3],)




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


sc.pl.spatial(adata,
    img_key="hires",color="clusters",
    groups=["1", "3"],crop_coord=[7000, 10000, 0, 6000],
    alpha=0.5,size=1.3,)


sc.tl.rank_genes_groups(adata, "clusters", method="t-test")
sc.pl.rank_genes_groups_heatmap(adata, groups="2", n_genes=10, groupby="clusters")


sc.pl.spatial(adata, img_key="hires", color=["clusters" ])

sc.pl.spatial(adata, img_key="hires", color=["COL1A2", "SYPL1"], alpha=0.7)


adata.write("GSE194329.data.anndata.h5ad")




```





