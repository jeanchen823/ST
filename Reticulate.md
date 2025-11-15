# Reticulate

```bash
2. 安装reticulate from CRAN
install.packages(“reticulate”)
## 加载R包
library(reticulate)
3. 配置python环境
方法一：直接指定使用python版本的执行程序

# 直接指定使用python版本的执行程序
use_python("D:/miniconda/envs/spyder/")
或者使用 Sys.setenv
Sys.setenv(RETICULATE_PYTHON = ".")

7.1 加载包
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

7.2 图片设置
sc.logging.print_versions()
sc.set_figure_params(facecolor="white", figsize=(8, 8))
sc.settings.verbosity = 3

7.3 设置路径
# 获取当前工作目录
print(os.getcwd())
os.chdir('GSE194329_RAW\\')
print(os.path )

7.4 读取数据
adata = sc.read_visium('GSE194329_RAW\\GBM1_spaceranger_out\\')



7.5 质控
adata.var_names_make_unique()
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

7.6 查看对象
adata


7.7 查看质控结果
fig, axs = plt.subplots(1, 4, figsize=(15, 4))
sns.histplot(adata.obs["total_counts"], kde=False, ax=axs[0])
sns.histplot(adata.obs["total_counts"][adata.obs["total_counts"] < 10000],
  kde=False, bins=40,ax=axs[1],)

sns.histplot(adata.obs["n_genes_by_counts"], kde=False, bins=60, ax=axs[2])
sns.histplot( adata.obs["n_genes_by_counts"][adata.obs["n_genes_by_counts"] < 4000],
  kde=False,bins=60,ax=axs[3],)
  
  
7.8  过滤低质量的细胞

sc.pp.filter_cells(adata, min_counts=5000)
sc.pp.filter_cells(adata, max_counts=35000)
adata = adata[adata.obs["pct_counts_mt"] < 20].copy()
print(f"#cells after MT filter: {adata.n_obs}")
sc.pp.filter_genes(adata, min_cells=10)

7.9 数据归一化和找高变基因
sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)


7.10 降维聚类
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata, key_added="clusters" )


7.11 画UMAP图 展示
plt.rcParams["figure.figsize"] = (20, 20)
sc.pl.umap(adata, color=["total_counts", "n_genes_by_counts", "clusters"], wspace=0.4)
sc.pl.umap(adata, color=["total_counts" ], wspace=0.4)
sc.pl.umap(adata, color=[ "n_genes_by_counts"], wspace=0.4)
sc.pl.umap(adata, color=[ "clusters"], wspace=0.4)


7.12 空间结果展示
plt.rcParams["figure.figsize"] = (8, 8)
sc.pl.spatial(adata, img_key="hires", color=["total_counts", "n_genes_by_counts"])
sc.pl.spatial(adata, img_key="hires", color="clusters", size=1.5)

sc.pl.spatial(adata,img_key="hires", color="clusters",
  groups=["1", "3"],crop_coord=[7000, 10000, 0, 6000],
  alpha=0.5,size=1.3,)
7.13 找marker基因

sc.tl.rank_genes_groups(adata, "clusters", method="t-test")
sc.pl.rank_genes_groups_heatmap(adata, groups="2", n_genes=10, groupby="clusters")


7.14  特定基因展示

sc.pl.spatial(adata, img_key="hires", color=["clusters" ])

sc.pl.spatial(adata, img_key="hires", color=["COL1A2", "SYPL1"], alpha=0.7)

7.15 保存数据

adata.write("GSE194329.data.anndata.h5ad")




```



```bash
library(reticulate)

use_python("D:/miniconda/envs/spyder/") 

py_config()  
py_available()
repl_python()

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
sns.histplot(adata.obs["total_counts"][adata.obs["total_counts"] < 10000],
             kde=False, bins=40,ax=axs[1],)

sns.histplot(adata.obs["n_genes_by_counts"], kde=False, bins=60, ax=axs[2])
sns.histplot( adata.obs["n_genes_by_counts"][adata.obs["n_genes_by_counts"] < 4000],
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


plt.rcParams["figure.figsize"] = (20, 20)
sc.pl.umap(adata, color=["total_counts", "n_genes_by_counts", "clusters"], wspace=0.4)

sc.pl.umap(adata, color=["total_counts" ], wspace=0.4)
sc.pl.umap(adata, color=[ "n_genes_by_counts"], wspace=0.4)
sc.pl.umap(adata, color=[ "clusters"], wspace=0.4)


plt.rcParams["figure.figsize"] = (8, 8)
sc.pl.spatial(adata, img_key="hires", color=["total_counts", "n_genes_by_counts"])


sc.pl.spatial(adata, img_key="hires", color="clusters", size=1.5)

sc.pl.spatial(adata,img_key="hires", color="clusters",
  groups=["1", "3"],crop_coord=[7000, 10000, 0, 6000],
  alpha=0.5,size=1.3,)


sc.tl.rank_genes_groups(adata, "clusters", method="t-test")
sc.pl.rank_genes_groups_heatmap(adata, groups="2", n_genes=10, groupby="clusters")

sc.pl.spatial(adata, img_key="hires", color=["clusters" ])

sc.pl.spatial(adata, img_key="hires", color=["COL1A2", "SYPL1"], alpha=0.7)

os.chdir('F:/spe.lesson\\s2\\')
adata.write("GSE194329.data.anndata.h5ad")
```





```bash
1. 首先打开一个新的Rstudio 指定miniconda的python环境

library(reticulate)
use_python("D:/miniconda/envs/spyder/") 
py_config()  
py_available()
repl_python()


5.空转数据运行
5.1 加载包
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

sc.logging.print_versions()
sc.set_figure_params(facecolor="white", figsize=(8, 8))
sc.settings.verbosity = 3
5.2  设置路径、读取数据

print(os.getcwd())
print(os.path )
adata = sc.read_visium('GBM1_spaceranger_out\\')
5.3 计算线粒体百分比

adata.var_names_make_unique()
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)


adata
5.4 质控展示

fig, axs = plt.subplots(1, 4, figsize=(15, 4))

fig, axs = plt.subplots(1, 4, figsize=(15, 4))
sns.histplot(adata.obs["total_counts"], kde=False, ax=axs[0])
sns.histplot(adata.obs["total_counts"][adata.obs["total_counts"] < 10000],
  kde=False, bins=40,ax=axs[1],)

sns.histplot(adata.obs["n_genes_by_counts"], kde=False, bins=60, ax=axs[2])
sns.histplot( adata.obs["n_genes_by_counts"][adata.obs["n_genes_by_counts"] < 4000],
  kde=False,bins=60,ax=axs[3],)
5.5 过滤细胞

sc.pp.filter_cells(adata, min_counts=5000)
sc.pp.filter_cells(adata, max_counts=35000)
adata = adata[adata.obs["pct_counts_mt"] < 20].copy()
print(f"#cells after MT filter: {adata.n_obs}")
sc.pp.filter_genes(adata, min_cells=10)



sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)
5.6 降维聚类

sc.pp.pca(adata)sc.pp.neighbors(adata)sc.tl.umap(adata)sc.tl.leiden(adata, key_added="clusters" )plt.rcParams["figure.figsize"] = (4, 4)sc.pl.umap(adata, color=["total_counts", "n_genes_by_counts", "clusters"], wspace=0.4)
5.7 画图展示

sc.pl.umap(adata, color=["total_counts" ], wspace=0.4)sc.pl.umap(adata, color=[ "n_genes_by_counts"], wspace=0.4)sc.pl.umap(adata, color=[ "clusters"], wspace=0.4)plt.rcParams["figure.figsize"] = (8, 8)sc.pl.spatial(adata, img_key="hires", color=["total_counts", "n_genes_by_counts"])
5.8 个性化展示

plt.rcParams["figure.figsize"] = (8, 8)
sc.pl.spatial(adata, img_key="hires", color=["total_counts", "n_genes_by_counts"])

sc.pl.spatial(adata, img_key="hires", color="clusters", size=1.5)

sc.pl.spatial(adata,img_key="hires", color="clusters",
  groups=["1", "3"],crop_coord=[7000, 10000, 0, 6000],
  alpha=0.5,size=1.3,)

sc.tl.rank_genes_groups(adata, "clusters", method="t-test")
sc.pl.rank_genes_groups_heatmap(adata, groups="2", n_genes=10, groupby="clusters")

sc.pl.spatial(adata, img_key="hires", color=["clusters" ])

sc.pl.spatial(adata, img_key="hires", color=["COL1A2", "SYPL1"], alpha=0.7)

adata.write("GSE194329.data.anndata.h5ad")


```


