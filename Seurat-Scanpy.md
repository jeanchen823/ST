# Seurat 和 scanpy 相互转换

```bash
1.加载R包
##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())
library(DropletUtils)
library(tibble)
library(Seurat)
library(jsonlite)
library(stringr)


2.设置路径读取数据

load("seurat.RData")

3.提取空间图像的缩放因子
Visium 数据通常包括高分辨率（hires）和低分辨率（lowres）的图像数据。这里，low_factor 和 hi_factor 分别存储了低分辨率和高分辨率图像的缩放因子。这些因子在分析空间转录组数据时非常重要，因为它们帮助调整图像大小以匹配转录组数据。
low_factor <- data@images[["anterior1"]]@scale.factors[["lowres"]]
hi_factor <- data@images[["anterior1"]]@scale.factors[["hires"]]


4.获取图像
## install.packages("magick")
images <- GetImage(data, mode = "raw")
images <- magick::image_read(images) 
image_height <- magick::image_info(images)$height


5.保存图像
#options(BioC_mirror="https://mirrors.westlake.edu.cn/bioconductor")
#BiocManager::install("EBImage")
EBImage::writeImage(magick::as_EBImage(images), file = './spatial/tissue_lowres_image.png')
images <- images %>% magick::image_resize(paste0('x',round(image_height/low_factor*hi_factor)))
EBImage::writeImage(magick::as_EBImage(images), file = './spatial/tissue_hires_image.png')


6.获取坐标信息
coordinates <- data@images$anterior1@coordinates
coordinates <- rownames_to_column(coordinates)
# coordinates[,1] <- str_split(coordinates[,1],'_',simplify = TRUE)[,2] #可以用这行修改barcodes前缀


7.保存坐标信息
write.table(coordinates,file = './spatial/tissue_positions_list.csv',quote = F,row.names = F,col.names = F,sep=',')


8.整理缩放因子
scale_factors <- data.frame(spot_diameter_fullres=data@images$anterior1@scale.factors$spot,
                            tissue_hires_scalef=data@images$anterior1@scale.factors$hires,
                            fiducial_diameter_fullres=data@images$anterior1@scale.factors$fiducial,
                            tissue_lowres_scalef=data@images$anterior1@scale.factors$lowres)

9.保存缩放因子
jsondata <- toJSON(scale_factors)
cat(str_remove_all(jsondata, '[\\[\\]]'), file = './spatial/scalefactors_json.json', fill = FALSE, labels = NULL, append = FALSE)


10.保存表达矩阵
write10xCounts(x = data@assays$Spatial$counts,
               #barcodes = str_split(colnames(P2_noninf@assays$Spatial$counts),'_',simplify = TRUE)[,2] , 
               # 如果barcodes已经被修改过带有前缀了，如果不想要前缀可以用上一行去掉，不过记得改coordinates里面的barcodes
               gene.id = rownames(data@assays$Spatial$counts),
               gene.symbol = rownames(data@assays$Spatial$counts),
               gene.type = "Gene Expression",
               overwrite = FALSE,
               type = 'HDF5',
               genome = "CRCh38",
               version = "3",
               chemistry = "Single Cell 3' v1",
               library.ids = "slice1",
               path = './filtered_feature_bc_matrix.h5')
               
保存为RDS
# 保存为RDS
saveRDS(data,file = './data.rds')


12.在spyder中读取数据
conda activate squidpy

打开spyder
spyder

读取数据
import matplotlib.pyplot as plt
import scanpy as sc
import pandas as pd
import seaborn as sns
import os
sc.settings.verbosity = 1             # verbosity errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
sc.settings.set_figure_params(dpi=80, frameon=False, figsize=(3, 3), facecolor='white')

print(os.getcwd())
os.chdir('seurat_scanpy\\')
## 更改完再检查一下有没有更改成功
print(os.getcwd())

 # PBMCs 3k(已经处理)
adata = sc.read_visium('seurat_scanpy\\',library_id="hauge")
adata.var_names_make_unique()

adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

adata

fig, axs = plt.subplots(1, 4, figsize=(15, 4))
sns.histplot(adata.obs["total_counts"], kde=False, ax=axs[0])
sns.histplot(
    adata.obs["total_counts"][adata.obs["total_counts"] < 10000],
    kde=False,
    bins=40,
    ax=axs[1],)
sns.histplot(adata.obs["n_genes_by_counts"], kde=False, bins=60, ax=axs[2])
sns.histplot(
    adata.obs["n_genes_by_counts"][adata.obs["n_genes_by_counts"] < 4000],
    kde=False,
    bins=60,
    ax=axs[3],)

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

```

```bash

conda activate squidpy
pip install pytables
conda install anaconda::pytables

1.读取数据
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

sc.logging.print_versions()
sc.set_figure_params(facecolor="white", figsize=(8, 8))
sc.settings.verbosity = 3


print(os.getcwd())

os.chdir('GBM1_spaceranger_out\\')

adata = sc.read_visium('GBM1_spaceranger_out\\')
adata.var_names_make_unique()
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

adata

fig, axs = plt.subplots(1, 4, figsize=(15, 4))
sns.histplot(adata.obs["total_counts"], kde=False, ax=axs[0])
sns.histplot(
    adata.obs["total_counts"][adata.obs["total_counts"] < 10000],
    kde=False,
    bins=40,
    ax=axs[1],)
sns.histplot(adata.obs["n_genes_by_counts"], kde=False, bins=60, ax=axs[2])
sns.histplot(
    adata.obs["n_genes_by_counts"][adata.obs["n_genes_by_counts"] < 4000],
    kde=False,
    bins=60,
    ax=axs[3],)

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
sc.tl.leiden(
    adata, key_added="clusters", directed=False, n_iterations=2
)

plt.rcParams["figure.figsize"] = (4, 4)
sc.pl.umap(adata, color=["total_counts", "n_genes_by_counts", "clusters"], wspace=0.4)

plt.rcParams["figure.figsize"] = (8, 8)
sc.pl.spatial(adata, img_key="hires", color=["total_counts", "n_genes_by_counts"])

sc.pl.spatial(adata, img_key="hires", color="clusters", size=1.5)


2.自定义函数
def huage_adata_info(adata):
    mat=pd.DataFrame(data=adata.X.todense(),index=adata.obs_names,columns=adata.var_names)
    mat.to_csv("mat.csv")
    meta=pd.DataFrame(data=adata.obs)
    meta.to_csv('metadata.tsv',sep="\t")
    cord=pd.DataFrame(data=adata.obsm['spatial'],index=adata.obs_names,columns=['x','y'])
    cord.to_csv('position_'+'spatial'+'.tsv',sep="\t")
    umap=pd.DataFrame(data=adata.obsm["X_umap"],index=adata.obs_names,columns=['x','y'])
    umap.to_csv('position_'+"X_umap"+'.tsv',sep="\t")

3.保存数据
os.chdir('scanpy_seurat\\')
huage_adata_info(adata)

4.查看结果

5.在R中读取数据
##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())
library(rhdf5)
library(dplyr)
library(data.table)
library(Matrix)
library(rjson)
library(Seurat)

##mydata =read.csv("exp1.csv")

mydata =read.csv("mat.csv")


6.整理表达矩阵
rownames(mydata)=mydata[,1]
mydata=mydata[,-1]
mat <- Matrix(t(mydata), sparse = TRUE)
exp=t(mydata)

7.创建seurat对象
obj <- CreateSeuratObject(exp, project ='Spatial', assay = 'Spatial',  meta.data=meta)


8.添加空间信息
tissue_lowres_image <- matrix(1, max(pos$y), max(pos$x))
tissue_positions_list <- data.frame(row.names = colnames(obj),
                                    tissue = 1,
                                    row = pos$y, col = pos$x,
                                    imagerow = pos$y, imagecol = pos$x)
scalefactors_json <- toJSON(list(fiducial_diameter_fullres = 1,
                                 tissue_hires_scalef = 1,
                                 tissue_lowres_scalef = 1))
seurat_spatialObj <- obj
generate_spatialObj <- function(image, scale.factors, tissue.positions, filter.matrix = TRUE)
{
  if (filter.matrix) {
    tissue.positions <- tissue.positions[which(tissue.positions$tissue == 1), , drop = FALSE]
  }
  
  unnormalized.radius <- scale.factors$fiducial_diameter_fullres * scale.factors$tissue_lowres_scalef
  spot.radius <- unnormalized.radius / max(dim(image))
  return(new(Class = 'VisiumV1',
             image = image,
             scale.factors = scalefactors(spot = scale.factors$tissue_hires_scalef,
                                          fiducial = scale.factors$fiducial_diameter_fullres,
                                          hires = scale.factors$tissue_hires_scalef,
                                          lowres = scale.factors$tissue_lowres_scalef),
             coordinates = tissue.positions,
             spot.radius = spot.radius))
}

spatialObj <- generate_spatialObj(image = tissue_lowres_image,
                                  scale.factors = fromJSON(scalefactors_json),
                                  tissue.positions = tissue_positions_list)

spatialObj <- spatialObj[Cells(seurat_spatialObj)]
DefaultAssay(spatialObj) <- 'Spatial'
seurat_spatialObj[['slice1']] <- spatialObj

9.结果展示
seurat_spatialObj=NormalizeData(seurat_spatialObj)

SpatialFeaturePlot(seurat_spatialObj, features = c( "MS4A1"))
SpatialPlot(seurat_spatialObj, features = "MS4A1")
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


print(os.getcwd())


adata = sc.read_visium('GBM1_spaceranger_out\\')
adata.var_names_make_unique()
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

adata

fig, axs = plt.subplots(1, 4, figsize=(15, 4))
sns.histplot(adata.obs["total_counts"], kde=False, ax=axs[0])
sns.histplot(
    adata.obs["total_counts"][adata.obs["total_counts"] < 10000],
    kde=False,
    bins=40,
    ax=axs[1],)
sns.histplot(adata.obs["n_genes_by_counts"], kde=False, bins=60, ax=axs[2])
sns.histplot(
    adata.obs["n_genes_by_counts"][adata.obs["n_genes_by_counts"] < 4000],
    kde=False,
    bins=60,
    ax=axs[3],)

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
sc.tl.leiden(
    adata, key_added="clusters", directed=False, n_iterations=2
)

plt.rcParams["figure.figsize"] = (4, 4)
sc.pl.umap(adata, color=["total_counts", "n_genes_by_counts", "clusters"], wspace=0.4)

plt.rcParams["figure.figsize"] = (8, 8)
sc.pl.spatial(adata, img_key="hires", color=["total_counts", "n_genes_by_counts"])

sc.pl.spatial(adata, img_key="hires", color="clusters", size=1.5)



def huage_adata_info(adata):
    mat=pd.DataFrame(data=adata.X.todense(),index=adata.obs_names,columns=adata.var_names)
    mat.to_csv("mat.csv")
    meta=pd.DataFrame(data=adata.obs)
    meta.to_csv('metadata.tsv',sep="\t")
    cord=pd.DataFrame(data=adata.obsm['spatial'],index=adata.obs_names,columns=['x','y'])
    cord.to_csv('position_'+'spatial'+'.tsv',sep="\t")
    umap=pd.DataFrame(data=adata.obsm["X_umap"],index=adata.obs_names,columns=['x','y'])
    umap.to_csv('position_'+"X_umap"+'.tsv',sep="\t")



huage_adata_info(adata)
```

```bash

##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())
library(rhdf5)
library(dplyr)
library(data.table)
library(Matrix)
library(rjson)
library(Seurat)


##mydata =read.csv("exp1.csv")

mydata =read.csv("mat.csv")

meta <- read.table("metadata.tsv",sep="\t",header=T,row.names=1)
pos <- read.table("position_spatial.tsv",sep="\t",header=T,row.names=1)
umap<- read.table("position_X_umap.tsv",sep="\t",header=T,row.names=1)
rownames(mydata)=mydata[,1]
mydata=mydata[,-1]
mat <- Matrix(t(mydata), sparse = TRUE)
exp=t(mydata)

obj <- CreateSeuratObject(exp, project ='Spatial', assay = 'Spatial',  meta.data=meta)


tissue_lowres_image <- matrix(1, max(pos$y), max(pos$x))
tissue_positions_list <- data.frame(row.names = colnames(obj),
                                    tissue = 1,
                                    row = pos$y, col = pos$x,
                                    imagerow = pos$y, imagecol = pos$x)
scalefactors_json <- toJSON(list(fiducial_diameter_fullres = 1,
                                 tissue_hires_scalef = 1,
                                 tissue_lowres_scalef = 1))
seurat_spatialObj <- obj
generate_spatialObj <- function(image, scale.factors, tissue.positions, filter.matrix = TRUE)
{
  if (filter.matrix) {
    tissue.positions <- tissue.positions[which(tissue.positions$tissue == 1), , drop = FALSE]
  }
  
  unnormalized.radius <- scale.factors$fiducial_diameter_fullres * scale.factors$tissue_lowres_scalef
  spot.radius <- unnormalized.radius / max(dim(image))
  return(new(Class = 'VisiumV1',
             image = image,
             scale.factors = scalefactors(spot = scale.factors$tissue_hires_scalef,
                                          fiducial = scale.factors$fiducial_diameter_fullres,
                                          hires = scale.factors$tissue_hires_scalef,
                                          lowres = scale.factors$tissue_lowres_scalef),
             coordinates = tissue.positions,
             spot.radius = spot.radius))
}

spatialObj <- generate_spatialObj(image = tissue_lowres_image,
                                  scale.factors = fromJSON(scalefactors_json),
                                  tissue.positions = tissue_positions_list)

spatialObj <- spatialObj[Cells(seurat_spatialObj)]
DefaultAssay(spatialObj) <- 'Spatial'
seurat_spatialObj[['slice1']] <- spatialObj





seurat_spatialObj=NormalizeData(seurat_spatialObj)

SpatialFeaturePlot(seurat_spatialObj, features = c( "MS4A1"))
SpatialPlot(seurat_spatialObj, features = "MS4A1")

```

