# Seurat 

```bash
https://github.com/KPMP/Cell-State-Atlas-2022/blob/develop/SourceByTechnology/Visium/start_spatial_objects.R

https://satijalab.org/seurat/articles/spatial_vignette



# 加载R包
#options(BioC_mirror="https://mirrors.westlake.edu.cn/bioconductor")
#BiocManager::install("DropletUtils")
##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(DropletUtils)
library(tibble)
library(Seurat)
library(jsonlite)
library(stringr)


# 读入10X Visium数据
data <- Load10X_Spatial('GBM1_spaceranger_out/')




2.数据处理

We note that the variance in molecular counts / spot can be substantial for spatial datasets, particularly if there are differences in cell density across the tissue. We see substantial heterogeneity here, which requires effective normalization.
plot1 <- VlnPlot(brain, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(brain, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)




3.基因可视化

我们拥有探索和交互空间数据固有的视觉特性的功能，Seurat中的SpatialFeaturePlot()函数扩展了FeaturePlot()，并且可以在组织组织学上覆盖分子数据。
SpatialFeaturePlot(brain, features = c("CD8A", "CD4"))

p1 <- SpatialFeaturePlot(brain, features = "GAPDH", pt.size.factor = 1)
p2 <- SpatialFeaturePlot(brain, features ="GAPDH", alpha = c(0.1, 1))
p1 + p2





4.降维、聚类、可视化
brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)
brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
brain <- FindClusters(brain, verbose = FALSE)
brain <- RunUMAP(brain, reduction = "pca", dims = 1:30)

然后，我们可以在UMAP空间中(使用DimPlot())或使用SpatialDimPlot()将聚类结果叠加在图像上，从而将聚类结果可视化。
p1 <- DimPlot(brain, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(brain, label = TRUE, label.size = 3)
p1 + p2



SpatialDimPlot(brain, cells.highlight = CellsByIdentities(object = brain, 
    idents = c(2, 1, 4, 3,5, 8)), facet.highlight = TRUE, ncol = 3)



5.交互式画图

SpatialDimPlot(brain, interactive = TRUE)


SpatialFeaturePlot(brain, features = "GAPDH", interactive = TRUE)



LinkedDimPlot(brain)






6.识别空间变异基因
Seurat提供了两种思路的空间变异基因识别方法，一种是用注释好的组织区域（手工注释或聚类）做差异分析；另一种方法是基于基因表达的空间自相关性分析，这种方法不需要提供组织注释或聚类信息。

de_markers <- FindMarkers(brain, ident.1 = 5, ident.2 = 6)
SpatialFeaturePlot(object = brain, features = rownames(de_markers)[1:3], alpha = c(0.1, 1), ncol = 3)


brain <- FindSpatiallyVariableFeatures(brain, assay = "SCT", features = VariableFeatures(brain)[1:1000],
                                       selection.method = "moransi")


top.features <- head(SpatiallyVariableFeatures(brain, selection.method = "moransi"), 6)
SpatialFeaturePlot(brain, features = top.features, ncol = 3, alpha = c(0.1, 1))



7.子集化切片区域

与单单元格对象一样，您可以对对象进行子集化，以便将重点放在数据的子集上。这里，我们大致对额叶皮层进行了子集分析。这一过程也有助于将这些数据与下一节的皮质scRNA-seq数据集整合。首先，我们取集群的一个子集，然后根据确切的位置进一步分割。在亚子集之后，我们可以在完整的图像上或裁剪后的图像上可视化皮质细胞。
##Subset out anatomical regions
cortex <- subset(brain, idents = c(1, 2, 3, 4, 6, 7))
# now remove additional cells, use SpatialDimPlots to visualize what to remove
# SpatialDimPlot(cortex,cells.highlight = WhichCells(cortex, expression = image_imagerow > 400
# | image_imagecol < 150))
#cortex <- subset(cortex, anterior1_imagerow > 400 | anterior1_imagecol < 150, invert = TRUE)
#cortex <- subset(cortex, anterior1_imagerow > 275 & anterior1_imagecol > 370, invert = TRUE)
#cortex <- subset(cortex, anterior1_imagerow > 250 & anterior1_imagecol > 440, invert = TRUE)

p1 <- SpatialDimPlot(cortex, crop = TRUE, label = TRUE)
p2 <- SpatialDimPlot(cortex, crop = FALSE, label = TRUE, pt.size.factor = 1, label.size = 3)
p1 + p2
```



```bash
1.处理空转数据
##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())
#InstallData("stxBrain")
#brain <- LoadData("stxBrain", type = "anterior1")
#save(brain,file = "brain.rdata")
load("brain.rdata")
plot1 <- VlnPlot(brain, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(brain, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)
brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)
brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
brain <- FindClusters(brain, verbose = FALSE)
brain <- RunUMAP(brain, reduction = "pca", dims = 1:30)
p1 <- DimPlot(brain, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(brain, label = TRUE, label.size = 3)
p1 + p2


allen_reference <- readRDS("allen_cortex.rds")


library(dplyr)
allen_reference <- SCTransform(allen_reference, ncells = 3000, verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30)
  
  

2.提取皮质区域子集

cortex <- subset(brain, idents = c(1, 2, 3, 4, 6, 7))
# now remove additional cells, use SpatialDimPlots to visualize what to remove
# SpatialDimPlot(cortex,cells.highlight = WhichCells(cortex, expression = image_imagerow > 400
# | image_imagecol < 150))
cortex <- subset(cortex, anterior1_imagerow > 400 | anterior1_imagecol < 150, invert = TRUE)
cortex <- subset(cortex, anterior1_imagerow > 275 & anterior1_imagecol > 370, invert = TRUE)
cortex <- subset(cortex, anterior1_imagerow > 250 & anterior1_imagecol > 440, invert = TRUE)

p1 <- SpatialDimPlot(cortex, crop = TRUE, label = TRUE)
p2 <- SpatialDimPlot(cortex, crop = FALSE, label = TRUE, pt.size.factor = 1, label.size = 3)
p1 + p2



3.在提取子集之后，我们重对皮质进行归一化
# After subsetting, we renormalize cortex
cortex <- SCTransform(cortex, assay = "Spatial", verbose = FALSE) %>%
    RunPCA(verbose = FALSE)
# 注释存储在metadata数据的  'subclass' 列中
DimPlot(allen_reference, group.by = "subclass", label = TRUE)


4.用单细胞来注释空转
anchors <- FindTransferAnchors(reference = allen_reference, query = cortex, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = allen_reference$subclass, prediction.assay = TRUE,
    weight.reduction = cortex[["pca"]], dims = 1:30)
cortex[["predictions"]] <- predictions.assay

DefaultAssay(cortex) <- "predictions"
SpatialFeaturePlot(cortex, features = c("L2/3 IT", "L4"), pt.size.factor = 1.6, ncol = 2, crop = TRUE)


cortex <- FindSpatiallyVariableFeatures(cortex, assay = "predictions", selection.method = "moransi",
                                        features = rownames(cortex), r.metric = 5, slot = "data")
top.clusters <- head(SpatiallyVariableFeatures(cortex, selection.method = "moransi"), 4)
SpatialPlot(object = cortex, features = top.clusters, ncol = 2)



SpatialFeaturePlot(cortex, features = c("Astro",   "Oligo"), 
      pt.size.factor = 1,ncol = 2, crop = FALSE, alpha = c(0.1, 1))
      

 
多个空转样本的整合

这个老鼠大脑的数据集包含了另一个与大脑另一半相对应的切片。这里我们读入它并执行相同的初始归一化。

#brain2 <- LoadData("stxBrain", type = "posterior1")
#save(brain2,file = "brain2.rdata")
load("brain2.rdata")
brain2 <- SCTransform(brain2, assay = "Spatial", verbose = FALSE)

In order to work with multiple slices in the same Seurat object, we provide the merge function.

brain.merge <- merge(brain, brain2)
这样就可以在RNA表达数据上进行联合降维和聚类。
DefaultAssay(brain.merge) <- "SCT"
VariableFeatures(brain.merge) <- c(VariableFeatures(brain), VariableFeatures(brain2))
brain.merge <- RunPCA(brain.merge, verbose = FALSE)
brain.merge <- FindNeighbors(brain.merge, dims = 1:30)
brain.merge <- FindClusters(brain.merge, verbose = FALSE)
brain.merge <- RunUMAP(brain.merge, dims = 1:30)
最后，数据可以在单个UMAP图中共同可视化。
DimPlot(brain.merge, reduction = "umap", group.by = c("ident", "orig.ident"))


SpatialDimPlot(brain.merge)



SpatialFeaturePlot(brain.merge, features = c("Hpca", "Plp1"))


 
画图小技巧

1.去除point黑色边圈

默认的可视化会把spot加上黑边，当spot点很密集的时候整个图片会显示的比较暗，可以设置stroke=NA来使图片变得明亮。
p1 <- SpatialDimPlot(brain)
p2 <- SpatialDimPlot(brain, stroke=NA) 
p2 <- p2 + guides(color=guide_legend(override.aes = list(size=8), ncol=2))
p2 <- p2 + theme_gray()+xlab("")+ylab("")
p2 <- p2 + theme(axis.text = element_blank(),axis.ticks=element_blank())
p1+p2

2.用seurat单细胞的函数实现空转spot分布
FeaturePlot, DimPlot是单细胞数据可视化的函数
SpatialDimPlot, SpatialFeaturePlot是空转数据可视化的函数
为了更加个性化的展示我们的空转数据，我们使用单细胞的函数(FeaturePlot, DimPlot)对空转的数据进行可视化。
spatial_corr <- brain@images$anterior1@coordinates[,c('col', 'row')]
colnames(spatial_corr) <- c('s_1', 's_2')
spatial_corr <- as.matrix(spatial_corr)
brain[["spatial"]] <- CreateDimReducObject(embeddings=spatial_corr, key = "s_", assay = "Spatial")

p1 <- FeaturePlot(rds, features = "nCount_Spatial", reduction='spatial') + 
  scale_y_reverse() + 
  scale_colour_viridis(option="inferno") + theme_void()
p2 <- FeaturePlot(rds, features = "Hpca", reduction='spatial') +
  scale_y_reverse() + 
  scale_colour_viridis(option="D") + theme_void()

mycolor <- c("#927A66","#DBAEA4","#97A4AB","#21A69A","#3A6688","#C98882","#A593A7","#CE8662","#B04929","#A59487","#747A87","#F2DCD5","#A7ABB3","#C1D1CD","#DCA5A5","#BC8B83")

p3 <- DimPlot(rds, reduction='spatial', cols=mycolor) + 
  guides(color=guide_legend(override.aes = list(size=6), ncol=2)) +
  theme_void() + scale_y_reverse()

p1+p2+p3+plot_layout(ncol=3, nrow=1)



https://www.nature.com/articles/s41588-023-01570-0
首先自定义画图函数
library(Seurat)
library(ggplot2)
library(viridis)
library(patchwork)
library(RColorBrewer)
library(paletteer)
ng.plot <- function(RDS, feature, ad=1, limits=1){
  viridis_plasma_light_high <- as.vector(x = paletteer_c(palette = "viridis::inferno", n = 250, direction = 1))
  viridis_plasma_light_high <- c( rep("black", ad), viridis_plasma_light_high)

  p <- FeaturePlot(RDS, features = feature, reduction='spatial')
  p <- p + theme_void()+ theme(
    axis.ticks=element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.line = element_blank(),
    panel.border = element_rect(color = "white", fill = NA, size =2),
  ) + 
    DarkTheme() +
    xlab(NULL) + 
    ylab(NULL)
  if (length(limits)==1){
    p <- p + scale_colour_gradientn(colours=viridis_plasma_light_high, na.value = "black") 
  }else{
    p <- p + scale_colour_gradientn(colours=viridis_plasma_light_high, na.value = "black",limits=limits)
  }

  return (p)
}

table(brain$region)
p1 <- ng.plot(brain,  c("Hpca"))
p2 <- ng.plot(brain,  c( "Plp1"))
p1+p2+ plot_layout(ncol=2, nrow=1)


4.SpatialDimPlot 自定义颜色
Idents(brain) = brain$seurat_clusters 

colors <- c("#927A66","#DBAEA4","#97A4AB","#21A69A","#3A6688","#C98882","#A593A7","#CE8662","#B04929","#A59487","#747A87","#F2DCD5","#A7ABB3","#C1D1CD","#DCA5A5","#BC8B83")


table(brain$seurat_clusters)
colors=colors[1:15]
names(colors) <- Idents(brain) %>% levels() ## 命名非常重要
colnames(brain@meta.data)
SpatialDimPlot(brain,
  images = "anterior1",
  group.by = "seurat_clusters",
  pt.size.factor = 1.15,
  label = TRUE, label.size = 6,
  repel = TRUE,combine = FALSE,
  cols = colors)
  



```



```bash

### 读取数据 ###
# 加载R包
#options(BioC_mirror="https://mirrors.westlake.edu.cn/bioconductor")
#BiocManager::install("DropletUtils")
##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(DropletUtils)
library(tibble)
library(Seurat)
library(jsonlite)
library(stringr)


# 读入10X Visium数据
data <- Load10X_Spatial('GBM1_spaceranger_out/')

brain=data
##Data preprocessing
plot1 <- VlnPlot(brain, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(brain, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)
#Gene expression visualization
SpatialFeaturePlot(brain, features = c("CD8A", "CD4"))


library(ggplot2)
plot <- SpatialFeaturePlot(brain, features = c("CD8A")) + theme(legend.text = element_text(size = 0),
                                                                 legend.title = element_text(size = 20), legend.key.size = unit(1, "cm"))
jpeg(filename = "./spatial_vignette_ttr.jpg", height = 700, width = 1200, quality = 50)
print(plot)
dev.off()


p1 <- SpatialFeaturePlot(brain, features = "GAPDH", pt.size.factor = 1)
p2 <- SpatialFeaturePlot(brain, features ="GAPDH", alpha = c(0.1, 1))
p1 + p2

##Dimensionality reduction, clustering, and visualization

brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)
brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
brain <- FindClusters(brain, verbose = FALSE)
brain <- RunUMAP(brain, reduction = "pca", dims = 1:30)


p1 <- DimPlot(brain, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(brain, label = TRUE, label.size = 3)
p1 + p2



SpatialDimPlot(brain, cells.highlight = CellsByIdentities(object = brain, 
                                                          idents = c(2, 1, 4, 3,5, 8)), facet.highlight = TRUE, ncol = 3)



#Interactive plotting
SpatialDimPlot(brain, interactive = TRUE)

SpatialFeaturePlot(brain, features = "GAPDH", interactive = TRUE)

LinkedDimPlot(brain)

##Identification of Spatially Variable Features
de_markers <- FindMarkers(brain, ident.1 = 5, ident.2 = 6)
SpatialFeaturePlot(object = brain, features = rownames(de_markers)[1:3], alpha = c(0.1, 1), ncol = 3)

brain <- FindSpatiallyVariableFeatures(brain, assay = "SCT", features = VariableFeatures(brain)[1:1000],
                                       selection.method = "moransi")

##Now we visualize the expression of the top 6 features identified by this measure.
SpatiallyVariableFeatures(brain )
top.features <- head(SpatiallyVariableFeatures(brain, selection.method = "moransi"), 6)
SpatialFeaturePlot(brain, features = top.features, ncol = 3, alpha = c(0.1, 1))


##Subset out anatomical regions
cortex <- subset(brain, idents = c(1, 2, 3, 4, 6, 7))
# now remove additional cells, use SpatialDimPlots to visualize what to remove
# SpatialDimPlot(cortex,cells.highlight = WhichCells(cortex, expression = image_imagerow > 400
# | image_imagecol < 150))
cortex <- subset(cortex, anterior1_imagerow > 400 | anterior1_imagecol < 150, invert = TRUE)
cortex <- subset(cortex, anterior1_imagerow > 275 & anterior1_imagecol > 370, invert = TRUE)
cortex <- subset(cortex, anterior1_imagerow > 250 & anterior1_imagecol > 440, invert = TRUE)

p1 <- SpatialDimPlot(cortex, crop = TRUE, label = TRUE)
p2 <- SpatialDimPlot(cortex, crop = FALSE, label = TRUE, pt.size.factor = 1, label.size = 3)
p1 + p2
```



```bash

##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())
#InstallData("stxBrain")
#brain <- LoadData("stxBrain", type = "anterior1")
#save(brain,file = "brain.rdata")
load("brain.rdata")
plot1 <- VlnPlot(brain, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(brain, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)
brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)
brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
brain <- FindClusters(brain, verbose = FALSE)
brain <- RunUMAP(brain, reduction = "pca", dims = 1:30)
p1 <- DimPlot(brain, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(brain, label = TRUE, label.size = 3)
p1 + p2


allen_reference <- readRDS("allen_cortex.rds")

library(dplyr)
allen_reference <- SCTransform(allen_reference, ncells = 3000, verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30)

cortex <- subset(brain, idents = c(1, 2, 3, 4, 6, 7))
# now remove additional cells, use SpatialDimPlots to visualize what to remove
# SpatialDimPlot(cortex,cells.highlight = WhichCells(cortex, expression = image_imagerow > 400
# | image_imagecol < 150))
cortex <- subset(cortex, anterior1_imagerow > 400 | anterior1_imagecol < 150, invert = TRUE)
cortex <- subset(cortex, anterior1_imagerow > 275 & anterior1_imagecol > 370, invert = TRUE)
cortex <- subset(cortex, anterior1_imagerow > 250 & anterior1_imagecol > 440, invert = TRUE)

p1 <- SpatialDimPlot(cortex, crop = TRUE, label = TRUE)
p2 <- SpatialDimPlot(cortex, crop = FALSE, label = TRUE, pt.size.factor = 1, label.size = 3)
p1 + p2


# After subsetting, we renormalize cortex
cortex <- SCTransform(cortex, assay = "Spatial", verbose = FALSE) %>%
  RunPCA(verbose = FALSE)
# the annotation is stored in the 'subclass' column of object metadata
DimPlot(allen_reference, group.by = "subclass", label = TRUE)

anchors <- FindTransferAnchors(reference = allen_reference, query = cortex, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = allen_reference$subclass, prediction.assay = TRUE,
                                  weight.reduction = cortex[["pca"]], dims = 1:30)
cortex[["predictions"]] <- predictions.assay

DefaultAssay(cortex) <- "predictions"
SpatialFeaturePlot(cortex, features = c("L2/3 IT", "L4"), pt.size.factor = 1.6, ncol = 2, crop = TRUE)


cortex <- FindSpatiallyVariableFeatures(cortex, assay = "predictions", selection.method = "moransi",
                                        features = rownames(cortex), r.metric = 5, slot = "data")
top.clusters <- head(SpatiallyVariableFeatures(cortex, selection.method = "moransi"), 4)
SpatialPlot(object = cortex, features = top.clusters, ncol = 2)


SpatialFeaturePlot(cortex, features = c("Astro",   "Oligo"), 
      pt.size.factor = 1,ncol = 2, crop = FALSE, alpha = c(0.1, 1))


#brain2 <- LoadData("stxBrain", type = "posterior1")
#save(brain2,file = "brain2.rdata")
load("brain2.rdata")
brain2 <- SCTransform(brain2, assay = "Spatial", verbose = FALSE)

brain.merge <- merge(brain, brain2)

DefaultAssay(brain.merge) <- "SCT"
VariableFeatures(brain.merge) <- c(VariableFeatures(brain), VariableFeatures(brain2))
brain.merge <- RunPCA(brain.merge, verbose = FALSE)
brain.merge <- FindNeighbors(brain.merge, dims = 1:30)
brain.merge <- FindClusters(brain.merge, verbose = FALSE)
brain.merge <- RunUMAP(brain.merge, dims = 1:30)

DimPlot(brain.merge, reduction = "umap", group.by = c("ident", "orig.ident"))

SpatialDimPlot(brain.merge)

SpatialFeaturePlot(brain.merge, features = c("Hpca", "Plp1"))




p1 <- SpatialFeaturePlot(brain, features ="Hpca")
p2 <- p1 + theme_gray()+xlab("")+ylab("")+theme(axis.text = element_blank(),axis.ticks=element_blank())
p1+p2


p1 <- SpatialDimPlot(brain)
p2 <- SpatialDimPlot(brain, stroke=NA) 
p2 <- p2 + guides(color=guide_legend(override.aes = list(size=8), ncol=2))
p2 <- p2 + theme_gray()+xlab("")+ylab("")
p2 <- p2 + theme(axis.text = element_blank(),axis.ticks=element_blank())
p1+p2


spatial_corr <- brain@images$anterior1@coordinates[,c('col', 'row')]
colnames(spatial_corr) <- c('s_1', 's_2')
spatial_corr <- as.matrix(spatial_corr)
brain[["spatial"]] <- CreateDimReducObject(embeddings=spatial_corr, key = "s_", assay = "Spatial")



p1 <- FeaturePlot(rds, features = "nCount_Spatial", reduction='spatial') + 
  scale_y_reverse() + 
  scale_colour_viridis(option="inferno") + theme_void()
p2 <- FeaturePlot(rds, features = "Hpca", reduction='spatial') +
  scale_y_reverse() + 
  scale_colour_viridis(option="D") + theme_void()

mycolor <- c("#927A66","#DBAEA4","#97A4AB","#21A69A","#3A6688","#C98882","#A593A7","#CE8662","#B04929","#A59487","#747A87","#F2DCD5","#A7ABB3","#C1D1CD","#DCA5A5","#BC8B83")

p3 <- DimPlot(rds, reduction='spatial', cols=mycolor) + 
  guides(color=guide_legend(override.aes = list(size=6), ncol=2)) +
  theme_void() + scale_y_reverse()

p1+p2+p3+plot_layout(ncol=3, nrow=1)




library(Seurat)
library(ggplot2)
library(viridis)
library(patchwork)
library(RColorBrewer)
library(paletteer)
ng.plot <- function(RDS, feature, ad=1, limits=1){
  viridis_plasma_light_high <- as.vector(x = paletteer_c(palette = "viridis::inferno", n = 250, direction = 1))
  viridis_plasma_light_high <- c( rep("black", ad), viridis_plasma_light_high)

  p <- FeaturePlot(RDS, features = feature, reduction='spatial')
  p <- p + theme_void()+ theme(
    axis.ticks=element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.line = element_blank(),
    panel.border = element_rect(color = "white", fill = NA, size =2),
  ) + 
    DarkTheme() +
    xlab(NULL) + 
    ylab(NULL)
  if (length(limits)==1){
    p <- p + scale_colour_gradientn(colours=viridis_plasma_light_high, na.value = "black") 
  }else{
    p <- p + scale_colour_gradientn(colours=viridis_plasma_light_high, na.value = "black",limits=limits)
  }

  return (p)
}

table(brain$region)
p1 <- Newf(brain,  c("Hpca"))
p2 <- Newf(brain,  c( "Plp1"))
p1+p2+ plot_layout(ncol=2, nrow=1)



Idents(brain) = brain$seurat_clusters 

colors <- c("#927A66","#DBAEA4","#97A4AB","#21A69A","#3A6688","#C98882","#A593A7","#CE8662","#B04929","#A59487","#747A87","#F2DCD5","#A7ABB3","#C1D1CD","#DCA5A5","#BC8B83")


table(brain$seurat_clusters)
colors=colors[1:15]
names(colors) <- Idents(brain) %>% levels() ## 命名非常重要
colnames(brain@meta.data)
SpatialDimPlot(brain,
  images = "anterior1",
  group.by = "seurat_clusters",
  pt.size.factor = 1.15,
  label = TRUE, label.size = 6,
  repel = TRUE,combine = FALSE,
  cols = colors)
```

