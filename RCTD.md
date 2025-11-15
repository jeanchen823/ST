# RCTD

```bash
# install.packages("devtools")
devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
library(spacexr)
library(Matrix)
```

```bash
##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())

library(spacexr)
library(Matrix)
library(doParallel)

library(ggplot2)

setwd("Doublet_mode/")

### Load in/preprocess your data, this might vary based on your file type
# directory for the reference
refdir <- system.file("extdata",'Reference/Vignette',package = 'spacexr') 
# load in counts matrix
counts <- read.csv(file.path(refdir,"dge.csv")) 
# Move first column to rownames
rownames(counts) <- counts[,1]
counts[,1] <- NULL
counts[1:4,1:4]
```

```bash
# load in meta_data (barcodes, clusters, and nUMI)
meta_data <- read.csv(file.path(refdir,"meta_data.csv")) 
head(meta_data)
```

```bash

# create cell_types named list
cell_types <- meta_data$cluster
names(cell_types) <- meta_data$barcode 
# convert to factor data type
cell_types <- as.factor(cell_types) 
head(cell_types)
```

```bash

# create nUMI named list
nUMI <- meta_data$nUMI
names(nUMI) <- meta_data$barcode 
head(nUMI)
### Create the Reference object
reference <- Reference(counts, cell_types, nUMI)
str(reference)
## Examine reference object (optional)
#observe Digital Gene Expression matrix
print(dim(reference@counts)) 
#> [1] 384 475

#number of occurences for each cell type
table(reference@cell_types) 

## Save RDS object (optional)
saveRDS(reference, file.path(refdir,'SCRef.rds'))
```

```bash

##  读取空间数据
# directory for sample Slide-seq dataset
datadir <- system.file("extdata",'SpatialRNA/Vignette',package = 'spacexr') 
counts <- read.csv(file.path(datadir,"MappedDGEForR.csv"))
rownames(counts) <- counts[,1]; counts[,1] <- NULL
counts[1:4,1:4]

coords <- read.csv(file.path(datadir,"BeadLocationsForR.csv"))
rownames(coords) <- coords$barcodes
coords$barcodes <- NULL
head(coords)

# In this case, total counts per pixel is nUMI
nUMI <- colSums(counts) 
head(nUMI)

### Create SpatialRNA object
puck <- SpatialRNA(coords, counts, nUMI)
str(puck)


### Create SpatialRNA object
puck <- SpatialRNA(coords, counts, nUMI)
str(puck)

## Examine SpatialRNA object (optional)
print(dim(puck@counts))

# histogram of log_2 nUMI
hist(log(puck@nUMI,2))
```



```bash
print(head(puck@coords)) # start of coordinate data.frame
barcodes <- colnames(puck@counts) # pixels to be used (a list of barcode names). 

# This list can be restricted if you want to crop the puck e.g. 
# puck <- restrict_puck(puck, barcodes) provides a basic plot of the nUMI of each pixel
# on the plot:
#使用x，y坐标展示每个像素中的UMI分布
p <- plot_puck_continuous(puck, barcodes, puck@nUMI, ylimit = c(0,round(quantile(puck@nUMI,0.9))), title ='plot of nUMI') 
p
```



```bash

myRCTD <- create.RCTD(puck, reference, max_cores = 10)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet')
str(myRCTD)
```



```bash

# normalize the cell type proportions to sum to 1.
norm_weights = normalize_weights(results$weights) 

#list of cell type names
cell_type_names <- myRCTD@cell_type_info$info[[2]] 
spatialRNA <- myRCTD@spatialRNA

## you may change this to a more accessible directory on your computer.
resultsdir <- 'RCTD_Plots' 
dir.create(resultsdir)
```



```bash
# make the plots 
# 绘制full_mode模式下每种细胞类型的可信权重 (saved as 
# 'results/cell_type_weights_unthreshold.pdf')
plot_weights(cell_type_names, spatialRNA, resultsdir, norm_weights) 

# 绘制full_mode模式下每种细胞类型的权重
# 这里每种细胞类型一幅图，点表示空间上的一个像素或者spot，颜色为权重 (saved as
# 'results/cell_type_weights.pdf')
plot_weights_unthreshold(cell_type_names, spatialRNA, resultsdir, norm_weights) 

# 绘制full_mode模式下每种细胞类型预测到的spots数 (saved as 
# 'results/cell_type_occur.pdf')
plot_cond_occur(cell_type_names, resultsdir, norm_weights, spatialRNA)

# 绘制doublet_mode模式下每种细胞类型的权重 (saved as 
# 'results/cell_type_weights_doublets.pdf')
plot_weights_doublet(cell_type_names, spatialRNA, resultsdir, results$weights_doublet, results$results_df) 
```

```bash

library(Giotto)
library(Seurat)
library(tidyverse)
library(patchwork)
setwd("NC/")

data <- Load10X_Spatial("10x_brain/")

brain=data
##Data preprocessing
plot1 <- VlnPlot(brain, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(brain, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)
#Gene expression visualization

SpatialFeaturePlot(brain, features = c("Sox2", "Gapdh"))

library(ggplot2)
plot <- SpatialFeaturePlot(brain, features = c("Sox2")) + theme(legend.text = element_text(size = 0),
                                                                legend.title = element_text(size = 20), legend.key.size = unit(1, "cm"))
jpeg(filename = "spatial_vignette_ttr.jpg", height = 700, width = 1200, quality = 50)
print(plot)
dev.off()


p1 <- SpatialFeaturePlot(brain, features = "Sox2", pt.size.factor = 1)
p2 <- SpatialFeaturePlot(brain, features ="Sox2", alpha = c(0.1, 1))
p1 + p2

##Dimensionality reduction, clustering, and visualization

brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)
brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
brain <- FindClusters(brain, verbose = FALSE)
brain <- RunUMAP(brain, reduction = "pca", dims = 1:30)

p1 <- DimPlot(brain, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(brain, label = TRUE, label.size = 3)
p1 + p2
```



```bash

library(data.table)
exp=fread("brain_sc_expression_matrix.txt",  data.table = F)
rownames(exp)=exp[,1]
exp=exp[,-1]

meta=read.csv("brain_sc_metadata.csv")
meta=meta[,c(1,2,3,4,12)] 
rownames(meta)=meta$X

sc=CreateSeuratObject(exp,meta.data = meta)

table(meta$Class)
colnames(meta)
scRNA=sc
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^mt-")
# 计算细胞中核糖体基因比例
#scRNA[["percent.rb"]] <- PercentageFeatureSet(scRNA, pattern = "^RP[LS]")
# 质控
scRNA <- subset(scRNA, nCount_RNA>1000&nFeature_RNA>200&percent.mt<25)
# 标准化
scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)
# 鉴定高变基因
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 3000)
# 尺度变换
scRNA <- ScaleData(scRNA, features = rownames(scRNA))
# 降维
library(harmony)
scRNA <- RunPCA(scRNA, features = VariableFeatures(object = scRNA))
scRNA <- RunHarmony(scRNA, group.by.vars="orig.ident", max.iter.harmony = 20)
scRNA <- RunUMAP(scRNA,reduction = "harmony", dims = 1:30, verbose = FALSE)
## 确定分辨率
scRNA <- FindNeighbors(scRNA, reduction = "harmony",dims = 1:30, verbose = FALSE)
for (res in c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5,0.8,1)) {
  print(res)
  scRNA <- FindClusters(scRNA, graph.name = "RNA_snn", 
                        resolution = res, algorithm = 1)
}

DimPlot(scRNA,label = T)
colnames(scRNA@meta.data)
DimPlot(scRNA,label = T,group.by ="Class")
```



```bash

scRNA$anno=scRNA$Class

ss=scRNA
library(spacexr)

#prepare the reference
counts<-GetAssayData(ss,slot="counts")
cell_types<-ss$anno
cell_types<-as.factor(cell_types)
nUMI<-colSums(counts)
reference<-Reference(counts,cell_types,nUMI)

spatialobj=brain
#prepare the spatial data
counts<-GetAssayData(spatialobj,assay="Spatial",slot="counts")
coords<-GetTissueCoordinates(spatialobj)
nUMI<-colSums(counts)
puck<-SpatialRNA(coords,counts,nUMI)

#cell type deconvolution
myRCTD<-create.RCTD(puck,reference,max_cores=8)

#The doublet_mode argument sets whether RCTD will be run in ‘doublet mode’ (at most 1-2  cell types per pixel),
#‘full mode’ (no restrictions on number of cell  types), or
#‘multi mode’ (finitely many cell types per pixel, e.g. 3 or  4).
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')
```



```bash
#save the deconvolution result to the spatialobj
spatialobj@misc$RCTD<-myRCTD
spatialobj@misc$weights<-myRCTD@results$weights

# 反卷积结果
anno=myRCTD@results$weights
anno=as.data.frame(anno)
spatialobj=AddMetaData(spatialobj,metadata = anno)
colnames(spatialobj@meta.data)
```



```bash
SpatialFeaturePlot(spatialobj, features = c("Oligos","Astrocytes"), pt.size.factor = 1.6, ncol = 2, crop = TRUE)

```



```bash

myRCTD.d <- run.RCTD(myRCTD, doublet_mode = "doublet")

# 分配细胞注释
results_df=myRCTD.d@results[["results_df"]]

spatialobj=AddMetaData(spatialobj,metadata =results_df)
colnames(spatialobj@meta.data)

SpatialFeaturePlot(spatialobj, features = c("first_type","second_type"), pt.size.factor = 1.6, ncol = 2, crop = TRUE)
SpatialDimPlot(spatialobj,group.by = "second_type" )
```

```bash

barcodes <- colnames(myRCTD@spatialRNA@counts)
weights <- myRCTD@results$weights
norm_weights <- normalize_weights(weights)
# observe weight values
head(norm_weights)
```



```bash

p <- plot_puck_continuous(myRCTD@spatialRNA, barcodes, norm_weights[,"Oligos"], ylimit = c(0,0.5), 
                          title ='plot of Dentate weights', size=1.3, alpha=0.8) 
p
ggsave(  "Spaital_weights.png" , width=8, height=6, plot=p,bg="white")
```



```bash

#options(BioC_mirror="https://mirrors.westlake.edu.cn/bioconductor")
#BiocManager::install("STdeconvolve")
library(STdeconvolve)
library(ggplot2)
library(ggsci)
packageVersion("STdeconvolve")

m <- as.matrix(norm_weights)
p <- coords
colnames(p)=c("x","y")
plt <- vizAllTopics(theta = m,
                    pos = p,
                    topicOrder=seq(ncol(m)),
                    topicCols=rainbow(ncol(m)),
                    groups = NA,
                    group_cols = NA,
                    r = 3, # size of scatterpies; adjust depending on the coordinates of the pixels
                    lwd = 0.3,
                    showLegend = TRUE,
                    plotTitle = "scatterpies")

## function returns a `ggplot2` object, so other aesthetics can be added on:
plt <- plt + ggplot2::guides(fill=ggplot2::guide_legend(ncol=2))
plt
ggsave( "Spaital_scatterpies.png" , width=12, height=6, plot=plt, bg="white")

```



```bash


##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())


library(spacexr)
library(Matrix)
library(doParallel)

library(ggplot2)

setwd("Doublet_mode/")

### Load in/preprocess your data, this might vary based on your file type
# directory for the reference
refdir <- system.file("extdata",'Reference/Vignette',package = 'spacexr') 
# load in counts matrix
counts <- read.csv(file.path(refdir,"dge.csv")) 
# Move first column to rownames
rownames(counts) <- counts[,1]
counts[,1] <- NULL
counts[1:4,1:4]


# load in meta_data (barcodes, clusters, and nUMI)
meta_data <- read.csv(file.path(refdir,"meta_data.csv")) 
head(meta_data)

# create cell_types named list
cell_types <- meta_data$cluster
names(cell_types) <- meta_data$barcode 
# convert to factor data type
cell_types <- as.factor(cell_types) 
head(cell_types)

# create nUMI named list
nUMI <- meta_data$nUMI
names(nUMI) <- meta_data$barcode 
head(nUMI)

### Create the Reference object
reference <- Reference(counts, cell_types, nUMI)
str(reference)

## Examine reference object (optional)
#observe Digital Gene Expression matrix
print(dim(reference@counts)) 
#> [1] 384 475

#number of occurences for each cell type
table(reference@cell_types) 

## Save RDS object (optional)
saveRDS(reference, file.path(refdir,'SCRef.rds'))




##  读取空间数据
# directory for sample Slide-seq dataset
datadir <- system.file("extdata",'SpatialRNA/Vignette',package = 'spacexr') 
counts <- read.csv(file.path(datadir,"MappedDGEForR.csv"))
rownames(counts) <- counts[,1]; counts[,1] <- NULL
counts[1:4,1:4]

coords <- read.csv(file.path(datadir,"BeadLocationsForR.csv"))
rownames(coords) <- coords$barcodes
coords$barcodes <- NULL
head(coords)

# In this case, total counts per pixel is nUMI
nUMI <- colSums(counts) 
head(nUMI)

### Create SpatialRNA object
puck <- SpatialRNA(coords, counts, nUMI)
str(puck)


### Create SpatialRNA object
puck <- SpatialRNA(coords, counts, nUMI)
str(puck)

## Examine SpatialRNA object (optional)
print(dim(puck@counts))

# histogram of log_2 nUMI
hist(log(puck@nUMI,2))

print(head(puck@coords)) # start of coordinate data.frame
barcodes <- colnames(puck@counts) # pixels to be used (a list of barcode names). 

# This list can be restricted if you want to crop the puck e.g. 
# puck <- restrict_puck(puck, barcodes) provides a basic plot of the nUMI of each pixel
# on the plot:
p <- plot_puck_continuous(puck, barcodes, puck@nUMI, ylimit = c(0,round(quantile(puck@nUMI,0.9))), title ='plot of nUMI') 
p

### 三.创建RCTD对象和运行RCTD
myRCTD <- create.RCTD(puck, reference, max_cores = 10)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet')
str(myRCTD)

 

results <- myRCTD@results

# normalize the cell type proportions to sum to 1.
norm_weights = normalize_weights(results$weights) 

#list of cell type names
cell_type_names <- myRCTD@cell_type_info$info[[2]] 
spatialRNA <- myRCTD@spatialRNA

## you may change this to a more accessible directory on your computer.
resultsdir <- 'RCTD_Plots' 
dir.create(resultsdir)
 

# make the plots 
# 绘制full_mode模式下每种细胞类型的可信权重 (saved as 
# 'results/cell_type_weights_unthreshold.pdf')
plot_weights(cell_type_names, spatialRNA, resultsdir, norm_weights) 

# 绘制full_mode模式下每种细胞类型的权重
# 这里每种细胞类型一幅图，点表示空间上的一个像素或者spot，颜色为权重 (saved as
# 'results/cell_type_weights.pdf')
plot_weights_unthreshold(cell_type_names, spatialRNA, resultsdir, norm_weights) 

# 绘制full_mode模式下每种细胞类型预测到的spots数 (saved as 
# 'results/cell_type_occur.pdf')
plot_cond_occur(cell_type_names, resultsdir, norm_weights, spatialRNA)

# 绘制doublet_mode模式下每种细胞类型的权重 (saved as 
# 'results/cell_type_weights_doublets.pdf')
plot_weights_doublet(cell_type_names, spatialRNA, resultsdir, results$weights_doublet, results$results_df) 

```



```bash


##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())

library(Giotto)
library(Seurat)
library(tidyverse)
library(patchwork)
setwd("NC/")

data <- Load10X_Spatial("10x_brain/")

brain=data
##Data preprocessing
plot1 <- VlnPlot(brain, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(brain, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)
#Gene expression visualization

SpatialFeaturePlot(brain, features = c("Sox2", "Gapdh"))

library(ggplot2)
plot <- SpatialFeaturePlot(brain, features = c("Sox2")) + theme(legend.text = element_text(size = 0),
                                                                legend.title = element_text(size = 20), legend.key.size = unit(1, "cm"))
jpeg(filename = "spatial_vignette_ttr.jpg", height = 700, width = 1200, quality = 50)
print(plot)
dev.off()


p1 <- SpatialFeaturePlot(brain, features = "Sox2", pt.size.factor = 1)
p2 <- SpatialFeaturePlot(brain, features ="Sox2", alpha = c(0.1, 1))
p1 + p2

##Dimensionality reduction, clustering, and visualization

brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)
brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
brain <- FindClusters(brain, verbose = FALSE)
brain <- RunUMAP(brain, reduction = "pca", dims = 1:30)


p1 <- DimPlot(brain, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(brain, label = TRUE, label.size = 3)
p1 + p2

 
library(data.table)
exp=fread("brain_sc_expression_matrix.txt",  data.table = F)
rownames(exp)=exp[,1]
exp=exp[,-1]

meta=read.csv("brain_sc_metadata.csv")
meta=meta[,c(1,2,3,4,12)] 
rownames(meta)=meta$X

sc=CreateSeuratObject(exp,meta.data = meta)

table(meta$Class)
colnames(meta)
scRNA=sc
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^mt-")
# 计算细胞中核糖体基因比例
#scRNA[["percent.rb"]] <- PercentageFeatureSet(scRNA, pattern = "^RP[LS]")
# 质控
scRNA <- subset(scRNA, nCount_RNA>1000&nFeature_RNA>200&percent.mt<25)
# 标准化
scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)
# 鉴定高变基因
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 3000)
# 尺度变换
scRNA <- ScaleData(scRNA, features = rownames(scRNA))
# 降维
library(harmony)
scRNA <- RunPCA(scRNA, features = VariableFeatures(object = scRNA))
scRNA <- RunHarmony(scRNA, group.by.vars="orig.ident", max.iter.harmony = 20)
scRNA <- RunUMAP(scRNA,reduction = "harmony", dims = 1:30, verbose = FALSE)
## 确定分辨率
scRNA <- FindNeighbors(scRNA, reduction = "harmony",dims = 1:30, verbose = FALSE)
for (res in c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5,0.8,1)) {
  print(res)
  scRNA <- FindClusters(scRNA, graph.name = "RNA_snn", 
                        resolution = res, algorithm = 1)
}
 

DimPlot(scRNA,label = T)
colnames(scRNA@meta.data)
DimPlot(scRNA,label = T,group.by ="Class")





scRNA$anno=scRNA$Class

ss=scRNA
library(spacexr)

#prepare the reference
counts<-GetAssayData(ss,slot="counts")
cell_types<-ss$anno
cell_types<-as.factor(cell_types)
nUMI<-colSums(counts)
reference<-Reference(counts,cell_types,nUMI)

spatialobj=brain
#prepare the spatial data
counts<-GetAssayData(spatialobj,assay="Spatial",slot="counts")
coords<-GetTissueCoordinates(spatialobj)
nUMI<-colSums(counts)
puck<-SpatialRNA(coords,counts,nUMI)

#cell type deconvolution
myRCTD<-create.RCTD(puck,reference,max_cores=16)

#The doublet_mode argument sets whether RCTD will be run in ‘doublet mode’ (at most 1-2  cell types per pixel),
#‘full mode’ (no restrictions on number of cell  types), or
#‘multi mode’ (finitely many cell types per pixel, e.g. 3 or  4).
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')

#save the deconvolution result to the spatialobj
spatialobj@misc$RCTD<-myRCTD
spatialobj@misc$weights<-myRCTD@results$weights


# 反卷积结果
anno=myRCTD@results$weights
anno=as.data.frame(anno)
spatialobj=AddMetaData(spatialobj,metadata = anno)
colnames(spatialobj@meta.data)

SpatialFeaturePlot(spatialobj, features = c("Oligos","Astrocytes"), pt.size.factor = 1.6, ncol = 2, crop = TRUE)
 

?run.RCTD
myRCTD.d <- run.RCTD(myRCTD, doublet_mode = "doublet")

# 分配细胞注释
results_df=myRCTD.d@results[["results_df"]]

spatialobj=AddMetaData(spatialobj,metadata =results_df)
colnames(spatialobj@meta.data)

SpatialFeaturePlot(spatialobj, features = c("first_type","second_type"), pt.size.factor = 1.6, ncol = 2, crop = TRUE)
SpatialDimPlot(spatialobj,group.by = "second_type" )

 #### 进一步探索

 
barcodes <- colnames(myRCTD@spatialRNA@counts)
weights <- myRCTD@results$weights
norm_weights <- normalize_weights(weights)
# observe weight values
head(norm_weights)

A# plot Dentate weights
p <- plot_puck_continuous(myRCTD@spatialRNA, barcodes, norm_weights[,"Oligos"], ylimit = c(0,0.5), 
                          title ='plot of Dentate weights', size=1.3, alpha=0.8) 
p
ggsave(  "Spaital_weights.png" , width=8, height=6, plot=p,bg="white")

#options(BioC_mirror="https://mirrors.westlake.edu.cn/bioconductor")
#BiocManager::install("STdeconvolve")
library(STdeconvolve)
library(ggplot2)
library(ggsci)
packageVersion("STdeconvolve")

m <- as.matrix(norm_weights)
p <- coords
colnames(p)=c("x","y")
plt <- vizAllTopics(theta = m,
                    pos = p,
                    topicOrder=seq(ncol(m)),
                    topicCols=rainbow(ncol(m)),
                    groups = NA,
                    group_cols = NA,
                    r = 3, # size of scatterpies; adjust depending on the coordinates of the pixels
                    lwd = 0.3,
                    showLegend = TRUE,
                    plotTitle = "scatterpies")

## function returns a `ggplot2` object, so other aesthetics can be added on:
plt <- plt + ggplot2::guides(fill=ggplot2::guide_legend(ncol=2))
plt
ggsave( "Spaital_scatterpies.png" , width=12, height=6, plot=plt, bg="white")
```


