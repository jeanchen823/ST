# CellChat

```bash
## 需要先安装依赖包##
install.packages("devtools")
BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
BiocManager::install("BiocNeighbors")
library(BiocNeighbors)
devtools::install_github("jinworks/CellChat")
library(CellChat)
```

```bash
##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())

library(CellChat)
library(patchwork)
setwd("cellchat/")
load("visium_mouse_cortex_annotated.RData")
```

```bash
visium.brain
```

```bash
SpatialDimPlot(visium.brain, label = T, pt.size.factor = 1000,label.size = 4 )
```

```bash
data.input = GetAssayData(visium.brain, slot = "data", assay = "SCT") # normalized data matrix
```

```bash
meta = data.frame(labels = Idents(visium.brain), row.names = names(Idents(visium.brain))) # manually create a dataframe consisting of the cell labels
# check the cell labels
unique(meta$labels)
```



```bash
# load spatial imaging information
# Spatial locations of spots from full (NOT high/low) resolution images are required
spatial.locs = GetTissueCoordinates(visium.brain, scale = NULL, cols = c("imagerow", "imagecol"))
```



```bash
scale.factors = jsonlite::fromJSON(txt = file.path("D:/spe.lesson/k15/cellchat/image/", 'scalefactors_json.json'))
spot.size = 65 #10X Visium spot大小为55μm，两个spot之间Gap为10μm
conversion.factor = spot.size/scalefactors$spot_diameter_fullres
spatial.factors = data.frame(ratio = conversion.factor, tol = spot.size/2)
d.spatial <- computeCellDistance(coordinates = spatial.locs, ratio = spatial.factors$ratio, tol = spatial.factors$tol)
```



```bash
cellchat <- createCellChat(object = data.input,                           
							meta = meta,                           
							group.by = "labels", #定义的名字是labels
							datatype = "spatial", #数据类型：空转
							coordinates = spatial.locs,
							spatial.factors = spatial.factors)
```



```bash
CellChatDB <- CellChatDB.mouse# use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
```

```bash
#使用CellChatDB的子集进行细胞间通信分析
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") #选择Secreted Signaling
cellchat@DB <- CellChatDB.use
```



```bash
https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/Update-CellChatDB.html
```



```bash
# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
```



```bash
future::plan("multisession", workers = 4) #多线程
```



```bash
#识别过表达基因
#devtools::install_github('immunogenomics/presto')
cellchat <- identifyOverExpressedGenes(cellchat)
```



```bash
#识别过表达配体受体对
cellchat <- identifyOverExpressedInteractions(cellchat, variable.both = F)
```

```bash
cellchat <- computeCommunProb(cellchat,
type = "truncatedMean", trim = 0.1,
distance.use = TRUE, contact.range=10,
scale.distance = 0.01)
```



```bash
#细胞间通信网络的推断
cellchat <- computeCommunProb(cellchat,                              
							type = "truncatedMean", trim = 0.1,
							distance.use = TRUE, contact.range =10, 
							scale.distance = 0.01)#默认情况下，每个细胞组中用于细胞间通信所需的最小细胞数为10
cellchat <- filterCommunication(cellchat, min.cells = 10)

#在信号通路水平上推断细胞间通讯
cellchat <- computeCommunProbPathway(cellchat)
#计算聚合的 cell-cell 通信网络
cellchat <- aggregateNet(cellchat)
```



```bash
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = rowSums(cellchat@net$count), weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = rowSums(cellchat@net$weight), weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")


#热图显示celltype间的通讯次数（左）或总通讯强度(右)
p1 <- netVisual_heatmap(cellchat, measure = "count", color.heatmap = "Blues")
p2 <- netVisual_heatmap(cellchat, measure = "weight", color.heatmap = "Blues")
p1 + p2


cellchat@netP$pathways


par(mfrow=c(1,1), xpd = TRUE)# xpd = TRUE以显示标题
pathways.show <- c("PTN")
#可视化 'PTN' 信号网络
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")



netVisual_aggregate(cellchat,
					signaling = pathways.show,
					layout = "spatial",
					edge.width.max = 2,
					vertex.size.max = 1,
					alpha.image = 0.2,
					vertex.label.cex = 3.5)


par(mfrow=c(1,1))

netVisual_aggregate(cellchat,                    
signaling = pathways.show,                    
layout = "spatial",                    
edge.width.max = 2,                    
alpha.image = 0.2,                    
vertex.weight = "incoming",                    
vertex.size.max = 4,                    
vertex.label.cex = 3.5)


cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")#“netP”是指推断的信号通路的细胞间通信网络
par(mfrow=c(1,1))
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)


spatialFeaturePlot(cellchat,                   
features = c("Mdk","Sdc1"),                   
point.size = 0.8,                   
color.heatmap = "Reds",                   
direction = 1)


spatialFeaturePlot(cellchat, pairLR.use = "MDK_SDC1", point.size = 3, do.binary = TRUE, cutoff = 0.05, enriched.only = F, color.heatmap = "Reds", direction = 1)


netVisual_bubble(cellchat, sources.use = c(3),                 signaling=cellchat@netP$pathways[1:6],                 
targets.use = c(1,2),                 
remove.isolate = FALSE)


netVisual_chord_gene(cellchat, sources.use = 3,                     
targets.use = c(1:2),                     
signaling=cellchat@netP$pathways[1:6],                     
lab.cex = 0.5,                     
legend.pos.y = 30)


p = plotGeneExpression(cellchat, signaling = "PTN")p

```



```bash
## install.packages("devtools")
BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
BiocManager::install("BiocNeighbors")
library(BiocNeighbors)
devtools::install_github("jinworks/CellChat")
library(CellChat)



##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())

library(CellChat)
library(patchwork)
library(Seurat)
load("visium_mouse_cortex_annotated.RData")

visium.brain
SpatialDimPlot(visium.brain, label = T, pt.size.factor = 1000,label.size = 4 )

# Prepare input data for CelChat analysis
data.input = GetAssayData(visium.brain, slot = "data", assay = "SCT") # normalized data matrix
meta = data.frame(labels = Idents(visium.brain), row.names = names(Idents(visium.brain))) # manually create a dataframe consisting of the cell labels
unique(meta$labels) # check the cell labels
#获取meta信息
meta = data.frame(labels = Idents(visium.brain),
row.names = names(Idents(visium.brain)))
unique(meta$labels)

#获取空间位置信息
spatial.locs = Seurat::GetTissueCoordinates(visium.brain, scale = NULL,                                            cols = c("imagerow", "imagecol"))

scalefactors = jsonlite::fromJSON(txt = file.path("cellchat/image/", 'scalefactors_json.json'))

spot.size = 65 #10X Visium spot大小为55μm，两个spot之间Gap为10μm

conversion.factor = spot.size/scalefactors$spot_diameter_fullresspatial.factors = data.frame(ratio = conversion.factor, tol = spot.size/2)

d.spatial <- computeCellDistance(coordinates = spatial.locs, ratio = spatial.factors$ratio, tol = spatial.factors$tol)


#创建CellChat对象
cellchat <- createCellChat(object = data.input,                           
							meta = meta,                           
							group.by = "labels", #定义的名字是labels
							datatype = "spatial", #数据类型：空转
							coordinates = spatial.locs,
							spatial.factors = spatial.factors)


#设置参考数据库
CellChatDB <- CellChatDB.mouse# use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

#使用CellChatDB的子集进行细胞间通信分析
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") #选择Secreted Signaling
cellchat@DB <- CellChatDB.use


#CellChat预处理
?subset
Datacellchat <- subsetData(cellchat, features = NULL) #即使使用整个数据库，此步骤也是必要的
future::plan("multisession", workers = 4) #多线程

#识别过表达基因
#devtools::install_github('immunogenomics/presto')
cellchat <- identifyOverExpressedGenes(cellchat)

#识别过表达配体受体对
cellchat <- identifyOverExpressedInteractions(cellchat, variable.both = F)

#细胞间通信网络的推断
cellchat <- computeCommunProb(cellchat,
type = "truncatedMean", trim = 0.1,
distance.use = TRUE, contact.range=10,
scale.distance = 0.01)

#默认情况下，每个细胞组中用于细胞间通信所需的最小细胞数为10
cellchat <- filterCommunication(cellchat, min.cells = 10)

#在信号通路水平上推断细胞间通讯
cellchat <- computeCommunProbPathway(cellchat)
#计算聚合的 cell-cell 通信网络
cellchat <- aggregateNet(cellchat)
#可视化交互次数或总交互次数
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = rowSums(cellchat@net$count), weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = rowSums(cellchat@net$weight), weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

#热图显示celltype间的通讯次数（左）或总通讯强度(右)
p1 <- netVisual_heatmap(cellchat, measure = "count", color.heatmap = "Blues")
p2 <- netVisual_heatmap(cellchat, measure = "weight", color.heatmap = "Blues")
p1 + p2



#展示显著通路结果
cellchat@netP$pathways
par(mfrow=c(1,1), xpd = TRUE)# xpd = TRUE以显示标题
pathways.show <- c("PTN")
#可视化 'PTN' 信号网络
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")


#在空间转录组上显示'PTN'信号网络
netVisual_aggregate(cellchat,
					signaling = pathways.show,
					layout = "spatial",
					edge.width.max = 2,
					vertex.size.max = 1,
					alpha.image = 0.2,
					vertex.label.cex = 3.5)



#在空间转录组学上显示信号网络时，可以以更大的圆圈表示更大的传入信号
par(mfrow=c(1,1))

netVisual_aggregate(cellchat,                    
signaling = pathways.show,                    
layout = "spatial",                    
edge.width.max = 2,                    
alpha.image = 0.2,                    
vertex.weight = "incoming",                    
vertex.size.max = 4,                    
vertex.label.cex = 3.5)



#计算和可视化网络中心性分数：
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")#“netP”是指推断的信号通路的细胞间通信网络
par(mfrow=c(1,1))
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)





#可视化组织上的基因表达空间分布
spatialFeaturePlot(cellchat,                   
features = c("Mdk","Sdc1"),                   
point.size = 0.8,                   
color.heatmap = "Reds",                   
direction = 1)

#取配体-受体对的输入，并以二进制形式显示表达
spatialFeaturePlot(cellchat, pairLR.use = "MDK_SDC1", point.size = 3, do.binary = TRUE, cutoff = 0.05, enriched.only = F, color.heatmap = "Reds", direction = 1)

#绘制配体-受体气泡图
netVisual_bubble(cellchat, sources.use = c(3),                 signaling=cellchat@netP$pathways[1:6],                 
targets.use = c(1,2),                 
remove.isolate = FALSE)


#显示前6条信号通路重的（L-R 对）:
netVisual_chord_gene(cellchat, sources.use = 3,                     
targets.use = c(1:2),                     
signaling=cellchat@netP$pathways[1:6],                     
lab.cex = 0.5,                     
legend.pos.y = 30)


#信号通路的所有基因在细胞群中的表达情况展示
p = plotGeneExpression(cellchat, signaling = "PTN")
p

```


