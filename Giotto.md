# Giotto 

```bash

https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02286-2

01
 
安装 Giotto

跟着官网学安装：

https://giottosuite.readthedocs.io/en/master/gettingstarted.html

Python Giotto Requirements
##To perform all potential steps and analysis in the Giotto spatial toolbox the user needs to have a number of python modules installed. To make this process as flexible and easy as possible two different strategies can be used
pandas
python-igraph (igraph)
networkx
leidenalg
python-louvain (community)
smfishHmrf
python.app (OSX only)
scikit-learn
##These are necessary to run all available analyses, but can be installed automatically 


1.创建pyhton环境
conda create -n giotto_env   python=3.10.2  -c https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge
conda activate giotto_env

2.安装依赖包

2.1 pip3安装---指定版本

https://drieslab.github.io/Giotto_website/articles/installation.html

版本指定
https://drieslab.github.io/GiottoClass/reference/giotto_python.html


pip3 install pandas==1.5.1




name: giotto_env
channels:
— defaults
dependencies:
— pip=3.4
— pandas
— networkx
- python-igraph
- leidenalg
- python-louvain
- python.app (!!only for OSX!!)
- scikit-learn
prefix: /Users/your_username/anaconda3/envs/giotto_env




library(devtools) # If not installed: install.packages('devtools')
library(remotes)  #If not installed: install.packages('remotes')

options(BioC_mirror="https://mirrors.westlake.edu.cn/bioconductor")
remotes::install_github("drieslab/Giotto@master")




https://github.com/drieslab/Giotto



remotes::install_local("Giotto-suite.zip",upgrade = F,dependencies = T)


library(Giotto)



checkGiottoEnvironment()


```

```bash
02
 
 创建Giotto Visium对象和可视化
 
 https://drieslab.github.io/Giotto_website/articles/installation.html
 
 
 
 1.加载R包、设置路径

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

2.可选：在conda或miniconda中指定Python可执行文件的路径

环境：如果设置为NULL(默认值)，则先前的Python可执行文件

或者将使用已安装的Giotto环境

python_path <- NULL 
# alternatively, "/local/python/path/python" if desired.
#python_path <- "D:\\miniconda\\envs\\squidpy/python.exe"
3.Create Giotto Instructions
instructions <- createGiottoInstructions(save_dir =  '10x_brain/out/',
                                         save_plot = TRUE,
                                         show_plot = FALSE,
                                         return_plot = FALSE,
                                         python_path = python_path)
                                         
           
4.提供visium数据所在的文件夹的路径并读取数据

data_path <- "10x_brain/"

## directly from visium folder
visium_brain <- createGiottoVisiumObject(visium_dir = data_path,
                                         expr_data = "filter",
                                         png_name = "tissue_lowres_image.png",
                                         gene_column_index = 2,
                                         instructions = instructions)

5.show associated images with giotto object
showGiottoImageNames(visium_brain) # "image" is the default name                               
6.check metadata

pDataDT(visium_brain)


7.画图

spatPlot2D(gobject = visium_brain, 
           cell_color = "in_tissue", 
           point_size = 2,
           cell_color_code = c("0" = "lightgrey", "1" = "blue"), 
           show_image = TRUE, 
           image_name = "image")
          

 
```



```bash
03
 
 数据处理

1.被组织覆盖的 spot 的子集
metadata <- pDataDT(visium_brain)
in_tissue_barcodes <- metadata[in_tissue == 1]$cell_ID

visium_brain <- subsetGiotto(visium_brain, 
                             cell_ids = in_tissue_barcodes)
2.过滤
visium_brain <- filterGiotto(gobject = visium_brain,
                             expression_threshold = 1,
                             feat_det_in_min_cells = 50,
                             min_det_feats_per_cell = 1000,
                             expression_values = "raw",
                             verbose = TRUE)

3.归一化

visium_brain <- normalizeGiotto(gobject = visium_brain, 
                                scalefactor = 6000, 
                                verbose = TRUE)
4.添加基因和细胞的统计信息
visium_brain <- addStatistics(gobject = visium_brain)

## visualize 画图
spatPlot2D(gobject = visium_brain, 
           show_image = TRUE, 
           point_alpha = 0.7,
           cell_color = "nr_feats", 
           color_as_factor = FALSE)
```

```bash
04
 
 降维
 
 
1.找高变基因

visium_brain <- calculateHVF(gobject = visium_brain, 
                             save_plot = TRUE)
                             
                             
2.PCA降维

gene_metadata <- fDataDT(visium_brain)
featgenes <- gene_metadata[hvf == "yes" & perc_cells > 3 & mean_expr_det > 0.4]$feat_ID

## run PCA on expression values (default)
visium_brain <- runPCA(gobject = visium_brain,
                       feats_to_use = featgenes)

screePlot(visium_brain, 
          ncp = 30)


dimPlot2D(gobject = visium_brain,
          dim_reduction_to_use = "pca")
          
3.Umap 和tsne 可视化降维

visium_brain <- runUMAP(visium_brain, 
                        dimensions_to_use = 1:10)
plotUMAP(gobject = visium_brain)

visium_brain <- runtSNE(visium_brain, 
                        dimensions_to_use = 1:10)
plotTSNE(gobject = visium_brain)


```

```bash
05
 
 聚类

1. sNN network

visium_brain <- createNearestNetwork(gobject = visium_brain, 
                                     dimensions_to_use = 1:10, 
                                     k = 15)
2.Leiden 聚类
visium_brain <- doLeidenCluster(gobject = visium_brain, 
                                resolution = 0.4, 
                                n_iterations = 1000)
## 可视化
3.画图展示

spatDimPlot(gobject = visium_brain, 
            cell_color = "leiden_clus",
            dim_point_size = 2, 
            spat_point_size = 2.5)

spatDimPlot(gobject = visium_brain, 
            cell_color = "nr_feats", 
            color_as_factor = FALSE,
            dim_point_size = 2, 
            spat_point_size = 2.5)

# dimension plots grouped by cluster
spatPlot2D(visium_brain, 
           cell_color = "leiden_clus",
           coord_fix_ratio = 1)

spatPlot2D(visium_brain, 
           cell_color = "leiden_clus",
           select_cell_groups = "8", 
           coord_fix_ratio = 1, 
           show_other_cells = TRUE,
           cell_color_code = c("8" = "red"), 
           other_cell_color = "grey", 
           other_point_size = 1.5)
           
Plot with group by:

spatPlot2D(visium_brain, 
           cell_color = "leiden_clus",
           group_by = "leiden_clus", 
           coord_fix_ratio = 1,
           cow_n_col = 6, 
           show_legend = FALSE)

```

```bash
06 
 数据提取

1.create subset

DG_subset <- subsetGiottoLocs(visium_brain,
                             x_max = 6500, x_min = 3000,
                             y_max = -2500, y_min = -5500,
                             return_gobject = TRUE)
2. show subset
spatDimPlot(gobject = DG_subset,
            cell_color = "leiden_clus", 
            spat_point_size = 5)

```



```bash
07
 
 找 marker 基因

1.Gini markers

markers_gini <- findMarkers_one_vs_all(gobject = visium_brain,
                                       method = "gini",
                                       expression_values = "normalized",
                                       cluster_column = "leiden_clus",
                                       min_feats = 20,
                                       min_expr_gini_score = 0.5,
                                       min_det_gini_score = 0.5)

topgenes_gini <- markers_gini[, head(.SD, 2), by = "cluster"]$feats

2.画热图和小提琴图展示

# violinplot
violinPlot(visium_brain, 
           feats = unique(topgenes_gini), 
           cluster_column = "leiden_clus",
           strip_text = 8, 
           strip_position = "right")

# cluster heatmap
plotMetaDataHeatmap(visium_brain, 
                    selected_feats = unique(topgenes_gini),
                    metadata_cols = "leiden_clus",
                    x_text_size = 10, 
                    y_text_size = 10)
                    
dimFeatPlot2D(visium_brain, 
              expression_values = "scaled",
              feats = markers_gini[, head(.SD, 1), by = "cluster"]$feats,
              cow_n_col = 4, 
              point_size = 0.75)

第二种算法找marker基因---Scran Markers


markers_scran <- findMarkers_one_vs_all(gobject = visium_brain,
                                        method = "scran",
                                        expression_values = "normalized",
                                        cluster_column = "leiden_clus")

topgenes_scran <- markers_scran[, head(.SD, 2), by = "cluster"]$feats
画图展示

# violinplot
violinPlot(visium_brain, 
           feats = unique(topgenes_scran), 
           cluster_column = "leiden_clus",
           strip_text = 10, 
           strip_position = "right")
# cluster heatmap
plotMetaDataHeatmap(visium_brain, 
                    selected_feats = topgenes_scran,
                    metadata_cols = "leiden_clus")

dimFeatPlot2D(visium_brain, 
              expression_values = "scaled",
              feats = markers_scran[, head(.SD, 1), by = "cluster"]$feats,
              cow_n_col = 3, 
              point_size = 1)



```



```bash
08
 
 细胞类型注释

1.创建单细胞数据对象
sc_expression <- file.path( "10x_brain/out/brain_sc_expression_matrix.txt.gz")
sc_metadata <- file.path( "10x_brain/out/brain_sc_metadata.csv")

giotto_SC <- createGiottoObject(expression = sc_expression,
                                instructions = instructions)

giotto_SC <- addCellMetadata(giotto_SC, 
                             new_metadata = data.table::fread(sc_metadata))

giotto_SC <- normalizeGiotto(giotto_SC)

1.1 创建 PAGE 矩阵

PAGE矩阵应为二进制矩阵，每行代表一个基因标记，每列代表一个细胞类型

有几种方法可以创建PAGE矩阵

1.1.1 用细胞signature genes创建二进制矩阵

gran_markers <- c("Nr3c2", "Gabra5", "Tubgcp2", "Ahcyl2",
                  "Islr2", "Rasl10a", "Tmem114", "Bhlhe22", 
                  "Ntf3", "C1ql2")

oligo_markers <- c("Efhd1", "H2-Ab1", "Enpp6", "Ninj2",
                   "Bmp4", "Tnr", "Hapln2", "Neu4",
                   "Wfdc18", "Ccp110")        

di_mesench_markers <- c("Cartpt", "Scn1a", "Lypd6b",  "Drd5",
                        "Gpr88", "Plcxd2", "Cpne7", "Pou4f1",
                        "Ctxn2", "Wnt4")

PAGE_matrix_1 <- makeSignMatrixPAGE(sign_names = c("Granule_neurons",
                                                   "Oligo_dendrocytes",
                                                   "di_mesenchephalon"),
                                    sign_list = list(gran_markers,
                                                     oligo_markers,
                                                     di_mesench_markers))
1.1.2 ：fully pre-prepared matrix for all cell types

sign_matrix_path <- system.file("extdata", "sig_matrix.txt", package = "GiottoData")

brain_sc_markers <- data.table::fread(sign_matrix_path) 

PAGE_matrix <- as.matrix(brain_sc_markers[,-1])
rownames(PAGE_matrix) <- brain_sc_markers$Event


1.1.3 ：make PAGE matrix from single cell dataset

markers_scran <- findMarkers_one_vs_all(gobject = giotto_SC, 
                                        method = "scran",
                                        expression_values = "normalized",
                                        cluster_column = "Class", 
                                        min_feats = 3)

topgenes_scran <- markers_scran[, head(.SD, 10), by = "cluster"]

celltypes <- levels(factor(markers_scran$cluster)) 

sign_list <- list()

for (i in 1:length(celltypes)){
  sign_list[[i]] <- topgenes_scran[which(topgenes_scran$cluster == celltypes[i]),]$feats
}

PAGE_matrix <- makeSignMatrixPAGE(sign_names = celltypes,
                                  sign_list = sign_list)

1.2 ：enrichment test with PAGE
runSpatialEnrich() can also be used as a wrapper for all currently provided enrichment options
visium_brain <- runPAGEEnrich(gobject = visium_brain, 
                              sign_matrix = PAGE_matrix)
1.3 ：heatmap of enrichment versus annotation (e.g. clustering result)
cell_types_PAGE <- colnames(PAGE_matrix)

plotMetaDataCellsHeatmap(gobject = visium_brain,
                         metadata_cols = "leiden_clus",
                         value_cols = cell_types_PAGE,
                         spat_enr_names = "PAGE",
                         x_text_size = 8,
                         y_text_size = 8)


1.4 可视化注释结果

spatCellPlot2D(gobject = visium_brain,
               spat_enr_names = "PAGE",
               cell_annotation_values = cell_types_PAGE[1:4],
               cow_n_col = 2,
               coord_fix_ratio = 1, 
               point_size = 1.25, 
               show_legend = TRUE)
               

spatDimCellPlot2D(gobject = visium_brain,
                  spat_enr_names = "PAGE",
                  cell_annotation_values = cell_types_PAGE[1:4],
                  cow_n_col = 1, 
                  spat_point_size = 1,
                  plot_alignment = "horizontal",
                  save_param = list(base_width = 7, base_height = 10))
                  
       
       
       

visium_brain <- runHyperGeometricEnrich(gobject = visium_brain,
                                        expression_values = "normalized",
                                        sign_matrix = PAGE_matrix)

cell_types_HyperGeometric <- colnames(PAGE_matrix)

spatCellPlot(gobject = visium_brain,
             spat_enr_names = "hypergeometric",
             cell_annotation_values = cell_types_HyperGeometric[1:4],
             cow_n_col = 2,
             coord_fix_ratio = NULL, 
             point_size = 1.75)




rank_matrix <- makeSignMatrixRank(sc_matrix = getExpression(giotto_SC,
                                                            values = "normalized",
                                                            output = "matrix"),
                                  sc_cluster_ids = pDataDT(giotto_SC)$Class)

colnames(rank_matrix) <- levels(factor(pDataDT(giotto_SC)$Class))

visium_brain <- runRankEnrich(gobject = visium_brain, 
                              sign_matrix = rank_matrix,
                              expression_values = "normalized")

画图

spatCellPlot2D(gobject = visium_brain,
               spat_enr_names = "rank",
               cell_annotation_values = colnames(rank_matrix)[1:4],
               cow_n_col = 2,
               coord_fix_ratio = 1, 
               point_size = 1)


方法四：DWLS 空间反卷积

DWLS_matrix <- makeSignMatrixDWLSfromMatrix(
  matrix = getExpression(giotto_SC,
                         values = "normalized",
                         output = "matrix"),
  cell_type = pDataDT(giotto_SC)$Class,
  sign_gene = topgenes_scran$feats)

visium_brain <- runDWLSDeconv(gobject = visium_brain, 
                              sign_matrix = DWLS_matrix)

画图展示反卷积的结果
# Plot DWLS deconvolution result
spatCellPlot2D(gobject = visium_brain,
               spat_enr_names = "DWLS",
               cell_annotation_values = levels(factor(pDataDT(giotto_SC)$Class))[1:4],
               cow_n_col = 2,
               coord_fix_ratio = 1, 
               point_size = 1)

Plot DWLS deconvolution result with Pie plots
spatDeconvPlot(visium_brain, 
               show_image = TRUE,
               radius = 50)    
               
09
 
 Spatial Grid--空间网格



visium_brain <- createSpatialGrid(gobject = visium_brain,
                                  sdimx_stepsize = 400,
                                  sdimy_stepsize = 400,
                                  minimum_padding = 0)

showGiottoSpatGrids(visium_brain)

spatPlot2D(visium_brain, 
           cell_color = "leiden_clus", 
           show_grid = TRUE,
           grid_color = "red", 
           spatial_grid_name = "spatial_grid")
           
           
           
10
 
 Spatial network---空间网络

创建一个空间网络，连接单个细胞基于它们彼此的物理距离。
visium_brain <- createSpatialNetwork(gobject = visium_brain,
                                     method = "kNN", 
                                     k = 5,
                                     maximum_distance_knn = 400,
                                     name = "spatial_network")

showGiottoSpatNetworks(visium_brain)

spatPlot2D(gobject = visium_brain,  
           show_network= TRUE,
           network_color = "blue", 
           spatial_network_name = "spatial_network")
        
        
        
11
 
 Spatial Genes---空间基因
 
ranktest <- binSpect(visium_brain, 
                     bin_method = "rank",
                     calc_hub = TRUE, 
                     hub_min_int = 5,
                     spatial_network_name = "spatial_network")

spatFeatPlot2D(visium_brain, 
               expression_values = "scaled",
               feats = ranktest$feats[1:6], 
               cow_n_col = 2, 
               point_size = 1.5)
12
 
 空间共表达模块



cluster the top 500 spatial genes into 20 cluster

my_spatial_genes <- ranktest[1:1500,]$feats
1.这里我们使用现有的detectSpatialCorGenes函数来计算基因之间的成对距离(但将network_smoothing=0设置为使用默认聚类)
spat_cor_netw_DT <- detectSpatialCorFeats(visium_brain,
                                          method = "network",
                                          spatial_network_name = "spatial_network",
                                          subset_feats = my_spatial_genes)

2.1 确定一个基因最相似的空间相关基因
top10_genes <- showSpatialCorFeats(spat_cor_netw_DT, 
                                   feats = "Mbp", 
                                   show_top_feats = 10)
2.2 可视化
spatFeatPlot2D(visium_brain, 
               expression_values = "scaled",
               feats = top10_genes$variable[1:4], 
               point_size = 3)


3.1cluster spatial genes
spat_cor_netw_DT <- clusterSpatialCorFeats(spat_cor_netw_DT, 
                                           name = "spat_netw_clus", 
                                           k = 20)
3.2 visualize clusters
heatmSpatialCorFeats(visium_brain,
                     spatCorObject = spat_cor_netw_DT,
                     use_clus_name = "spat_netw_clus",
                     heatmap_legend_param = list(title = NULL))



## 查看每个cluster的基因
netw_ranks = rankSpatialCorGroups(visium_brain, spatCorObject = spat_cor_netw_DT, 
                                  use_clus_name = 'spat_netw_clus')

cluster_genes_DT = showSpatialCorGenes(spat_cor_netw_DT, 
                                       use_clus_name = 'spat_netw_clus', 
                                       show_top_genes = 1)
cluster_genes = cluster_genes_DT$clus  
names(cluster_genes) = cluster_genes_DT$feat_ID

4. rank spatial correlated clusters and show genes for selected clusters

netw_ranks <- rankSpatialCorGroups(visium_brain,
                                   spatCorObject = spat_cor_netw_DT, 
                                   use_clus_name = "spat_netw_clus")

top_netw_spat_cluster <- showSpatialCorFeats(spat_cor_netw_DT, 
                                             use_clus_name = "spat_netw_clus",
                                             selected_clusters = 6, 
                                             show_top_feats = 1)


5. create metagene enrichment score for clusters
cluster_genes_DT <- showSpatialCorFeats(spat_cor_netw_DT, 
                                        use_clus_name = "spat_netw_clus", 
                                        show_top_feats = 1)

cluster_genes <- cluster_genes_DT$clus 
names(cluster_genes) <- cluster_genes_DT$feat_ID

visium_brain <- createMetafeats(visium_brain, 
                                feat_clusters = cluster_genes, 
                                name = "cluster_metagene")

spatCellPlot(visium_brain,
             spat_enr_names = "cluster_metagene",
             cell_annotation_values = netw_ranks$clusters,
             point_size = 1, 
             cow_n_col = 5)


13
 
 Spatially informed clusters
 
 
 
table(spat_cor_netw_DT$cor_clusters$spat_netw_clus)

coexpr_dt <- data.table::data.table(
  genes = names(spat_cor_netw_DT$cor_clusters$spat_netw_clus),
  cluster = spat_cor_netw_DT$cor_clusters$spat_netw_clus)

data.table::setorder(coexpr_dt, cluster)
top30_coexpr_dt <- coexpr_dt[, head(.SD, 30), by = cluster]

my_spatial_genes <- top30_coexpr_dt$genes

visium_brain <- runPCA(gobject = visium_brain,
                       feats_to_use = my_spatial_genes,
                       name = "custom_pca")

visium_brain <- runUMAP(visium_brain, 
                        dim_reduction_name = "custom_pca", 
                        dimensions_to_use = 1:20,
                        name = "custom_umap")

visium_brain <- createNearestNetwork(gobject = visium_brain,
                                     dim_reduction_name = "custom_pca",
                                     dimensions_to_use = 1:20, 
                                     k = 5,
                                     name = "custom_NN")

visium_brain <- doLeidenCluster(gobject = visium_brain, 
                                network_name = "custom_NN",
                                resolution = 0.15, 
                                n_iterations = 1000,
                                name = "custom_leiden")

cell_metadata <- pDataDT(visium_brain)
cell_clusters <- unique(cell_metadata$custom_leiden)

giotto_colors <- getDistinctColors(length(cell_clusters))
names(giotto_colors) <- cell_clusters

spatPlot2D(visium_brain, 
           cell_color = "custom_leiden", 
           cell_color_code = giotto_colors, 
           coord_fix_ratio = 1)

plotUMAP(gobject = visium_brain, 
         cell_color = "custom_leiden", 
         cell_color_code = giotto_colors, 
         point_size = 1.5)
         
4
 
 Spatial domains with HMRF

每个空间共表达模块前30个基因上具有不同β的HMRF
hmrf_folder <- file.path(data_path, "HMRF")

if(!file.exists(hmrf_folder)) dir.create(hmrf_folder, recursive = TRUE)
library(smfishHmrf)
HMRF_spatial_genes <- doHMRF(gobject = visium_brain,
                             expression_values = "scaled",
                             spatial_genes = my_spatial_genes, 
                             k = 20,
                             spatial_network_name="spatial_network",
                             betas = c(0, 10, 5),
                             output_folder = file.path(hmrf_folder, "Spatial_genes/SG_topgenes_k20_scaled"))

pip install smfishHmrf

再次运行

HMRF_spatial_genes <- doHMRF(gobject = visium_brain,
                             expression_values = "scaled",
                             spatial_genes = my_spatial_genes, 
                             k = 20,
                             spatial_network_name="spatial_network",
                             betas = c(0, 10, 5),
                             output_folder = file.path(hmrf_folder, "Spatial_genes/SG_topgenes_k20_scaled"))
                             


conda install conda-forge::openjdk


再次运行

HMRF_spatial_genes <- doHMRF(gobject = visium_brain,
                             expression_values = "scaled",
                             spatial_genes = my_spatial_genes, 
                             k = 20,
                             spatial_network_name="spatial_network",
                             betas = c(0, 10, 5),
                             output_folder = file.path(hmrf_folder, "Spatial_genes/SG_topgenes_k20_scaled"))

 
 
再次运行

HMRF_spatial_genes <- doHMRF(gobject = visium_brain,
                             expression_values = "scaled",
                             spatial_genes = my_spatial_genes, 
                             k = 20,
                             spatial_network_name="spatial_network",
                             betas = c(0, 10, 5),
                             output_folder = file.path(hmrf_folder, "Spatial_genes/SG_topgenes_k20_scaled"))



```



```bash




library(devtools) # If not installed: install.packages('devtools')
library(remotes)  #If not installed: install.packages('remotes')

options(BioC_mirror="https://mirrors.westlake.edu.cn/bioconductor")
remotes::install_github("drieslab/Giotto@master")

install.packages(c("smfishHmrf", "trendsceek", "multinet", "RTriangle", "FactoMineR"))

remotes::install_local("Giotto-suite.zip",upgrade = F,dependencies = T)

library(Giotto)




# Ensure the Python environment for Giotto has been installed.
genv_exists <- Giotto::checkGiottoEnvironment()
if(!genv_exists){
  # The following command need only be run once to install the Giotto environment.
  Giotto::installGiottoEnvironment()
}

Sys.setenv()

installGiottoEnvironment(
  packages_to_install = c("pandas==1.5.1", "networkx==2.8.8", "python-igraph==0.10.2",
                          "leidenalg==0.9.0", "python-louvain==0.16", "python.app==1.4", "scikit-learn==1.1.3"),
  python_version = "3.10.2",
  mini_install_path = NULL,
  confirm = TRUE,
  envname = "giotto_env",
  conda = "D:/miniconda/envs/giotto_env/",
  force_miniconda = FALSE,
  force_environment = FALSE,
  verbose = NULL
)

checkGiottoEnvironment()



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



checkGiottoEnvironment()
?checkGiottoEnvironment()

# 3. Create Giotto Instructions
python_path <- NULL
instructions <- createGiottoInstructions(save_dir =  '10x_brain/out/',
                                         save_plot = TRUE,
                                         show_plot = FALSE,
                                         return_plot = FALSE,
                                         python_path = python_path)


## provide path to visium folder
data_path <- "10x_brain/"

## directly from visium folder
visium_brain <- createGiottoVisiumObject(visium_dir = data_path,
                                         expr_data = "filter",
                                         png_name = "tissue_lowres_image.png",
                                         gene_column_index = 2,
                                         instructions = instructions)

## show associated images with giotto object
showGiottoImageNames(visium_brain) # "image" is the default name

## check metadata
pDataDT(visium_brain)

## show plot
spatPlot2D(gobject = visium_brain, 
           cell_color = "in_tissue", 
           point_size = 2,
           cell_color_code = c("0" = "lightgrey", "1" = "blue"), 
           show_image = TRUE, 
           image_name = "image")




## subset on spots that were covered by tissue
metadata <- pDataDT(visium_brain)
in_tissue_barcodes <- metadata[in_tissue == 1]$cell_ID

visium_brain <- subsetGiotto(visium_brain, 
                             cell_ids = in_tissue_barcodes)

## filter
visium_brain <- filterGiotto(gobject = visium_brain,
                             expression_threshold = 1,
                             feat_det_in_min_cells = 50,
                             min_det_feats_per_cell = 1000,
                             expression_values = "raw",
                             verbose = TRUE)

## normalize
visium_brain <- normalizeGiotto(gobject = visium_brain, 
                                scalefactor = 6000, 
                                verbose = TRUE)

## add gene & cell statistics
visium_brain <- addStatistics(gobject = visium_brain)

## visualize
spatPlot2D(gobject = visium_brain, 
           show_image = TRUE, 
           point_alpha = 0.7,
           cell_color = "nr_feats", 
           color_as_factor = FALSE)


## 5 highly variable features / genes (HVF)
visium_brain <- calculateHVF(gobject = visium_brain, 
                             save_plot = TRUE)


## run PCA on expression values (default)
gene_metadata <- fDataDT(visium_brain)
featgenes <- gene_metadata[hvf == "yes" & perc_cells > 3 & mean_expr_det > 0.4]$feat_ID

## run PCA on expression values (default)
visium_brain <- runPCA(gobject = visium_brain,
                       feats_to_use = featgenes)

screePlot(visium_brain, 
          ncp = 30)

dimPlot2D(gobject = visium_brain,
          dim_reduction_to_use = "pca")

## run UMAP and tSNE on PCA space (default)
visium_brain <- runUMAP(visium_brain, 
                        dimensions_to_use = 1:10)

plotUMAP(gobject = visium_brain)


visium_brain <- runtSNE(visium_brain, 
                        dimensions_to_use = 1:10)

plotTSNE(gobject = visium_brain)

### 6 Clustering
## sNN network (default)
visium_brain <- createNearestNetwork(gobject = visium_brain, 
                                     dimensions_to_use = 1:10, 
                                     k = 15)

## Leiden clustering
visium_brain <- doLeidenCluster(gobject = visium_brain, 
                                resolution = 0.4, 
                                n_iterations = 1000)
edit(doLeidenCluster)




plotUMAP(gobject = visium_brain,
         cell_color = "leiden_clus", 
         show_NN_network = TRUE, 
         point_size = 2.5)

# spatial and dimension plots
spatDimPlot(gobject = visium_brain, 
            cell_color = "leiden_clus",
            dim_point_size = 2, 
            spat_point_size = 2.5)


spatDimPlot(gobject = visium_brain, 
            cell_color = "nr_feats", 
            color_as_factor = FALSE,
            dim_point_size = 2, 
            spat_point_size = 2.5)

# dimension plots grouped by cluster
spatPlot2D(visium_brain, 
           cell_color = "leiden_clus",
           coord_fix_ratio = 1)



spatPlot2D(visium_brain, 
           cell_color = "leiden_clus",
           group_by = "leiden_clus", 
           coord_fix_ratio = 1,
           cow_n_col = 6, 
           show_legend = FALSE)

spatPlot2D(visium_brain, 
           cell_color = "leiden_clus",
           select_cell_groups = "8", 
           coord_fix_ratio = 1, 
           show_other_cells = TRUE,
           cell_color_code = c("8" = "red"), 
           other_cell_color = "grey", 
           other_point_size = 1.5)

## 7 Subset data
# create and show subset
DG_subset <- subsetGiottoLocs(visium_brain,
                              x_max = 6500, x_min = 3000,
                              y_max = -2500, y_min = -5500,
                              return_gobject = TRUE)

spatDimPlot(gobject = DG_subset,
            cell_color = "leiden_clus", 
            spat_point_size = 5)


#8 Marker gene detection for clusters

## Gini markers
markers_gini <- findMarkers_one_vs_all(gobject = visium_brain,
                                       method = "gini",
                                       expression_values = "normalized",
                                       cluster_column = "leiden_clus",
                                       min_feats = 20,
                                       min_expr_gini_score = 0.5,
                                       min_det_gini_score = 0.5)

topgenes_gini <- markers_gini[, head(.SD, 2), by = "cluster"]$feats

# violinplot
violinPlot(visium_brain, 
           feats = unique(topgenes_gini), 
           cluster_column = "leiden_clus",
           strip_text = 8, 
           strip_position = "right")



# cluster heatmap
plotMetaDataHeatmap(visium_brain, 
                    selected_feats = unique(topgenes_gini),
                    metadata_cols = "leiden_clus",
                    x_text_size = 10, 
                    y_text_size = 10)


# umap plots
dimFeatPlot2D(visium_brain, 
              expression_values = "scaled",
              feats = markers_gini[, head(.SD, 1), by = "cluster"]$feats,
              cow_n_col = 4, 
              point_size = 0.75)



# Scran Markers
markers_scran <- findMarkers_one_vs_all(gobject = visium_brain,
                                        method = "scran",
                                        expression_values = "normalized",
                                        cluster_column = "leiden_clus")

topgenes_scran <- markers_scran[, head(.SD, 2), by = "cluster"]$feats

# violinplot
violinPlot(visium_brain, 
           feats = unique(topgenes_scran), 
           cluster_column = "leiden_clus",
           strip_text = 10, 
           strip_position = "right")
# cluster heatmap
plotMetaDataHeatmap(visium_brain, 
                    selected_feats = topgenes_scran,
                    metadata_cols = "leiden_clus")


# umap plots
dimFeatPlot2D(visium_brain, 
              expression_values = "scaled",
              feats = markers_scran[, head(.SD, 1), by = "cluster"]$feats,
              cow_n_col = 3, 
              point_size = 1)


# 9 Cell type enrichment
# download data to results directory ####
# if wget is installed, set method = "wget"
# if you run into authentication issues with wgeTRUE, then add " extra = "--no-check-certificate" "
#GiottoData::getSpatialDataset(dataset = "scRNA_mouse_brain", 
#                              directory = "10x_brain/out/")
getwd()
sc_expression <- file.path( "brain_sc_expression_matrix.txt.gz")
sc_metadata <- file.path( "brain_sc_metadata.csv")

giotto_SC <- createGiottoObject(expression = sc_expression,
                                instructions = instructions)

giotto_SC <- addCellMetadata(giotto_SC, 
                             new_metadata = data.table::fread(sc_metadata))

giotto_SC <- normalizeGiotto(giotto_SC)

## 9.1 PAGE enrichment

# Create PAGE matrix
# PAGE matrix should be a binary matrix with each row represent a gene marker and each column represent a cell type
# There are several ways to create PAGE matrix
# 1.1 create binary matrix of cell signature genes
# small example #
gran_markers <- c("Nr3c2", "Gabra5", "Tubgcp2", "Ahcyl2",
                  "Islr2", "Rasl10a", "Tmem114", "Bhlhe22", 
                  "Ntf3", "C1ql2")

oligo_markers <- c("Efhd1", "H2-Ab1", "Enpp6", "Ninj2",
                   "Bmp4", "Tnr", "Hapln2", "Neu4",
                   "Wfdc18", "Ccp110")        

di_mesench_markers <- c("Cartpt", "Scn1a", "Lypd6b",  "Drd5",
                        "Gpr88", "Plcxd2", "Cpne7", "Pou4f1",
                        "Ctxn2", "Wnt4")

PAGE_matrix_1 <- makeSignMatrixPAGE(sign_names = c("Granule_neurons",
                                                   "Oligo_dendrocytes",
                                                   "di_mesenchephalon"),
                                    sign_list = list(gran_markers,
                                                     oligo_markers,
                                                     di_mesench_markers))

# ----

# 1.2 [shortcut] fully pre-prepared matrix for all cell types
sign_matrix_path <- system.file("extdata", "sig_matrix.txt", package = "GiottoData")

brain_sc_markers <- data.table::fread(sign_matrix_path) 

PAGE_matrix <- as.matrix(brain_sc_markers[,-1])
rownames(PAGE_matrix) <- brain_sc_markers$Event

# ---

# 1.3 make PAGE matrix from single cell dataset
markers_scran <- findMarkers_one_vs_all(gobject = giotto_SC, 
                                        method = "scran",
                                        expression_values = "normalized",
                                        cluster_column = "Class", 
                                        min_feats = 3)

topgenes_scran <- markers_scran[, head(.SD, 10), by = "cluster"]

celltypes <- levels(factor(markers_scran$cluster)) 

sign_list <- list()

for (i in 1:length(celltypes)){
  sign_list[[i]] <- topgenes_scran[which(topgenes_scran$cluster == celltypes[i]),]$feats
}

PAGE_matrix <- makeSignMatrixPAGE(sign_names = celltypes,
                                  sign_list = sign_list)

# 1.4 enrichment test with PAGE

# runSpatialEnrich() can also be used as a wrapper for all currently provided enrichment options
visium_brain <- runPAGEEnrich(gobject = visium_brain, 
                              sign_matrix = PAGE_matrix)

# 1.5 heatmap of enrichment versus annotation (e.g. clustering result)
cell_types_PAGE <- colnames(PAGE_matrix)

plotMetaDataCellsHeatmap(gobject = visium_brain,
                         metadata_cols = "leiden_clus",
                         value_cols = cell_types_PAGE,
                         spat_enr_names = "PAGE",
                         x_text_size = 8,
                         y_text_size = 8)



# 1.6 visualizations
spatCellPlot2D(gobject = visium_brain,
               spat_enr_names = "PAGE",
               cell_annotation_values = cell_types_PAGE[1:4],
               cow_n_col = 2,
               coord_fix_ratio = 1, 
               point_size = 1.25, 
               show_legend = TRUE)

spatDimCellPlot2D(gobject = visium_brain,
                  spat_enr_names = "PAGE",
                  cell_annotation_values = cell_types_PAGE[1:4],
                  cow_n_col = 1, 
                  spat_point_size = 1,
                  plot_alignment = "horizontal",
                  save_param = list(base_width = 7, base_height = 10))

## 9.2 HyperGeometric test
visium_brain <- runHyperGeometricEnrich(gobject = visium_brain,
                                        expression_values = "normalized",
                                        sign_matrix = PAGE_matrix)

cell_types_HyperGeometric <- colnames(PAGE_matrix)

spatCellPlot(gobject = visium_brain,
             spat_enr_names = "hypergeometric",
             cell_annotation_values = cell_types_HyperGeometric[1:4],
             cow_n_col = 2,
             coord_fix_ratio = NULL, 
             point_size = 1.75)


## 9.3 Rank Enrichment
# Create rank matrix, not that rank matrix is different from PAGE
# A count matrix and a vector for all cell labels will be needed
rank_matrix <- makeSignMatrixRank(sc_matrix = getExpression(giotto_SC,
                                                            values = "normalized",
                                                            output = "matrix"),
                                  sc_cluster_ids = pDataDT(giotto_SC)$Class)

colnames(rank_matrix) <- levels(factor(pDataDT(giotto_SC)$Class))

visium_brain <- runRankEnrich(gobject = visium_brain, 
                              sign_matrix = rank_matrix,
                              expression_values = "normalized")

# Plot Rank enrichment result
spatCellPlot2D(gobject = visium_brain,
               spat_enr_names = "rank",
               cell_annotation_values = colnames(rank_matrix)[1:4],
               cow_n_col = 2,
               coord_fix_ratio = 1, 
               point_size = 1)


# 9.4 DWLS spatial deconvolution
# Create DWLS matrix, not that DWLS matrix is different from PAGE and rank
# A count matrix a vector for a list of gene signatures and a vector for all cell labels will be needed
DWLS_matrix <- makeSignMatrixDWLSfromMatrix(
  matrix = getExpression(giotto_SC,
                         values = "normalized",
                         output = "matrix"),
  cell_type = pDataDT(giotto_SC)$Class,
  sign_gene = topgenes_scran$feats)

visium_brain <- runDWLSDeconv(gobject = visium_brain, 
                              sign_matrix = DWLS_matrix)

# Plot DWLS deconvolution result
spatCellPlot2D(gobject = visium_brain,
               spat_enr_names = "DWLS",
               cell_annotation_values = levels(factor(pDataDT(giotto_SC)$Class))[1:4],
               cow_n_col = 2,
               coord_fix_ratio = 1, 
               point_size = 1)


# Plot DWLS deconvolution result with Pie plots
spatDeconvPlot(visium_brain, 
               show_image = TRUE,
               radius = 50)



## 10 Spatial Grid

visium_brain <- createSpatialGrid(gobject = visium_brain,
                                  sdimx_stepsize = 400,
                                  sdimy_stepsize = 400,
                                  minimum_padding = 0)

showGiottoSpatGrids(visium_brain)

spatPlot2D(visium_brain, 
           cell_color = "leiden_clus", 
           show_grid = TRUE,
           grid_color = "red", 
           spatial_grid_name = "spatial_grid")



##11 Spatial network
visium_brain <- createSpatialNetwork(gobject = visium_brain,
                                     method = "kNN", 
                                     k = 5,
                                     maximum_distance_knn = 400,
                                     name = "spatial_network")

showGiottoSpatNetworks(visium_brain)

spatPlot2D(gobject = visium_brain,  
           show_network= TRUE,
           network_color = "blue", 
           spatial_network_name = "spatial_network")

## 12 Spatial Genes
## rank binarization
ranktest <- binSpect(visium_brain, 
                     bin_method = "rank",
                     calc_hub = TRUE, 
                     hub_min_int = 5,
                     spatial_network_name = "spatial_network")

spatFeatPlot2D(visium_brain, 
               expression_values = "scaled",
               feats = ranktest$feats[1:6], 
               cow_n_col = 2, 
               point_size = 1.5)

## 13 Spatial Co-Expression modules
# cluster the top 500 spatial genes into 20 clusters
my_spatial_genes <- ranktest[1:1500,]$feats

# here we use existing detectSpatialCorGenes function to calculate pairwise distances between genes (but set network_smoothing=0 to use default clustering)
spat_cor_netw_DT <- detectSpatialCorFeats(visium_brain,
                                          method = "network",
                                          spatial_network_name = "spatial_network",
                                          subset_feats = my_spatial_genes)

# 2. identify most similar spatially correlated genes for one gene
top10_genes <- showSpatialCorFeats(spat_cor_netw_DT, 
                                   feats = "Mbp", 
                                   show_top_feats = 10)

spatFeatPlot2D(visium_brain, 
               expression_values = "scaled",
               feats = top10_genes$variable[1:4], 
               point_size = 3)

# cluster spatial genes
spat_cor_netw_DT <- clusterSpatialCorFeats(spat_cor_netw_DT, 
                                           name = "spat_netw_clus", 
                                           k = 20)

# visualize clusters
heatmSpatialCorFeats(visium_brain,
                     spatCorObject = spat_cor_netw_DT,
                     use_clus_name = "spat_netw_clus",
                     heatmap_legend_param = list(title = NULL))


# 4. rank spatial correlated clusters and show genes for selected clusters
netw_ranks <- rankSpatialCorGroups(visium_brain,
                                   spatCorObject = spat_cor_netw_DT, 
                                   use_clus_name = "spat_netw_clus")

top_netw_spat_cluster <- showSpatialCorFeats(spat_cor_netw_DT, 
                                             use_clus_name = "spat_netw_clus",
                                             selected_clusters = 6, 
                                             show_top_feats = 1)


# 5. create metagene enrichment score for clusters
cluster_genes_DT <- showSpatialCorFeats(spat_cor_netw_DT, 
                                        use_clus_name = "spat_netw_clus", 
                                        show_top_feats = 1)

cluster_genes <- cluster_genes_DT$clus 
names(cluster_genes) <- cluster_genes_DT$feat_ID

visium_brain <- createMetafeats(visium_brain, 
                                feat_clusters = cluster_genes, 
                                name = "cluster_metagene")

spatCellPlot(visium_brain,
             spat_enr_names = "cluster_metagene",
             cell_annotation_values = netw_ranks$clusters,
             point_size = 1, 
             cow_n_col = 5)


## 14 Spatially informed clusters
# top 30 genes per spatial co-expression cluster
table(spat_cor_netw_DT$cor_clusters$spat_netw_clus)

coexpr_dt <- data.table::data.table(
  genes = names(spat_cor_netw_DT$cor_clusters$spat_netw_clus),
  cluster = spat_cor_netw_DT$cor_clusters$spat_netw_clus)

data.table::setorder(coexpr_dt, cluster)
top30_coexpr_dt <- coexpr_dt[, head(.SD, 30), by = cluster]

my_spatial_genes <- top30_coexpr_dt$genes

visium_brain <- runPCA(gobject = visium_brain,
                       feats_to_use = my_spatial_genes,
                       name = "custom_pca")

visium_brain <- runUMAP(visium_brain, 
                        dim_reduction_name = "custom_pca", 
                        dimensions_to_use = 1:20,
                        name = "custom_umap")

visium_brain <- createNearestNetwork(gobject = visium_brain,
                                     dim_reduction_name = "custom_pca",
                                     dimensions_to_use = 1:20, 
                                     k = 5,
                                     name = "custom_NN")

visium_brain <- doLeidenCluster(gobject = visium_brain, 
                                network_name = "custom_NN",
                                resolution = 0.15, 
                                n_iterations = 1000,
                                name = "custom_leiden")

cell_metadata <- pDataDT(visium_brain)
cell_clusters <- unique(cell_metadata$custom_leiden)

giotto_colors <- getDistinctColors(length(cell_clusters))
names(giotto_colors) <- cell_clusters

spatPlot2D(visium_brain, 
           cell_color = "custom_leiden", 
           cell_color_code = giotto_colors, 
           coord_fix_ratio = 1)

plotUMAP(gobject = visium_brain, 
         cell_color = "custom_leiden", 
         cell_color_code = giotto_colors, 
         point_size = 1.5)

##15 Spatial domains with HMRF
# do HMRF with different betas on top 30 genes per spatial co-expression module
hmrf_folder <- file.path(data_path, "HMRF")

if(!file.exists(hmrf_folder)) dir.create(hmrf_folder, recursive = TRUE)
library(smfishHmrf)
HMRF_spatial_genes <- doHMRF(gobject = visium_brain,
                             expression_values = "scaled",
                             spatial_genes = my_spatial_genes, 
                             k = 20,
                             spatial_network_name="spatial_network",
                             betas = c(0, 10, 5),
                             output_folder = file.path(hmrf_folder, "Spatial_genes/SG_topgenes_k20_scaled"))

visium_brain <- addHMRF(gobject = visium_brain, 
                        HMRFoutput = HMRF_spatial_genes,
                        k = 20, 
                        betas_to_add = c(0, 10, 20, 30, 40),
                        hmrf_name = "HMRF")

spatPlot2D(gobject = visium_brain, 
           cell_color = "HMRF_k20_b.40")

```



