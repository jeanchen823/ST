# SPATA2

```bash
# packages from biocmanager
if(!requireNamespace("BiocManager", quietly = TRUE)){
  
  install.packages("BiocManager")
      
}

BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'Matrix.utils', 'EBImage'))

# packages from github
if(!requireNamespace("devtools", quietly = TRUE)){
  
  install.packages("devtools")
  
}


library(BiocGenerics)
library(DelayedArray)
library(DelayedMatrixStats)
library(limma)
library(S4Vectors)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(batchelor)
library(Matrix.utils)
library(EBImage)


devtools::install_github(repo = "kueckelj/confuns")
devtools::install_github(repo = "theMILOlab/SPATAData")
devtools::install_github(repo = "theMILOlab/SPATA2")

library(SPATA2)
```



```bash
##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())

library(SPATA2)
library(SPATAData)
library(tidyverse)

setwd("spata/")

# load SPATA2 inbuilt example data
object_t269 <- loadExampleObject(sample_name = "UKF269T", process = TRUE, meta = TRUE)

save(object_t269,file = "object_t269.rdata")

# histology only
plotSurface(object_t269, pt_alpha = 0)

# colored by histological classification
plotSurface(object_t269, color_by = "histology")

object_t269 <- createSpatialTrajectories(object = object_t269)


getTrajectoryWidth(object_t269, id = "huage")

getTrajectoryLength(object_t269, id = "huage")

# created with code 
plotSpatialTrajectories(
  object = object_t269, 
  ids = "huage", 
  color_by = "histology"
)




# created with code 
plotSpatialTrajectories(
  object = object_t269, 
  ids = "horizontal_mid", 
  color_by = "histology"
)


# created with code 
plotSpatialTrajectories(
  object = object_t269, 
  ids = "huage", 
  color_by = "histology"
)

# define start and end positions of the trajectory directly
# by default, the width equals the trajectory length
object_t269 <-
addSpatialTrajectory(
object = object_t269,
id = "horizontal_mid",
start = c("1.5mm", "4mm"),
end = c("6.5mm", "4mm"),
overwrite = TRUE
)


# created with code 
plotSpatialTrajectories(
  object = object_t269, 
  ids = "huage", 
  color_by = "histology"
)
colnames(object_t269@meta_obs)
plotSpatialTrajectories(
  object = object_t269, 
  ids = "huage", 
  color_by = "seurat_clusters"
)

object_t269 <- runSparkx(object = object_t269)

spark_df <- getSparkxGeneDf(object = object_t269, threshold_pval = 0.05)

# show results
spark_df

# `getSparkxGenes()` would work, too
input_genes <- spark_df[["genes"]]


sts_out <- 
  spatialTrajectoryScreening(
    object = object_t269, 
    id = "horizontal_mid", # ID of the spatial trajectory
    variables = input_genes # the variables/genes to include in the screening 
  )

class(sts_out)


sign_df <- 
  sts_out@results$significance %>% 
  filter(fdr < 0.05)

# show significance data.frame
sign_df


# extract variables names
non_random <- getSgsResultsVec(sts_out) %>% head(4)

trajectory_add_on <- 
  ggpLayerSpatialTrajectories(object = object_t269, ids = "horizontal_mid")

plotSurfaceComparison(
  object = object_t269,
  color_by = non_random,
  pt_clrsp = "Reds 3",
  outline = T,
  nrow = 1
) + 
  trajectory_add_on

plotStsLineplot(object_t269, variables = non_random, id = "horizontal_mid", line_color = "red", nrow = 1) 




# extract random variable names 
random <- 
  sts_out@results$significance %>% 
  filter(fdr > 0.05) %>% 
  slice_max(tot_var, n = 4) %>% 
  pull(variables) %>% 
  head(4)

plotSurfaceComparison(
  object = object_t269,
  color_by = random,
  pt_clrsp = "BuPu",
  outline = T, 
  nrow = 1
) + 
  trajectory_add_on


plotStsLineplot(object_t269, variables = random, id = "horizontal_mid", line_color = "blue", nrow = 1)



# left plot
plotStsLineplot(object_t269, variables = "LEPROT", line_color = "blue", id = "horizontal_mid") + 
  ggplot2::geom_point()

# right plot
plotStsLineplot(object_t269, variables = "SHISA5", line_color = "red", id = "horizontal_mid") + 
  ggplot2::geom_point()


showModels(nrow = 3) + 
  labs(x = "Distance along Trajectory [%]")


best_fits <- 
  sts_out@results$model_fits %>% 
  filter(variables %in% sign_df[["variables"]]) %>% 
  group_by(variables) %>% 
  slice_min(mae, n = 1)
best_fits

best_fits_by_model <- 
  group_by(best_fits, models) %>% 
  slice_min(mae, n = 1) %>% 
  filter(rmse < 0.2) # threshold suggestions for root mean squared error

best_fits_by_model

plotSurfaceComparison(
  object = object_t269, 
  color_by = best_fits_by_model[["variables"]], 
  outline = TRUE, 
  display_image = FALSE, 
  pt_clrsp = "Reds 3", 
  nrow = 2
) + 
  trajectory_add_on

plotStsLineplot(
  object = object_t269, 
  variables = best_fits_by_model[["variables"]], 
  id = "horizontal_mid", 
  line_color = "red", 
  nrow = 2
)
```



```bash

##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())

library(SPATA2)
library(SPATAData)
library(tidyverse)

setwd("spata/")


object  <- 
  initiateSpataObjectVisium(
    sample_name = "UKF269T", 
    directory_visium = "GBM1_spaceranger_out/" # adjust to your liking 
  )


object <- identifyPixelContent(object)

plotImageMask(object)

plotPixelContent(object)

plotImage(object)

#plotImage(object, outline = TRUE, line_size = 1)

object <- identifyTissueOutline(object, method = "obs", eps = "125um", minPts = 3)

plotSurface(object, color_by = "tissue_section", pt_clrp = "tab20")



## 4.过滤低质量的spot
# uses the results of identifyTissueOutline() to create a logical variable called sp_outlier
object <- identifySpatialOutliers(object, method = "obs")

plot_with_outliers <- plotSurface(object, color_by = "sp_outlier", clrp_adjust = c("TRUE" = "blue"))

# remove where sp_outlier == TRUE
object <- removeSpatialOutliers(object)

plot_without_outliers <- plotSurface(object, color_by = "sp_outlier")

# left plot
plot_with_outliers

# right plot
plot_without_outliers


# before
nGenes(object)
## [1] 36601

# removes stress genes
object <- removeGenesStress(object)

# removes genes that were not detected in any of the observations
object <- removeGenesZeroCounts(object)

# afterwards
nGenes(object)
## [1] 25055


object <- computeMetaFeatures(object)

# plot left
plotSurface(object, color_by = "n_counts_gene")

# plot right
plotSurface(object, color_by = "n_distinct_gene")


# obtain matrix names prior to normalization
getMatrixNames(object)
## [1] "counts"

plot_before <- 
  plotSurface(object, color_by = "MAG") + labs(color = "MAG\n(Counts)")

# create log normalized matrix
object <- normalizeCounts(object, method = "LogNormalize")
## Normalizing layer: counts
## 14:20:34 Active matrix in assay 'gene': 'LogNormalize'

# obtain matrix names after normalization
getMatrixNames(object)
## [1] "counts"       "LogNormalize"

# check active matrix 
activeMatrix(object)
## [1] "LogNormalize"

plot_afterwards <- 
  plotSurface(object, color_by = "MAG") + labs(color = "MAG\n(logNorm)")

# left plot
plot_before

# right plot
plot_afterwards


# results are stored inside the SPATA2 object
object <- runSPARKX(object, verbose = FALSE)


# get genes with a p-value < 0.01
sparkx_genes <- getSparkxGenes(object, threshold_pval = 0.01)
str(sparkx_genes)

# visualize in space
plotSurfaceComparison(object, color_by = head(sparkx_genes, 6), nrow = 2)


# total number of genes in this (subsetted) object
nGenes(object)
## [1] 25055

# identify most variable ones (using Seurat in the background)
object <- identifyVariableMolecules(object, n_mol = 2500, method = "vst")

# variable mols
vm <- getVariableMolecules(object, method = "vst")

head(vm)

length(vm)

# run the algorithm
object <- runPCA(object, variables = vm, n_pcs = 20)

plotPcaElbow(object)

###########################################################################
##  10.聚类
# current grouping options 
getGroupingOptions(object )

# run PCA based on which clustering is conducted
object <- runPCA(object, n_pcs = 20)

object <- 
  runKmeansClustering(
    object = object, 
    ks = c(7, 8), 
    methods_kmeans = "Lloyd"
  )

# results are immediately stored in the objects feature data
getGroupingOptions(object)


plotSurface(
  object = object, 
  color_by = "Lloyd_k8", 
  pt_clrp = "uc"
)

# right plot
plotSurface(
  object = object, 
  color_by = "Lloyd_k7",
  pt_clrp = "jco"
)


###########################################################################
##  11  tsne和umap可视化降维

# run dimensional reduction
object <- runTSNE(object, n_pcs = 10)
object <- runUMAP(object, n_pcs = 10)

# left plot
plotTSNE(object, color_by = "Lloyd_k7")

# right plot
plotUMAP(object, color_by = "Lloyd_k7")




object_t269 <- createSpatialTrajectories(object = object)

object_t269 <- 
  addSpatialTrajectory(
    object =  object,
    id = "huage",
    start = c("1.5mm", "4mm"),
    end = c("6.5mm", "4mm"),
    overwrite = TRUE
  )


# created with code 
plotSpatialTrajectories(
  object = object_t269, 
  ids = "huage", 
  color_by = "Lloyd_k7"
)




object_t269 <- runSPARKX(object = object_t269)
#保留sparkx值为0.01或更低的基因
spark_df <- getSparkxGeneDf(object = object_t269, threshold_pval = 0.01)
# 提取基因
input_genes <- spark_df[["genes"]]

#筛选沿着空间轨迹的过程跟随特定表达式变化的数值变量样本
sts_out <- 
  spatialTrajectoryScreening(
    object = object_t269, 
    id = "huage", # ID of the spatial trajectory
    variables = input_genes # the variables/genes to include in the screening 
  )



sign_df <- 
  sts_out@results$significance %>% 
  filter(fdr < 0.05)

# extract variables names
non_random <- getSgsResultsVec(sts_out) %>% head(4)

trajectory_add_on <- 
  ggpLayerSpatialTrajectories(object = object_t269, ids = "huage")

plotSurfaceComparison(
  object = object_t269,
  color_by = non_random,
  pt_clrsp = "Reds 3",
  outline = T,
  nrow = 1
) + 
  trajectory_add_on

plotStsLineplot(object_t269, variables = non_random, id = "huage", line_color = "red", nrow = 1) 




# extract random variable names 
random <- 
  sts_out@results$significance %>% 
  filter(fdr > 0.05) %>% 
  slice_max(tot_var, n = 4) %>% 
  pull(variables) %>% 
  head(4)

plotSurfaceComparison(
  object = object_t269,
  color_by = random,
  pt_clrsp = "BuPu",
  outline = T, 
  nrow = 1
) + 
  trajectory_add_on


plotStsLineplot(object_t269, variables = random, id = "huage", line_color = "blue", nrow = 1)



# left plot
plotStsLineplot(object_t269, variables = "LEPROT", line_color = "blue", id = "huage") + 
  ggplot2::geom_point()

# right plot
plotStsLineplot(object_t269, variables = "SHISA5", line_color = "red", id ="huage") + 
  ggplot2::geom_point()


showModels(nrow = 3) + 
  labs(x = "Distance along Trajectory [%]")


best_fits <- 
  sts_out@results$model_fits %>% 
  filter(variables %in% sign_df[["variables"]]) %>% 
  group_by(variables) %>% 
  slice_min(mae, n = 1)
best_fits

best_fits_by_model <- 
  group_by(best_fits, models) %>% 
  slice_min(mae, n = 1) %>% 
  filter(rmse < 0.2) # threshold suggestions for root mean squared error

best_fits_by_model

plotSurfaceComparison(
  object = object_t269, 
  color_by = best_fits_by_model[["variables"]], 
  outline = TRUE, 
  display_image = FALSE, 
  pt_clrsp = "Reds 3", 
  nrow = 2
) + 
  trajectory_add_on

plotStsLineplot(
  object = object_t269, 
  variables = best_fits_by_model[["variables"]], 
  id = "huage", 
  line_color = "red", 
  nrow = 2
)

```




