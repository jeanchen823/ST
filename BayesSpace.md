# BayesSpace

```bash

https://github.com/LieberInstitute/spatialDLPFC/blob/main/code/analysis/03_BayesSpace/01_BayesSpace.R

https://github.com/edward130603/BayesSpace/blob/master/vignettes/BayesSpace.Rmd

https://www.bioconductor.org/packages/release/bioc/vignettes/BayesSpace/inst/doc/BayesSpace.html

options(BioC_mirror = "https://mirrors.westlake.edu.cn/bioconductor")
BiocManager::install("BayesSpace")
BiocManager::install("SpatialExperiment")

##系统报错改为英文
Sys.setenv(LANGUAGE = "en")

##禁止转化为因子
options(stringsAsFactors = FALSE)

##清空环境
rm(list = ls())

library(SpatialExperiment)
library(SingleCellExperiment)
library(ggplot2)
library(BayesSpace)
library(cowplot)
library(patchwork)
library(dplyr)
library(ggpubr)
library(pheatmap)
library(Seurat)
library(SeuratObject)
library(tidyverse)

set.seed(7788)

sp1 <- readVisium("GBM1_spaceranger_out/")
sp1

set.seed(102)
sp1 <- scater::logNormCounts(sp1)
sp1 <- spatialPreprocess(sp1, platform = "Visium", 
                         n.PCs = 15, n.HVGs = 2000, log.normalize = FALSE)

sp1 <- qTune(sp1, qs = seq(2, 10), platform = "Visium", d = 7)
qPlot(sp1)

set.seed(149)
sp1 <- spatialCluster(sp1, q = 7, platform = "Visium", d = 10, 
                      init.method = "mclust", model = "t", gamma = 2, 
                      nrep = 1000, burn.in = 100, save.chain = TRUE)

clusterPlot(sp1)
clusterPlot(sp1, palette = c('#3cb44b', '#1f77b4', '#000075', '#ff7f0e', 
                             '#aa40fc', '#d62728', '#c49c94'), color = "black") +  
  theme_bw() + xlab("Column") + ylab("Row") + 
  labs(fill = "BayesSpace\ncluster", title = "Spatial clustering of ST_mel1_rep2")

sp1.enhanced <- spatialEnhance(sp1, q = 7, platform = "Visium", d = 10, 
                               model = "t", gamma = 2, jitter_prior = 0.3, 
                               jitter_scale = 3.5, nrep = 1000, burn.in = 100, 
                               save.chain = TRUE)

mcmcChain(sp1.enhanced, "Ychange")
clusterPlot(sp1.enhanced)

markers <- c("PMEL", "CD2", "CD19", "COL1A1")
sp1.enhanced <- enhanceFeatures(sp1.enhanced, sp1, feature_names = markers, nrounds = 0)

featurePlot(sp1.enhanced, "PMEL")

enhanced.plots <- purrr::map(markers, function(x) featurePlot(sp1.enhanced, x))
patchwork::wrap_plots(enhanced.plots, ncol = 2)

spot.plots <- purrr::map(markers, function(x) featurePlot(sp1, x))
patchwork::wrap_plots(c(enhanced.plots, spot.plots), ncol = 4)

```




