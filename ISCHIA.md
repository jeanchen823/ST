# ISCHIA

```bash
https://www.embopress.org/doi/full/10.1038/s44320-023-00006-5

##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())

## Install package
devtools::install_github("ati-lz/ISCHIA")

## Load required packages
library(ISCHIA)
library(robustbase)
library(data.table)
library(ggplot2)
library(Seurat)
library(dplyr)

## Set working directory
setwd("")

## Load Seurat object
load("IBD_visium_SeuratObj_small.RData")
IBD.visium.P4

## Prepare deconvolution matrix
deconv.mat <- as.matrix(IBD.visium.P4@meta.data[, 9:28])
colnames(deconv.mat) <- sapply(colnames(deconv.mat), function(x) unlist(strsplit(x, split = "_"))[2])
head(deconv.mat)

Composition.cluster.k(deconv.mat, 20)

## Perform composition clustering
IBD.visium.P4 <- Composition.cluster(IBD.visium.P4, deconv.mat, 8)
table(IBD.visium.P4$CompositionCluster_CC)

## Plot SpatialDimPlot
SpatialDimPlot(IBD.visium.P4, group.by = c("CompositionCluster_CC")) +
  scale_fill_manual(values = c("cyan", "orange", "purple", "green", "yellow", "blue", "red", "black"))

## Plot enriched cell types
Composition_cluster_enrichedCelltypes(IBD.visium.P4, "CC4", deconv.mat)
Composition_cluster_enrichedCelltypes(IBD.visium.P4, "CC7", deconv.mat)

## Perform UMAP on composition clusters
IBD.visium.P4.umap <- Composition_cluster_umap(IBD.visium.P4, deconv.mat)

## Plot UMAP results
IBD.visium.P4.umap$umap.cluster.gg
IBD.visium.P4.umap$umap.deconv.gg

## Calculate spatial co-occurrence for CC4 and CC7
CC4.celltype.cooccur <- spatial.celltype.cooccurence(
  spatial.object = IBD.visium.P4,
  deconv.prob.mat = deconv.mat,
  COI = "CC4",
  prob.th = 0.05,
  Condition = unique(IBD.visium.P4$orig.ident)
)
plot.celltype.cooccurence(CC4.celltype.cooccur)

CC7.celltype.cooccur <- spatial.celltype.cooccurence(
  spatial.object = IBD.visium.P4,
  deconv.prob.mat = deconv.mat,
  COI = "CC7",
  prob.th = 0.05,
  Condition = unique(IBD.visium.P4$orig.ident)
)
plot.celltype.cooccurence(CC7.celltype.cooccur)

## Read ligand-receptor network
lr_network = readRDS("lr_network.rds")
all.LR.network <- cbind(lr_network[, c("from", "to")], LR_Pairs = paste(lr_network$from, lr_network$to, sep = "_"))
all.LR.network.exp <- all.LR.network[which(all.LR.network$from %in% rownames(IBD.visium.P4) & all.LR.network$to %in% rownames(IBD.visium.P4)), ]

## Sample ligand-receptor interactions
all.LR.network.exp <- sample_n(all.LR.network.exp, 500)
all.LR.genes <- unique(c(all.LR.network.exp$from, all.LR.network.exp$to))
all.LR.genes.comm <- intersect(all.LR.genes, rownames(IBD.visium.P4))
LR.pairs <- all.LR.network.exp$LR_Pairs
LR.pairs.AllCombos <- combn(all.LR.genes.comm, 2, paste0, collapse = "_")
LR.pairs.AllCombos

## Calculate enriched ligand-receptor pairs for CC4 and CC7
CC4.Enriched.LRs <- Enriched.LRs(
  IBD.visium.P4, c("CC4"), unique(IBD.visium.P4$orig.ident), 
  all.LR.genes.comm, LR.pairs, 1, 0.2
)
CC7.Enriched.LRs <- Enriched.LRs(
  IBD.visium.P4, c("CC7"), unique(IBD.visium.P4$orig.ident), 
  all.LR.genes.comm, LR.pairs, 1, 0.2
)

## Compare enriched ligand-receptor pairs between CC4 and CC7
CC1vsCC4.Enriched.LRs.Specific <- Diff.cooc.LRs(CC4.Enriched.LRs, CC7.Enriched.LRs, 0.05, 0.1)

## Plot enriched ligand-receptor pairs
ChordPlot.Enriched.LRs(CC4.Enriched.LRs$COI.enrcihed.LRs[1:20, ])
SankeyPlot.Diff.LRs(CC4.Enriched.LRs$COI.enrcihed.LRs, CC7.Enriched.LRs$COI.enrcihed.LRs)

```







