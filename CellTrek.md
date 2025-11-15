# CellTrek

```bash
1.安装R包 
library(devtools)
install_github("navinlabcode/CellTrek")
2.加载R包、设置路径
##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())

#library(devtools)
#install_github("navinlabcode/CellTrek")
library("CellTrek")
library("dplyr")
library("Seurat")
library("viridis")
library("ConsensusClusterPlus")
library(SeuratData)
library(ggplot2)
library(patchwork)
library(readr)
library(ConsensusClusterPlus)
library(SeuratObject)

setwd("CellTrek/")

brain_st_cortex <- read_rds("brain_st_cortex.rds")
brain_sc <- read_rds("brain_sc.rds")

3.使用语法上有效的名称重命名细胞和spot
## Rename the cells/spots with syntactically valid names
brain_st_cortex <- RenameCells(brain_st_cortex, new.names=make.names(Cells(brain_st_cortex)))
brain_sc <- RenameCells(brain_sc, new.names=make.names(Cells(brain_sc)))

4.可视化
## Visualize the ST data
SpatialDimPlot(brain_st_cortex)

## Visualize the scRNA-seq data
DimPlot(brain_sc, label = T, label.size = 4.5)
```

```bash
1.我们首先使用训练方法共同嵌入ST和scRNA-seq数据集

brain_traint <- CellTrek::traint(st_data=brain_st_cortex, 
                                 sc_data=brain_sc, 
                                 sc_assay='RNA', 
                                 cell_names='cell_type')

2.查看共嵌入的结果
## 我们可以检查共嵌入的结果，以查看这两种数据模态之间是否存在重叠。
DimPlot(brain_traint, group.by = "type")
    
3.共在共嵌入之后，我们可以将单个细胞映射到它们的空间位置。在这里，我们使用非线性插值（intp = T，intp_lin=F）方法来增强ST的spots

brain_celltrek <- celltrek(st_sc_int=brain_traint, int_assay='traint', sc_data=brain_sc, sc_assay = 'RNA', 
                                     reduction='pca', intp=T, intp_pnt=5000, intp_lin=F, nPCs=30, ntree=1000, 
                                     dist_thresh=0.55, top_spot=5, spot_n=5, repel_r=20, repel_iter=20, keep_model=T)$celltrek
                                    
4.细胞映射完成后，可以使用celltrek_vis交互式可视化CellTrek的结果
brain_celltrek$cell_type <- factor(brain_celltrek$cell_type, 
                                   levels=sort(unique(brain_celltrek$cell_type)))

celltrek_vis(brain_celltrek@meta.data %>% 
                         dplyr::select(coord_x, coord_y, cell_type:id_new),
                       brain_celltrek@images$anterior1@image, 
                       brain_celltrek@images$anterior1@scale.factors$lowres)

```

```bash
1.提取子集
基于CellTrek的结果，我们可以使用SColoc模块总结不同细胞类型之间的共定位模式。在这里，我们以谷氨酸能神经元细胞类型为例(建议去除一些细胞很少的细胞类型，例如n<20)。我们首先从我们的图表结果中划分出谷氨酸能神经元细胞类型的子集。

glut_cell <- c('L2/3 IT', 'L4', 'L5 IT', 'L5 PT', 'NP', 'L6 IT', 'L6 CT',  'L6b')
names(glut_cell) <- make.names(glut_cell)
brain_celltrek_glut <- subset(brain_celltrek, subset=cell_type %in% glut_cell)
brain_celltrek_glut$cell_type <- factor(brain_celltrek_glut$cell_type, levels=glut_cell)
2.然后使用scoloc进行共定位分析
brain_sgraph_KL <-   scoloc(brain_celltrek_glut, col_cell='cell_type', use_method='KL', eps=1e-50)

3.从图中提取最小生成树（MST）的结果
brain_sgraph_KL_mst_cons <- brain_sgraph_KL$mst_cons
rownames(brain_sgraph_KL_mst_cons) <- colnames(brain_sgraph_KL_mst_cons) <- glut_cell[colnames(brain_sgraph_KL_mst_cons)]
4.提取meta.data数据
brain_cell_class <- brain_celltrek@meta.data %>% dplyr::select(id=cell_type) %>% unique
brain_celltrek_count <- data.frame(freq = table(brain_celltrek$cell_type))
brain_cell_class_new <- merge(brain_cell_class, brain_celltrek_count, by.x ="id", by.y = "freq.Var1")

5.可视化共定位结果。请随意调整边缘值截止。
CellTrek::scoloc_vis(brain_sgraph_KL_mst_cons,
                     meta_data=brain_cell_class_new )
```

```bash
1.提取感兴趣的细胞类型

brain_celltrek_l5 <- subset(brain_celltrek, subset=cell_type=='L5 IT')
brain_celltrek_l5[["RNA"]] <- as(object =brain_celltrek_l5[["RNA"]], Class = "Assay")
brain_celltrek_l5@assays$RNA@scale.data <- matrix(NA, 1, 1)
brain_celltrek_l5$cluster <- gsub('L5 IT VISp ', '', brain_celltrek_l5$cluster)
DimPlot(brain_celltrek_l5, group.by = 'cluster')

2.选择前2000个可变基因(不包括线粒体、核糖体和高零基因)

brain_celltrek_l5 <- FindVariableFeatures(brain_celltrek_l5)
vst_df <- brain_celltrek_l5@assays$RNA@meta.features %>% data.frame %>% mutate(id=rownames(.))
nz_test <- apply(as.matrix(brain_celltrek_l5[['RNA']]@data), 1, function(x) mean(x!=0)*100)
hz_gene <- names(nz_test)[nz_test<20]
mt_gene <- grep('^Mt-', rownames(brain_celltrek_l5), value=T)
rp_gene <- grep('^Rpl|^Rps', rownames(brain_celltrek_l5), value=T)
vst_df <- vst_df %>% dplyr::filter(!(id %in% c(mt_gene, rp_gene, hz_gene))) %>% arrange(., -vst.variance.standardized)
feature_temp <- vst_df$id[1:2000]
3.使用 scoexp 进行空间加权基因共表达分析
brain_celltrek_l5_scoexp_res_cc <- CellTrek::scoexp(celltrek_inp=brain_celltrek_l5, 
                             assay='RNA', 
                               approach='cc', 
                                gene_select = feature_temp, 
                                 sigm=140, 
                                  avg_cor_min=.4, 
                                    zero_cutoff=3, 
                                  min_gen=40, max_gen=400)

4.使用热图可视化共表达式模块

brain_celltrek_l5_k <- rbind(data.frame(gene=c(brain_celltrek_l5_scoexp_res_cc$gs[[1]]), G='K1'), 
                           data.frame(gene=c(brain_celltrek_l5_scoexp_res_cc$gs[[2]]), G='K2')) %>% 
                           magrittr::set_rownames(.$gene) %>% dplyr::select(-1)
pheatmap::pheatmap(brain_celltrek_l5_scoexp_res_cc$wcor[rownames(brain_celltrek_l5_k), rownames(brain_celltrek_l5_k)], 
                   clustering_method='ward.D2', annotation_row=brain_celltrek_l5_k, show_rownames=F, show_colnames=F, 
                   treeheight_row=10, treeheight_col=10, annotation_legend = T, fontsize=8,
                   color=viridis(10), main='L5 IT spatial co-expression')
                  
5.确定了两个不同的模块。基于我们识别的共表达模块，我们可以计算模块的分数
brain_celltrek_l5 <- AddModuleScore(brain_celltrek_l5, 
                                    features=brain_celltrek_l5_scoexp_res_cc$gs, 
                                    name='CC_', nbin=10, ctrl=50, seed=42)
## 可视化1
FeaturePlot(brain_celltrek_l5, 
            grep('CC_', colnames(brain_celltrek_l5@meta.data), 
                 value=T), ncol = 1)


## 可视化2
SpatialFeaturePlot(brain_celltrek_l5, 
                   grep('CC_', colnames(brain_celltrek_l5@meta.data), value=T))


```

```bash
06
 
 细胞类型空间距离的计算

1.提取空间坐标
inp_df <- brain_celltrek_glut@meta.data %>% dplyr::select(cell_names = dplyr::one_of('cell_type'), 
                                                          coord_x, coord_y)
inp_df$coord_x = 270-inp_df$coord_x
head(inp_df)
2.距离计算
output <- kdist(inp_df = inp_df, 
                ref = "L2/3 IT", #目标细胞类型
                ref_type = 'all', 
                que = glut_cell,  #其余的细胞
                k = 10, 
                new_name = "L23ITvs.Others",
                keep_nn = F)
## 查看计算结果
head(output$kdist_df)
3.metda信息与距离结果合并

res = output$kdist_df
res$barcode = row.names(res)
inp_df$barcode = row.names(inp_df)
res = left_join(res, inp_df)
head(res)
4.画图展示亚群空间距离差异

library(ggpubr)
ggboxplot(data = res, 
          x = "cell_names",
          y = "L23ITvs.Others", 
          fill = "cell_names", 
          title = "K-distance to L2/3 IT cells")+ 
  stat_compare_means(method = "kruskal.test") +
  theme(plot.title = element_text(color="black",hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5), #,vjust = 0.5
        legend.position = "none") + labs(y = "Distance")
```

```bash
https://www.nature.com/articles/s41586-023-06252-9

https://github.com/navinlabcode/HumanBreastCellAtlas/blob/main/R/figure2_andOtherSpatial.r

## Spatial cell proximity ##
resolve_hbca_srt_merge_sp_df_dt <- CellTrek:::DT_boot_mst(resolve_hbca_srt_merge_sp_df[, c(3, 4)], 
                                                          coord_df=resolve_hbca_srt_merge_sp_df[, c(1, 2)], col_cell='cell_type', boot_n=20, dist_cutoff=1000)

resolve_hbca_srt_merge_sp_df_meta <- data.frame(id=rownames(resolve_hbca_srt_merge_sp_df_dt$mst_cons), 
                                                type=c('epi', 'strm', 'epi', 'epi', 'endo', 'imm', 'strm', 'imm', 'endo', 'imm'), 
                                                freq=c(table(resolve_hbca_merge$celltype)[rownames(resolve_hbca_srt_merge_sp_df_dt$mst_cons)]))

resolve_hbca_srt_merge_sp_df_dt_boot <- apply(resolve_hbca_srt_merge_sp_df_dt$boot_array, c(1, 2), mean)
diag(resolve_hbca_srt_merge_sp_df_dt_boot) <- NA
resolve_hbca_srt_merge_sp_df_dt_boot <- max(resolve_hbca_srt_merge_sp_df_dt_boot, na.rm = T) - resolve_hbca_srt_merge_sp_df_dt_boot
CellTrek::scoloc_vis(resolve_hbca_srt_merge_sp_df_dt_boot, meta_data = resolve_hbca_srt_merge_sp_df_meta)


```



```bash
https://www.nature.com/articles/s41588-024-01802-x

https://github.com/amin69upenn/Human_Kidney_Multiomics_and_Spatial_Atlas_/blob/main/spRNA-seq_Seurat


#-----------------------------------Mapping Back snRNA-seq, scRNA-seq, and snATAC-seq to Spatial dataset----------------------#
HK.SC.SN.ATAC <- HK.SC.SN.ATAC[, sample(colnames(HK.SN), size =20000, replace=F)] #Downsample the annotated snRNA-seq matrix
Control1.ST = Load10X_Spatial("/Your_10X-Folder/outs", filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "slice1")
Control1.ST <- RenameCells(Control1.ST, new.names=make.names(Cells(Control1.ST)))
HK.SC.SN.ATAC <- RenameCells(HK.SC.SN.ATAC, new.names=make.names(Cells(HK.SC.SN.ATAC)))
Coltrol1.ST.Samples_traint <- CellTrek::traint(st_data=Control1.ST, sc_data=HK.SC.SN.ATAC, sc_assay='RNA', cell_names='seurat_clusters')
Control1.ST.HK.SN_celltrek <- CellTrek::celltrek(st_sc_int=Coltrol1.ST.Samples_traint, int_assay='traint', sc_data=HK.SC.SN.ATAC, sc_assay = 'RNA', reduction='pca', intp=T, intp_pnt=5000, intp_lin=F, nPCs=30, ntree=1000, dist_thresh=0.55, top_spot=5, spot_n=5, repel_r=20, repel_iter=20, keep_model=T)$celltrek
Control1.ST.HK.SN_celltrek$seurat_clusters <- factor(Control1.ST.HK.SN_celltrek$seurat_clusters, levels=sort(Control1.ST.HK.SN_celltrek$seurat_clusters)))
Idents = Control1.ST.HK.SN_celltrek[["Idents2"]]
Idents (Control1.ST.HK.SN_celltrek) = Idents
#Do this for all samples

```



