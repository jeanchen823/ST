# ingest&BBKNN

```bash
ingest
04 使用ingest映射到reference批次

import scanpy as sc
import pandas as pd
import seaborn as sns
import os
sc.settings.verbosity = 1             # verbosity errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
sc.settings.set_figure_params(dpi=80, frameon=False, figsize=(3, 3), facecolor='white')
print(os.getcwd())

## 更改完再检查一下有没有更改成功
print(os.getcwd())

 # PBMCs 3k(已经处理)
adata_ref = sc.read_h5ad('pbmc3k.h5ad')
adata_ref
adata = sc.read_h5ad('pbmc3k_1.h5ad')


# 求交集
var_names = adata_ref.var_names.intersection(adata.var_names)
var_names
#对数据提取
adata_ref = adata_ref[:, var_names]
adata = adata[:, var_names]


sc.pp.pca(adata_ref)
sc.pp.neighbors(adata_ref)
sc.tl.umap(adata_ref)
adata_ref


sc.pl.umap(adata_ref, color='louvain')



sc.tl.ingest(adata, adata_ref, obs='louvain')


参数embedding_method：
使用adata_ref的umap或pca参数（即adata_ref.uns['pca']和adata_ref.uns['umap']）去映射得到adata的embedding。

参数labeling_method：
使用knn去映射得到adata的label（训练数据或参考数据为adata_ref的embedding，测试数据为映射后的adata的embedding，根据相似度赋予adata标签）。

参数obs：
映射使用的标签键，比如cell_type，louvain等（从adata_ref映射到adata中）。
在调用ingest前，需要对adatat_ref运行neighbors()。（neighbors()用于计算观测的neighborhood graph，这是ingest计算前的处理）


adata.uns['louvain_colors'] = adata_ref.uns['louvain_colors']
adata.uns['louvain_colors']



sc.pl.umap(adata, color=['louvain', 'bulk_labels'], wspace=0.5)

df1=adata_ref.obs


adata_concat = adata_ref.concatenate(adata, batch_categories=['ref', 'new'])
df2=adata_concat.obs
## 多了一列信息--batch 区分了数据来源



adata_concat = adata_ref.concatenate(adata, batch_categories=['ref', 'new'])
df2=adata_concat.obs
## 多了一列信息--batch 区分了数据来源



adata_concat.obs.louvain=adata_concat.obs.louvain.cat.reorder_categories(adata_ref.obs.louvain.cat.categories )


adata_concat.obs.louvain.cat.categories



adata_concat.uns['louvain_colors'] = adata_ref.uns['louvain_colors']
## 画图展示
sc.pl.umap(adata_concat, color=['batch', 'louvain'])



sc.tl.pca(adata_concat)

sc.external.pp.bbknn(adata_concat, batch_key='batch')

sc.tl.umap(adata_concat)
sc.pl.umap(adata_concat, color=['batch', 'louvain'])

pip install bbknn==1.6.0


```



```bash
bbknn


import scanpy as sc
import pandas as pd
import seaborn as sns
sc.settings.verbosity = 1             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
sc.settings.set_figure_params(dpi=80, frameon=False, figsize=(3, 3), facecolor='white')

adata_all = sc.read('pancreas.h5ad' )

adata_all.shape


counts = adata_all.obs.celltype.value_counts()
counts


minority_classes = counts.index[-5:].tolist()        # get the minority classes
adata_all = adata_all[                               # actually subset
    ~adata_all.obs.celltype.isin(minority_classes)]  # ~表示取反
adata_all.obs.celltype.cat.reorder_categories(       # reorder according to abundance
    counts.index[:-5].tolist(), inplace=True)

counts = adata_all.obs.celltype.value_counts()
counts
## 少了五个细胞类型

sc.pp.pca(adata_all)
sc.pp.neighbors(adata_all)
sc.tl.umap(adata_all)
# palette用于选择颜色
sc.pl.umap(adata_all, color=['batch', 'celltype'], palette=sc.pl.palettes.vega_20_scanpy)


sc.external.pp.bbknn(adata_all, batch_key='batch')
sc.tl.umap(adata_all)
sc.pl.umap(adata_all, color=['batch', 'celltype'])


```



```bash
bbknn
---ingest


adata_ref = adata_all[adata_all.obs.batch == '0']
计算PCA，邻居图 UMAP：
sc.pp.pca(adata_ref)
sc.pp.neighbors(adata_ref)
sc.tl.umap(adata_ref)
sc.pl.umap(adata_ref, color='celltype')

adatas = [adata_all[adata_all.obs.batch == i].copy() for i in ['1', '2', '3']]


sc.settings.verbosity = 2  # a bit more logging
for iadata, adata in enumerate(adatas):
    print(f'... integrating batch {iadata+1}')
    # celltype_orig保存批次原始的细胞类型
    adata.obs['celltype_orig'] = adata.obs.celltype
    sc.tl.ingest(adata, adata_ref, obs='celltype')

adata_concat = adata_ref.concatenate(adatas)

adata_concat.obs.celltype = adata_concat.obs.celltype.astype('category')
adata_concat.obs.celltype.cat.reorder_categories(adata_ref.obs.celltype.cat.categories )  # fix category ordering
adata_concat.uns['celltype_colors'] = adata_ref.uns['celltype_colors']  # fix category coloring

sc.pl.umap(adata_concat, color=['batch', 'celltype'])

09

adata_query = adata_concat[adata_concat.obs.batch.isin(['1', '2', '3'])]
sc.pl.umap(adata_query, color=['batch', 'celltype', 'celltype_orig'], wspace=0.4)


obs_query = adata_query.obs
conserved_categories = obs_query.celltype.cat.categories.intersection(obs_query.celltype_orig.cat.categories)  # intersected categories
obs_query_conserved = obs_query.loc[obs_query.celltype.isin(conserved_categories) & obs_query.celltype_orig.isin(conserved_categories)]  # intersect categories
obs_query_conserved.celltype=obs_query_conserved.celltype.cat.remove_unused_categories()  # remove unused categoriyes

obs_query_conserved.celltype_orig=obs_query_conserved.celltype_orig.cat.remove_unused_categories()  # remove unused categoriyes
obs_query_conserved.celltype_orig=obs_query_conserved.celltype_orig.cat.reorder_categories(obs_query_conserved.celltype.cat.categories)  # fix category ordering
pd.crosstab(obs_query_conserved.celltype, obs_query_conserved.celltype_orig)


10

sc.tl.embedding_density(adata_concat, groupby='batch')
sc.pl.embedding_density(adata_concat, groupby='batch')


for batch in ['1', '2', '3']:
    sc.pl.umap(adata_concat, color='batch', groups=[batch])


```





```bash
11 -harmony整合
https://github.com/slowkow/harmonypy


import scanpy as sc
import pandas as pd
import seaborn as sns
import os
sc.settings.verbosity = 1             # verbosity errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
sc.settings.set_figure_params(dpi=80, frameon=False, figsize=(3, 3), facecolor='white')

print(os.getcwd())
os.chdir('F://spe.lesson//k7')
## 更改完再检查一下有没有更改成功
print(os.getcwd())

data_s10=sc.read_10x_mtx('BC10\\')
data_s21=sc.read_10x_mtx('BC21\\')
data_s2=sc.read_10x_mtx('BC2\\')

data_s10.obs["orig"]='s10'
data_s21.obs["orig"]='s21'
data_s2.obs["orig"]='s2'

adatas = sc.AnnData.concatenate(data_s10,data_s21,data_s2)
adatas.obs

sc.pp.filter_cells(adatas, min_genes=200)
sc.pp.filter_genes(adatas, min_cells=3)

adatas.var['mt'] = adatas.var_names.str.startswith('MT-')
adatas.var['mt']
sc.pp.calculate_qc_metrics(adatas, qc_vars=['mt'], 
    percent_top=None, log1p=False, inplace=True)
adatas.obs

sc.pl.violin(adatas, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)

adatas = adatas[adatas.obs.n_genes_by_counts > 500, :]
adatas = adatas[adatas.obs.pct_counts_mt < 20, :]

sc.pp.normalize_total(adatas, target_sum=1e4)
sc.pp.log1p(adatas)
sc.pp.highly_variable_genes(adatas, min_mean=0.0125, max_mean=3, min_disp=0.5)

adatas.raw = adatas
adatas = adatas[:, adatas.var.highly_variable]

sc.pp.regress_out(adatas, ['total_counts', 'pct_counts_mt'])

sc.pp.scale(adatas, max_value=10)
sc.pp.pca(adatas)
sc.pp.neighbors(adatas)
sc.tl.umap(adatas)
sc.pl.umap(adatas, color=['batch', 'orig'], palette=sc.pl.palettes.vega_20_scanpy)

import scanpy.external as sce
sce.pp.harmony_integrate(adatas, 'orig')
sc.pp.neighbors(adatas, use_rep="X_pca_harmony")
sc.tl.umap(adatas,init_pos='X_pca_harmony')
sc.pl.umap(adatas, color=['batch', 'orig'], legend_fontsize=8)


adatas.obsm

sc.tl.leiden(adatas )
sc.pl.umap(adatas, color=['leiden', 'CD3D' ])


sc.external.pp.bbknn(adatas, batch_key='batch')
sc.tl.umap(adatas)
sc.pl.umap(adatas, color=['batch', 'orig'])
adatas.obsm
sc.tl.leiden(adatas )

sc.pl.umap(adatas, color=['leiden', 'CD3D', 'NKG7'])


sc.tl.leiden(adatas )
sc.pl.umap(adatas, color=['leiden', 'CD3D', 'NKG7'])

sc.tl.rank_genes_groups(adatas, 'leiden', method='t-test')
sc.pl.rank_genes_groups(adatas, n_genes=25, sharey=False)
pd.DataFrame(adatas.uns['rank_genes_groups']['names']).head(5)


new_cluster_names = [
    '0', '1', '2: CD14+Monocytes','3 : epi', '4:endo',
    '5: NK', '6;T', '7;T', '8: Megakaryocytes','9',
    '10','11','12','13','14','15','16','17','18',
    '19','20','21','22','23','24']
adatas.rename_categories('leiden', new_cluster_names)


sc.pl.umap(adatas, color='leiden', legend_loc='on data', title='', frameon=False)

sc.pl.rank_genes_groups_tracksplot(adatas, n_genes=3)

import scanpy.external as sce
sce.pp.harmony_integrate(adatas, 'batch')

```



```bash
harmony

import scanpy as sc
import pandas as pd
import seaborn as sns
import os
sc.settings.verbosity = 1             # verbosity errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
sc.settings.set_figure_params(dpi=80, frameon=False, figsize=(3, 3), facecolor='white')

print(os.getcwd())
os.chdir('F://spe.lesson//k7')
## 更改完再检查一下有没有更改成功
print(os.getcwd())

data_s10=sc.read_10x_mtx('BC10\\')
data_s21=sc.read_10x_mtx('BC21\\')
data_s2=sc.read_10x_mtx('BC2\\')

data_s10.obs["orig"]='s10'
data_s21.obs["orig"]='s21'
data_s2.obs["orig"]='s2'

adatas = sc.AnnData.concatenate(data_s10,data_s21,data_s2)
#可以加一个batch_key='sample' 列名就改为sample
adatas1 = sc.AnnData.concatenate(data_s10,data_s21,data_s2,batch_key='sample')
adatas1.obs

sc.pp.filter_cells(adatas, min_genes=200)
sc.pp.filter_genes(adatas, min_cells=3)

adatas.var['mt'] = adatas.var_names.str.startswith('MT-')
adatas.var['mt']
sc.pp.calculate_qc_metrics(adatas, qc_vars=['mt'], 
    percent_top=None, log1p=False, inplace=True)
adatas.obs


sc.pl.violin(adatas, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)

adatas = adatas[adatas.obs.n_genes_by_counts > 500, :]
adatas = adatas[adatas.obs.pct_counts_mt < 20, :]

sc.pp.normalize_total(adatas, target_sum=1e4)
sc.pp.log1p(adatas)

sc.pp.highly_variable_genes(adatas, min_mean=0.0125, max_mean=3, min_disp=0.5)

adatas.raw = adatas

adatas = adatas[:, adatas.var.highly_variable]

sc.pp.regress_out(adatas, ['total_counts', 'pct_counts_mt'])


sc.pp.scale(adatas, max_value=10)
sc.tl.pca(adatas)

import scanpy.external as sce
sce.pp.harmony_integrate(adatas, 'orig')

sc.pp.neighbors(adatas, use_rep="X_pca_harmony")
sc.tl.umap(adatas,init_pos='X_pca_harmony')
sc.pl.umap(adatas, color=['batch', 'orig'], legend_fontsize=8)
adatas.obsm


sc.tl.leiden(adatas )
sc.pl.umap(adatas, color=['leiden', 'CD3D' ])
```





```bash
bbknn

import scanpy as sc
import pandas as pd
import seaborn as sns
import os
sc.settings.verbosity = 1             # verbosity errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
sc.settings.set_figure_params(dpi=80, frameon=False, figsize=(3, 3), facecolor='white')

print(os.getcwd())
os.chdir('F://spe.lesson//k7')
## 更改完再检查一下有没有更改成功
print(os.getcwd())

data_s10=sc.read_10x_mtx('BC10\\')
data_s21=sc.read_10x_mtx('BC21\\')
data_s2=sc.read_10x_mtx('BC2\\')

data_s10.obs["orig"]='s10'
data_s21.obs["orig"]='s21'
data_s2.obs["orig"]='s2'

adatas = sc.AnnData.concatenate(data_s10,data_s21,data_s2)
#可以加一个batch_key='sample' 列名就改为sample
adatas1 = sc.AnnData.concatenate(data_s10,data_s21,data_s2,batch_key='sample')
adatas1.obs

sc.pp.filter_cells(adatas, min_genes=200)
sc.pp.filter_genes(adatas, min_cells=3)

adatas.var['mt'] = adatas.var_names.str.startswith('MT-')
adatas.var['mt']
sc.pp.calculate_qc_metrics(adatas, qc_vars=['mt'], 
    percent_top=None, log1p=False, inplace=True)
adatas.obs


sc.pl.violin(adatas, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)

adatas = adatas[adatas.obs.n_genes_by_counts > 500, :]
adatas = adatas[adatas.obs.pct_counts_mt < 20, :]

sc.pp.normalize_total(adatas, target_sum=1e4)
sc.pp.log1p(adatas)

sc.pp.highly_variable_genes(adatas, min_mean=0.0125, max_mean=3, min_disp=0.5)

adatas.raw = adatas

adatas = adatas[:, adatas.var.highly_variable]

sc.pp.regress_out(adatas, ['total_counts', 'pct_counts_mt'])


sc.pp.scale(adatas, max_value=10)
sc.pp.pca(adatas)
sc.pp.neighbors(adatas)
sc.tl.umap(adatas)
sc.pl.umap(adatas, color=['batch', 'orig'], palette=sc.pl.palettes.vega_20_scanpy)

 
sc.external.pp.bbknn(adatas, batch_key='batch')
sc.tl.umap(adatas)
sc.pl.umap(adatas, color=['batch', 'orig'])
adatas.obsm
sc.tl.leiden(adatas )

sc.pl.umap(adatas, color=['leiden', 'CD3D', 'NKG7'])

sc.tl.rank_genes_groups(adatas, 'leiden', method='t-test')

sc.pl.rank_genes_groups(adatas, n_genes=25, sharey=False)

pd.DataFrame(adatas.uns['rank_genes_groups']['names']).head(5)


new_cluster_names = [
    '0', '1', '2: CD14+Monocytes','3 : epi', '4:endo',
    '5: NK', '6;T', '7;T', '8: Megakaryocytes','9',
    '10','11','12','13','14','15','16','17','18',
    '19','20','21','22','23','24']
adatas.rename_categories('leiden', new_cluster_names)


sc.pl.umap(adatas, color='leiden', legend_loc='on data', title='', frameon=False)

sc.pl.rank_genes_groups_tracksplot(adatas, n_genes=3)

```

