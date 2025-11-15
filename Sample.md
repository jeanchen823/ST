# Sample

```bash

#加载包
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np
import pooch
from scipy.sparse import csr_matrix
from scipy.io import mmwrite

#初始化设置
sc.logging.print_versions()
sc.set_figure_params(facecolor="white", figsize=(10, 10))
sc.settings.verbosity = 3

#设置路径
print(os.getcwd())
os.chdir('')
## 更改完再检查一下有没有更改成功
print(os.getcwd())
#声明h5ad用于存储分析结果：
results_file = 'smaple1.h5ad'

#读取数据
adata = sc.read_10x_mtx(
    '', # `.mtx`文件所在的目录
    var_names='gene_symbols', # 用 gene 作为var
    cache=True)# 开启缓存读写

#消除重复的数据
adata.var_names_make_unique()


#计算每一个基因在所有细胞中的平均表达量，并绘制了平均表达量前10的基因箱型图。箱型图是对每个基因在所有细胞中表达量分布的更详细描述。
sc.pl.highest_expr_genes(adata, n_top=10)

##质量控制
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

adata.var['mt'] = adata.var_names.str.startswith('MT-')
adata.var['mt']
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'],
                           percent_top=None, log1p=False, inplace=True)
adata.obs


##小提琴图看看质控效果
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)

##散点图查看相关性
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')

##进行过滤
adata = adata[adata.obs.n_genes_by_counts > 500, :]
adata = adata[adata.obs.pct_counts_mt < 20, :]

##归一化

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

##找高变基因
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata)

##在进一步筛选数据前，我们将当前adata保存到adata的元素raw下：
adata.raw = adata

##提取高变基因的表达矩阵
adata = adata[:, adata.var.highly_variable]

##我们对两个注释进行回归，即：total_counts ：每个细胞的基因总计数（总表达量）；
#pct_counts_mt ：每个细胞中，线粒体基因表达量占该细胞所有基因表达量的百分比；
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])

##按零均值单位方差标准化数据，并剪裁值超过标准差10的细胞。
sc.pp.scale(adata, max_value=10)


##主成分分析
sc.tl.pca(adata, svd_solver='arpack')
##检查单个PC（主成分）对数据总方差的贡献
sc.pl.pca_variance_ratio(adata, log=True)

adata.write(results_file)

##计算细胞邻域
##使用数据矩阵的PCA表示来计算细胞的邻域图
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

##对neighborhood graph进行embedding
sc.tl.umap(adata)
sc.pl.umap(adata, color=['CD3D', 'NKG7', 'MS4A1','EPCAM'])
sc.pl.umap(adata, color=['CD3D', 'NKG7', 'MS4A1'], use_raw=False)

##对neighborhood graph进行聚类
sc.tl.leiden(adata)
##用umap对leiden的聚类结果可视化：
sc.pl.umap(adata, color=["leiden"])
sc.pl.umap(adata, color=['leiden', 'CD3D', 'NKG7'])
adata.write(results_file)




##找到marker基因
#计算每个簇中存在较高程度差异的基因的排名。默认情况下，使用 adata.raw 进行计算，在计算差异时，最简单最快的方式是 t-检验 ：
sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)

sc.tl.rank_genes_groups(adata, 'leiden', method='logreg')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
#以表格形式查看
#相当于把cluster的marker变成了dataframe
pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(5)




#显示差异基因的名称和差异度：
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
pd.DataFrame(
    {group + '_' + key[:1]: result[key][group]
     for group in groups for key in ['names', 'pvals']}).head(5)


##marker基因展示
from matplotlib.pyplot import rc_context
with rc_context({'figure.figsize': (9, 1.5)}):
    sc.pl.rank_genes_groups_violin(adata, n_genes=20, jitter=False)

sc.pl.violin(adata, ['CD3D', 'NKG7' ], groupby='leiden')


##细胞类型注释
new_cluster_names = [
    '0', '1', '2: CD14+Monocytes','3 : epi', '4:endo',
    '5: NK', '6;T', '7;T', '8: Megakaryocytes','9',
    '10','11','12','13','14','15','16','17','18',
    '19','20','21','22','23','24']
adata.rename_categories('leiden', new_cluster_names)


sc.pl.umap(adata, color='leiden', legend_loc='on data', title='', frameon=False)
sc.pl.rank_genes_groups_tracksplot(adata, n_genes=3)


```



```bash
多个样本读取
函数的帮助文档：https://scanpy.readthedocs.io/en/latest/generated/scanpy.read_10x_mtx.html

#导入scanpy库
import scanpy as sc
#读取tab文件
data=sc.read_10x_mtx('GSM6567952')
#查看数据
data

函数的帮助文档：https://scanpy.readthedocs.io/en/latest/generated/scanpy.read_10x_h5.html

#导入scanpy库
import scanpy as sc
#读取h5文件
data1=sc.read_10x_h5('GSM5344021_WT_tumor_filtered_feature_bc_matrix.h5')
data1


https://scanpy.readthedocs.io/en/latest/generated/scanpy.read_h5ad.html

data2=sc.read_h5ad('GSM4648564_adipose_raw_counts.h5ad')
data2


meta2=pd.read_csv('GSM4648564_adipose_celltypes.csv',header=None)

meta2.rename(columns={0: 'bc',1: 'celltype'}, inplace=True)

meta2.new=meta2.set_index([ 'bc'],inplace=False)

data2.var = meta2_new
data2.var


https://scanpy.readthedocs.io/en/latest/generated/scanpy.read_csv.html


data4=sc.read_csv('GSE130148_raw_counts.csv')
data4


https://scanpy.readthedocs.io/en/latest/generated/scanpy.read_text.html#scanpy.read_text

data5=sc.read_text('GSE116481_all_samples_raw_counts_matrix.txt')
data5




总结

ata=sc.read_10x_mtx('BC21\\')
data

data1=sc.read_10x_h5('GSM5344021_WT_tumor_filtered_feature_bc_matrix.h5')
data1

data2=sc.read_h5ad('GSM4648564_adipose_raw_counts.h5ad')
data2

data2.var

meta2=pd.read_csv('GSM4648564_adipose_celltypes.csv',header=None)


meta2.rename(columns={0: 'bc',1: 'celltype'}, inplace=True)

meta2_new=meta2.set_index([ 'bc'],inplace=False)

data2.var = meta2_new
data2.var

data4=sc.read_csv('GSE130148_raw_counts.csv')
data4

data5=sc.read_text('GSE116481_all_samples_raw_counts_matrix.txt')
data5

```


