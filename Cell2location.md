# cell2location

```bash
创建conda环境，指定python=3.10.14
conda create -n cell2location python=3.10.14

激活环境
conda activate cell2location

先安装cell2location
pip install cell2location

加载包
python
import cell2location

exit()

再把个别依赖包的版本调整一下：
pip install --force-reinstall anndata==0.10.8
pip install --force-reinstall attrs==23.2.0
pip install --force-reinstall scipy==1.13.1
pip install --force-reinstall scvi-tools==1.1.5
pip install --force-reinstall NumPy==1.26.4

继续检查是否安装成功
python
import cell2location


(cell2location) C:\Users\admin>python
Python 3.10.14 | packaged by conda-forge | (main, Mar 20 2024, 12:40:08) [MSC v.1938 64 bit (AMD64)] on win32
Type "help", "copyright", "credits" or "license" for more information.
>>> import cell2location
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File "D:\miniconda\envs\cell2location\lib\site-packages\cell2location\__init__.py", line 9, in <module>
    from . import models
  File "D:\miniconda\envs\cell2location\lib\site-packages\cell2location\models\__init__.py", line 1, in <module>
    from ._cell2location_model import Cell2location
  File "D:\miniconda\envs\cell2location\lib\site-packages\cell2location\models\_cell2location_model.py", line 29, in <module>
    from cell2location.models.base._pyro_base_loc_module import Cell2locationBaseModule
  File "D:\miniconda\envs\cell2location\lib\site-packages\cell2location\models\base\_pyro_base_loc_module.py", line 6, in <module>
    from ._pyro_mixin import AutoGuideMixinModule, init_to_value
  File "D:\miniconda\envs\cell2location\lib\site-packages\cell2location\models\base\_pyro_mixin.py", line 20, in <module>
    from scvi.model._utils import parse_use_gpu_arg
ImportError: cannot import name 'parse_use_gpu_arg' from 'scvi.model._utils' (D:\miniconda\envs\cell2location\lib\site-packages\scvi\model\_utils.py)
>>>

我们需要找到路径修改源代码
D:\miniconda\envs\cell2location\lib\site-packages\scvi\model\_utils.py

拖到文件的最后，将下面代码复制进来
def parse_use_gpu_arg(
    use_gpu: Optional[Union[str, int, bool]] = None,
    return_device=True,
):
    """Parses the use_gpu arg in codebase.
    Returned gpus are is compatible with PytorchLightning's gpus arg.
    If return_device is True, will also return the device.
    Parameters
    ----------
    use_gpu
        Use default GPU if available (if None or True), or index of GPU to use (if int),
        or name of GPU (if str, e.g., `'cuda:0'`), or use CPU (if False).
    return_device
        If True, will return the torch.device of use_gpu.
    Returns
    -------
    Arguments for lightning trainer, including the accelerator (str), devices
    (int or sequence of int), and optionally the torch device.
    """
    # Support Apple silicon
    cuda_available = torch.cuda.is_available()
    # If using an older version of torch.
    try:
        mps_available = torch.backends.mps.is_available()
    except AttributeError:
        mps_available = False
    gpu_available = cuda_available
    lightning_devices = None
    if (use_gpu is None and not gpu_available) or (use_gpu is False):
        accelerator = "cpu"
        device = torch.device("cpu")
        lightning_devices = "auto"
    elif (use_gpu is None and gpu_available) or (use_gpu is True):
        current = torch.cuda.current_device() if cuda_available else "mps"
        if current != "mps":
            lightning_devices = [current]
            accelerator = "gpu"
        else:
            accelerator = "mps"
            lightning_devices = 1
        device = torch.device(current)
    # Also captures bool case
    elif isinstance(use_gpu, int):
        device = torch.device(use_gpu) if not mps_available else torch.device("mps")
        accelerator = "gpu" if not mps_available else "mps"
        lightning_devices = [use_gpu] if not mps_available else 1
    elif isinstance(use_gpu, str):
        device = torch.device(use_gpu)
        accelerator = "gpu"
        # changes "cuda:0" to "0,"
        lightning_devices = [int(use_gpu.split(":")[-1])]
    else:
        raise ValueError("use_gpu argument not understood.")

    if return_device:
        return accelerator, lightning_devices, device
    else:
        return accelerator, lightning_devices
        
      
最后安装 spyder-notebook
conda install spyder-notebook  -c https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge

pip3 install igraph
pip3 install leidenalg

打开在终端输入 spyder
spyder

```

```bash

import sys
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib as mpl
import cell2location
from cell2location.utils.filtering import filter_genes
import os

plt.rcParams['pdf.fonttype'] = 42

#定义存储路径
results_folder = 'results/'
ref_run_name = f'{results_folder}/reference_signatures'
run_name = f'{results_folder}/cell2location_map'

print(os.getcwd())
adata_vis = sc.read_h5ad("V1_Human_Lymph_Node")


adata_vis.obs.columns
df=adata_vis.obs
sc.pl.spatial(adata_vis, color= ['in_tissue', 'array_col'] )

adata_vis.obs['sample'] = list(adata_vis.uns['spatial'].keys())[0]

adata_vis.var_names

adata_vis.var['SYMBOL'] = adata_vis.var_names
adata_vis.var.set_index('gene_ids', drop=True, inplace=True)
adata_vis.var_names

adata_vis.var['MT_gene'] = [gene.startswith('MT-') for gene in adata_vis.var['SYMBOL']]
# 去除线粒体基因
adata_vis.obsm['MT'] = adata_vis[:, adata_vis.var['MT_gene'].values].X.toarray()
adata_vis = adata_vis[:, ~adata_vis.var['MT_gene'].values]


#adata_ref = sc.read(
#    f'./data/sc.h5ad',
#    backup_url='https://cell2location.cog.sanger.ac.uk/paper/integrated_lymphoid_organ_scrna/RegressionNBV4Torch_57covariates_73260cells_10237genes/sc.h5ad'
#)
# adata_ref.write('sc.h5ad')
# 或者先把数据下载下来再读入
adata_ref = sc.read_h5ad('sc.h5ad')

adata_ref.var['SYMBOL'] = adata_ref.var.index
adata_ref.var 

adata_ref.var.set_index('GeneID-2', drop=True, inplace=True)
adata_ref.var

# 删除不需要的数据
del adata_ref.raw

selected = filter_genes(adata_ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)

# 过滤
adata_ref = adata_ref[:, selected].copy()

# 设置Anndata对象，以便后续在cell2location模型中使用
cell2location.models.RegressionModel.setup_anndata(adata=adata_ref,
                        # 指定Anndata对象中用于区分不同批次或样本的键值
                        batch_key='Sample',
                        # 指定Anndata对象中用于区分不同细胞类型的键值
                        labels_key='Subset',
                        # 指定Anndata对象中用于建立细胞类型签名的其他类别型协变量
                        categorical_covariate_keys=['Method']
                       )


from cell2location.models import RegressionModel
mod = RegressionModel(adata_ref)
# view anndata_setup as a sanity check
mod.view_anndata_setup()

mod.train(max_epochs=10)

mod.plot_history(2)


# num_samples：指定要抽样的后验样本的数量
# batch_size：控制每次推断所用的数据批次的大小
adata_ref = mod.export_posterior(
    adata_ref, sample_kwargs={'num_samples': 500, 'batch_size': 2500}
)

# Save model
mod.save(f"{ref_run_name}", overwrite=True)

# Save anndata object with results
adata_file = f"{ref_run_name}/sc.h5ad"
adata_ref.write(adata_file)
adata_file


#计算后验分布的5％、50％和95％分位数，而不必使用来自分布的1000个样本（或任何其他分位数）。这样可以加快对大型数据集的应用速度，并且需要更少的内存 - 但是，不能使用这种方法计算后验均值和标准差。
adata_ref = mod.export_posterior(
    adata_ref, use_quantiles=True,
    # choose quantiles
    add_to_varm=["q05","q50", "q95", "q0001"],
    sample_kwargs={'batch_size': 2500}
)
# 查看模型的效果
mod.plot_QC(summary_name="q50")


adata_file = f"{ref_run_name}/sc.h5ad"
adata_ref = sc.read_h5ad(adata_file)
mod = cell2location.models.RegressionModel.load(f"{ref_run_name}", adata_ref)

if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
    inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
inf_aver.columns = adata_ref.uns['mod']['factor_names']
inf_aver.iloc[0:5, 0:5]


##  Cell2location: spatial mapping
intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
adata_vis = adata_vis[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()

# prepare anndata for cell2location model
cell2location.models.Cell2location.setup_anndata(adata=adata_vis, batch_key="sample")

# create and train the model
mod = cell2location.models.Cell2location(
    adata_vis, cell_state_df=inf_aver,
    # 设置每个空间spot有多少个细胞
N_cells_per_location=30,
    # hyperparameter controlling normalisation of
    # 控制RNA检测中实验内变异归一化的超参数
detection_alpha=20
)
mod.view_anndata_setup()


mod.train(max_epochs=100,
          # train using full data (batch_size=None)
          batch_size=None,
          # use all data points in training because
          # we need to estimate cell abundance at all locations
          train_size=1, use_gpu=True,
         )

mod.plot_history(10)
plt.legend(labels=['full data training']);


adata_vis = mod.export_posterior(
    adata_vis, sample_kwargs={'num_samples': 100, 'batch_size': mod.adata.n_obs}
)

# Save model
mod.save(f"{run_name}", overwrite=True)

# mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)

# Save anndata object with results
adata_file = f"{run_name}/sp.h5ad"
adata_vis.write(adata_file)
adata_file

mod.plot_QC()
fig = mod.plot_spatial_QC_across_batches()

adata_vis.obsm
adata_vis.obsm['q05_cell_abundance_w_sf']
pd.DataFrame(adata_vis.obsm['q05_cell_abundance_w_sf']).to_csv(f"{results_folder}/st_cell2location_res.csv"  )

adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']

from cell2location.utils import select_slide
slide = select_slide(adata_vis, 'V1_Human_Lymph_Node')

with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):

    sc.pl.spatial(slide, cmap='magma',
                  # show first 8 cell types
                  color=['B_Cycling', 'B_GC_LZ', 'T_CD4+_TfH_GC', 'FDC',
                         'B_naive', 'T_CD4+_naive', 'B_plasma', 'Endo'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  # limit color scale at 99.2% quantile of cell abundance
                  vmin=0, vmax='p99.2'
                 )


# Now we use cell2location plotter that allows showing multiple cell types in one panel
from cell2location.plt import plot_spatial

# select up to 6 clusters
clust_labels = ['T_CD4+_naive', 'B_naive', 'FDC']
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'V1_Human_Lymph_Node')

with mpl.rc_context({'figure.figsize': (15, 15)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
        show_img=True,
        # 'fast' (white background) or 'dark_background'
        style='fast',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=6,
        colorbar_position='right'
    )


# Now we use cell2location plotter that allows showing multiple cell types in one panel
from cell2location.plt import plot_spatial

# select up to 6 clusters
clust_labels = ['T_CD4+_naive', 'B_naive', 'FDC']
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

#slide = select_slide(adata_st, 'V1_Human_Lymph_Node')

with mpl.rc_context({'figure.figsize': (10, 10)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
        show_img=False,
        img_key='lowres',
        # 'fast' (white background) or 'dark_background'
        style='fast',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=6,
        colorbar_position='right'
    )

#####################################################
###  下游分析
adata_vis.obs.columns
# compute KNN using the cell2location output stored in adata.obsm
sc.pp.neighbors(adata_vis, use_rep='q05_cell_abundance_w_sf',
                n_neighbors = 15)

# Cluster spots into regions using scanpy
sc.tl.leiden(adata_vis, resolution=1.1)

# add region as categorical variable
adata_vis.obs["region_cluster"] = adata_vis.obs["leiden"].astype("category")

# compute UMAP using KNN graph based on the cell2location output
sc.tl.umap(adata_vis, min_dist = 0.3, spread = 1)

# show regions in UMAP coordinates
with mpl.rc_context({'axes.facecolor':  'white',
                     'figure.figsize': [8, 8]}):
    sc.pl.umap(adata_vis, color=['region_cluster'], size=30,
               color_map = 'RdPu', ncols = 2, legend_loc='on data',
               legend_fontsize=20)
    sc.pl.umap(adata_vis, color=['sample'], size=30,
               color_map = 'RdPu', ncols = 2,
               legend_fontsize=20)

# plot in spatial coordinates
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(adata_vis, color=['region_cluster'],
                  size=1.3, img_key='hires', alpha=0.5)


#########################################################################
####### NMF

from cell2location import run_colocation
res_dict, adata_vis = run_colocation(
    adata_vis,
    model_name='CoLocatedGroupsSklearnNMF',
    train_args={
      'n_fact': np.arange(11, 13), # IMPORTANT: use a wider range of the number of factors (5-30)
      'sample_name_col': 'sample', # columns in adata_vis.obs that identifies sample
      'n_restarts': 3 # number of training restarts
    },
    # the hyperparameters of NMF can be also adjusted:
    model_kwargs={'alpha': 0.01, 'init': 'random', "nmf_kwd_args": {"tol": 0.000001}},
    export_args={'path': f'{run_name}/CoLocatedComb/'}
)

# Here we plot the NMF weights (Same as saved to `cell_type_fractions_heatmap`)
res_dict['n_fact12']['mod'].plot_cell_type_loadings()



############################################################################
### Estimate cell-type specific expression of every gene
# Compute expected expression per cell type
expected_dict = mod.module.model.compute_expected_per_cell_type(
    mod.samples["post_sample_q05"], mod.adata_manager
)

# Add to anndata layers
for i, n in enumerate(mod.factor_names_):
    adata_vis.layers[n] = expected_dict['mu'][i]

# Save anndata object with results
adata_file = f"{run_name}/sp.h5ad"
adata_vis.write(adata_file)
adata_file

# list cell types and genes for plotting
ctypes = ['T_CD4+_TfH_GC', 'T_CD4+_naive', 'B_GC_LZ']
genes = ['CD3D', 'CR2']

with mpl.rc_context({'axes.facecolor':  'black'}):
    # select one slide
    slide = select_slide(adata_vis, 'V1_Human_Lymph_Node')

    from tutorial_utils import plot_genes_per_cell_type
    plot_genes_per_cell_type(slide, genes, ctypes);


def huage(slide, genes, ctypes):
    slide.var['SYMBOL'] = slide.var.index
    n_genes = len(genes)
    n_ctypes = len(ctypes)
    fig, axs = plt.subplots(
        nrows=n_genes, ncols=n_ctypes + 1, figsize=(4.5 * (n_ctypes + 1) + 2, 5 * n_genes + 1), squeeze=False
    )
    # axs = axs.reshape((n_genes, n_ctypes+1))

    # plots of every gene
    for j in range(n_genes):
        # limit color scale at 99.2% quantile of gene expression (computed across cell types)
        quantile_across_ct = np.array(
            [
                np.quantile(slide.layers[n][:, slide.var["SYMBOL"] == genes[j]].toarray(), 0.992)
                for n in slide.uns["mod"]["factor_names"]
            ]
        )
        quantile_across_ct = np.partition(quantile_across_ct.flatten(), -2)[-2]
        sc.pl.spatial(
            slide,
            cmap="magma",
            color=genes[j],
            # layer=ctypes[i],
            #gene_symbols="SYMBOL",
            ncols=4,
            size=1.3,
            #img_key="hires",
            # limit color scale at 99.2% quantile of gene expression
            vmin=0,
            vmax="p99.2",
            ax=axs[j, 0],
            show=False,
        )

        # plots of every cell type
        for i in range(n_ctypes):
            sc.pl.spatial(
                slide,
                cmap="magma",
                color=genes[j],
                layer=ctypes[i],
                gene_symbols="SYMBOL",
                ncols=4,
                size=1.3,
                #img_key="hires",
                # limit color scale at 99.2% quantile of gene expression
                vmin=0,
                vmax=quantile_across_ct,
                ax=axs[j, i + 1],
                show=False,
            )
            axs[j, i + 1].set_title(f"{genes[j]} {ctypes[i]}")

    return fig, axs





# list cell types and genes for plotting
ctypes = ['T_CD4+_TfH_GC', 'T_CD4+_naive', 'B_GC_LZ']
genes = ['CD3D', 'CR2']

with mpl.rc_context({'axes.facecolor':  'black'}):
    # select one slide
    slide = select_slide(adata_vis, 'V1_Human_Lymph_Node')
    huage(slide, genes,ctypes);
```


