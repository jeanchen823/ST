# Stereo



```bash
官网：https://stereopy.readthedocs.io/en/latest/Tutorials/CellBin_Clustering.html

conda create -n  stereo  python=3.8  -c https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge
conda activate stereo

pip install stereopy  -i https://mirrors.aliyun.com/pypi/simple/

conda install spyder-kernels=2.5  -c https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge
conda install spyder-notebook  -c https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge


import stereo as st
import warnings
import os

warnings.filterwarnings('ignore')


st.__version__


data_path = './Demo_MouseBrain/SS200000135TL_D1.tissue.gef'
st.io.read_gef_info(data_path)


data = st.io.read_gef(file_path=data_path, bin_size=50)

data

data.tl.cal_qc()
data.plt.violin()
data.plt.spatial_scatter()


data.plt.genes_count()

data.tl.filter_cells(
        min_counts=20,
        min_genes=3,
        pct_counts_mt=5,
        inplace=True)
data

data.tl.raw_checkpoint()
data.tl.raw


data.tl.normalize_total(target_sum=10000)
data.tl.log1p()

data.tl.highly_variable_genes(
        min_mean=0.0125,
        max_mean=3,
        min_disp=0.5,
        n_top_genes=2000,
        res_key='highly_variable_genes')

data.plt.highly_variable_genes(res_key='highly_variable_genes')

data.tl.scale(max_value=10, zero_center=True)


data.tl.pca(
        use_highly_genes=False,
        n_pcs=30,
        res_key='pca'
        )

data.tl.key_record

data.plt.elbow(pca_res_key='pca')


data.tl.neighbors(
        pca_res_key='pca',
        n_pcs=30,
        res_key='neighbors'
        )

# compute spatial neighbors
data.tl.spatial_neighbors(
        neighbors_res_key='neighbors',
        res_key='spatial_neighbors'
        )


data.tl.umap(pca_res_key='pca', neighbors_res_key='neighbors', res_key='umap')

data.tl.leiden(neighbors_res_key='neighbors', res_key='leiden')
data.plt.cluster_scatter(res_key='leiden')

data.plt.cluster_scatter(res_key='leiden', groups=['1', '2'])


data.plt.umap(res_key='umap', cluster_key='leiden')

data.tl.leiden(neighbors_res_key='spatial_neighbors', res_key='spatial_leiden')
data.plt.cluster_scatter(res_key='spatial_leiden')


data.tl.louvain(neighbors_res_key='neighbors', res_key='louvain')
data.plt.cluster_scatter(res_key='louvain')

data.tl.phenograph(phenograph_k=30, pca_res_key='pca', res_key='phenograph')
data.plt.cluster_scatter(res_key='phenograph')

data.tl.find_marker_genes(
        cluster_res_key='leiden',
        method='t_test',
        use_highly_genes=False,
        use_raw=True
        )

data.plt.marker_genes_text(
        res_key='marker_genes',
        markers_num=10,
        sort_key='scores'
        )

data.plt.marker_genes_scatter(res_key='marker_genes', markers_num=5)

data.plt.marker_gene_volcano(group_name='2.vs.rest', vlines=False)

data.tl.filter_marker_genes(
    marker_genes_res_key='marker_genes',
    min_fold_change=1,
    min_in_group_fraction=0.25,
    max_out_group_fraction=0.5,
    res_key='marker_genes_filtered'
)


annotation_dict = {
    '1':'a', '2':'b',
    '3':'c', '4':'d',
    '5':'e', '6':'f',
    '7':'g', '8':'h',
    '9':'i', '10':'j',
    '11':'k', '12': 'l',
    '13': 'm', '14': 'n',
    '15': 'o', '16': 'p',
    '17': 'q', '18': 'r',
    '19': 's', '20': 't',
    '21': 'u', '22': 'v',
    '23': 'w', '24': 'x',
    '25': 'y', '26': 'z'
    }
data.tl.annotation(
        annotation_information=annotation_dict,
        cluster_res_key='leiden',
        res_key='anno_leiden'
        )

data.plt.cluster_scatter(res_key='anno_leiden')
```

