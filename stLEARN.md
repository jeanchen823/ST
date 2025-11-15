# stLEARN



```bash
conda create -n stlearn.new  python=3.9  -c https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge
## 激活环境
conda activate stlearn.new


pip install -i https://pypi.tuna.tsinghua.edu.cn/simple -U stlearn


conda install NumPy=1.21.0
pip install --force-reinstall rpy2==3.5.1


conda install spyder-notebook  -c https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge


pip install --force-reinstall networkx==2.8.6
pip install --force-reinstall matplotlib==3.7.0
pip install --force-reinstall numba==0.57.1
pip install --force-reinstall scipy==1.10.0
pip install --force-reinstall numpy==1.21.0 --user





import stlearn as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
os.chdir('stlearn\\' )


data=st.Read10X("./data/BCBA")


data.var_names_make_unique()st.add.image(adata=data,
imgpath="./data/BCBA/spatial/tissue_hires_image.png",
library_id="V1_Breast_Cancer_Block_A_Section_1",visium=True)


st.pp.filter_genes(data, min_cells=3)
st.pp.normalize_total(data) # NOTE: no log1p


# Adding the label transfer results,  #
spot_mixtures = pd.read_csv(data_dir+'../../Brad/label_transfer_bc.csv', index_col=0, sep='\t')
labels = spot_mixtures.loc[:,'predicted.id'].values.astype(str)spot_mixtures = spot_mixtures.drop(['predicted.id','prediction.score.max'],                                   axis=1)
spot_mixtures.columns = [col.replace('prediction.score.', '')                         for col in spot_mixtures.columns]
# Note the format! #
print(labels)



print(spot_mixtures)


# Check is in correct order
print('Spot mixture order correct?: ',      np.all(spot_mixtures.index.values==data.obs_names.values)) # Check is in correct order



# NOTE: using the same key in data.obs & data.uns
data.obs['cell_type'] = labels # Adding the dominant cell type labels per spot
data.obs['cell_type'] = data.obs['cell_type'].astype('category')
data.uns['cell_type'] = spot_mixtures # Adding the cell type scores



st.pl.cluster_plot(data, use_label='cell_type')

lrs = st.tl.cci.load_lrs(['connectomeDB2020_lit'], species='human')
print(len(lrs))

# Running the analysis #
st.tl.cci.run(data, lrs,
					min_spots = 20, #Filter out any LR pairs with no scores for less than min_spots                  
					distance=None, # None defaults to spot+immediate neighbours; distance=0 for within-spot mode                  
					n_pairs=100, # Number of random pairs to generate; low as example, recommend ~10,000                  
					n_cpus=4, # Number of CPUs for parallel. If None, detects & use all available.                  )



# A dataframe detailing the LR pairs ranked by number of significant spots.
lr_info = data.uns['lr_summary']
print('\n', lr_info)


st.tl.cci.adj_pvals(data, correct_axis='spot',pval_adj_cutoff=0.05, adj_method='fdr_bh')


st.pl.lr_summary(data, n_top=500)
st.pl.lr_summary(data, n_top=50, figsize=(10,3))


st.pl.lr_diagnostics(data, figsize=(10,2.5))


st.pl.lr_n_spots(data, n_top=50, figsize=(11, 3), max_text=100)
st.pl.lr_n_spots(data, n_top=500, figsize=(11, 3), max_text=100)



options(BioC_mirror="https://mirrors.westlake.edu.cn/bioconductor")
BiocManager::install('clusterProfiler')
BiocManager::install("org.Hs.eg.db")
BiocManager::install("org.Mm.eg.db")
#检查是否安装成功
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(clusterProfiler)


r_path = "D:/R-4.3.1"
st.tl.cci.run_lr_go(data, r_path)


R[write to console]:  cannot open file
D:\miniconda\envs\stlearn.new\Lib\site-packages\stlearn\tools\microenv\cci/go.R
No such file or directory



r_path = "D:/R-4.3.1"
st.tl.cci.run_lr_go(data, r_path)



st.pl.lr_go(data, lr_text_fp={'weight': 'bold', 'size': 10}, rot=15,               figsize=(12,3.65), n_top=15, show=False)


```



```bash
best_lr = data.uns['lr_summary'].index.values[0]
stats = ['lr_scores', 'p_vals', 'p_adjs', '-log10(p_adjs)']
fig, axes = plt.subplots(ncols=len(stats), figsize=(16, 6))
for i, stat in enumerate(stats):
    st.pl.lr_result_plot(data, use_result=stat, use_lr=best_lr, show_color_bar=False, ax=axes[i])
    axes[i].set_title(f'{best_lr} {stat}')


fig, axes = plt.subplots(ncols=2, figsize=(8, 6))
st.pl.lr_result_plot(data, use_result='-log10(p_adjs)', use_lr=best_lr, show_color_bar=False, ax=axes[0])
st.pl.lr_result_plot(data, use_result='lr_sig_scores', use_lr=best_lr, show_color_bar=False, ax=axes[1])
axes[0].set_title(f'{best_lr} -log10(p_adjs)')
axes[1].set_title(f'{best_lr} lr_sig_scores')


st.pl.lr_plot(data, best_lr, inner_size_prop=0.1, outer_mode='binary', pt_scale=5,
              use_label=None, show_image=True, sig_spots=False)


st.pl.lr_plot(data, best_lr, outer_size_prop=1, outer_mode='binary', pt_scale=20,
              use_label=None, show_image=True, sig_spots=True)


st.pl.lr_plot(data, best_lr, inner_size_prop=0.04, middle_size_prop=0.07,
			outer_size_prop=0.4,outer_mode='continuous', pt_scale=60, use_label=None,
			show_image=True, sig_spots=False)



st.pl.lr_plot(data, best_lr, 
			inner_size_prop=0.04, middle_size_prop=0.07, outer_size_prop=0.4,
			outer_mode='continuous', pt_scale=60, use_label=None, show_image=True,
			sig_spots=True)


st.pl.lr_plot(data, best_lr, inner_size_prop=0.08, middle_size_prop=0.3, outer_size_prop=0.5,
              outer_mode='binary', pt_scale=50, show_image=True, arrow_width=10, arrow_head_width=10,
              sig_spots=True, show_arrows=True)

st.pl.lr_plot(data, best_lr, inner_size_prop=0.08, middle_size_prop=0.3, outer_size_prop=0.5,
              outer_mode='binary', pt_scale=150, use_label='cell_type', show_image=True, sig_spots=True)


st.tl.cci.run_cci(data, 'cell_type', min_spots=3, spot_mixtures=True, cell_prop_cutoff=0.2,
                  sig_spots=True, n_perms=100)

st.tl.cci.run_cci(data, 'cell_type', # Spot cell information either in data.obs or data.uns                  
						min_spots=3, # Minimum number of spots for LR to be tested.    							
						spot_mixtures=True, # If True will use the label transfer scores,                                      # so spots can have multiple cell types if score>cell_prop_cutoff                  
						cell_prop_cutoff=0.2, # Spot considered to have cell type if score>0.2                  
						sig_spots=True, # Only consider neighbourhoods of spots which had significant LR scores.                  
						n_perms=100 # Permutations of cell information to get background, recommend ~1000                 )


st.pl.cci_check(data, 'cell_type')

pos_1 = st.pl.ccinet_plot(data, 'cell_type', return_pos=True)

# Just examining the cell type interactions between selected pairs
# lrs = data.uns['lr_summary'].index.values[0:3]
for best_lr in lrs[0:3]:
    st.pl.ccinet_plot(data, 'cell_type', best_lr, min_counts=2, figsize=(10, 7.5), pos=pos_1)


st.pl.lr_chord_plot(data, 'cell_type')



for lr in lrs:
    st.pl.lr_chord_plot(data, 'cell_type', lr)


st.pl.lr_cci_map(data, 'cell_type', lrs=None, min_total=100, figsize=(20, 4))

st.pl.lr_cci_map(data, 'cell_type', lrs=lrs, min_total=100, figsize=(20, 4))

st.pl.cci_map(data, 'cell_type')

lrs = data.uns['lr_summary'].index.values[0:3]
for lr in lrs[0:3]:
    st.pl.cci_map(data, 'cell_type', lr)



best_lr = lrs[0]

# This will plot with simple black arrows
st.pl.lr_plot(data, best_lr, outer_size_prop=1, outer_mode=None, pt_scale=40, use_label='cell_type',
              show_arrows=True, show_image=True, sig_spots=False, sig_cci=True, arrow_head_width=4,
              arrow_width=1, cell_alpha=0.8)

# This will colour the spot by the mean LR expression in the spots connected by arrows
st.pl.lr_plot(data, best_lr, outer_size_prop=1, outer_mode=None, pt_scale=10, use_label='cell_type',
              show_arrows=True, show_image=True, sig_spots=False, sig_cci=True, arrow_head_width=4,
              arrow_width=2, arrow_cmap='YlOrRd', arrow_vmax=1.5)

```



```bash
import stlearn as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os



data = st.Read10X("stlearn\\BCBA")

data.var_names_make_unique()
st.add.image(
    adata=data,
    imgpath="stlearn\\BCBA/spatial/tissue_hires_image.png",
    library_id="V1_Breast_Cancer_Block_A_Section_1",
    visium=True
)


st.pp.filter_genes(data, min_cells=3)
st.pp.normalize_total(data)  # NOTE: no log1p

# Adding the label transfer results
spot_mixtures = pd.read_csv('stlearn\\label_transfer_bc.csv', index_col=0, sep='\t')
labels = spot_mixtures.loc[:, 'predicted.id'].values.astype(str)
spot_mixtures = spot_mixtures.drop(['predicted.id', 'prediction.score.max'], axis=1)
spot_mixtures.columns = [col.replace('prediction.score.', '') for col in spot_mixtures.columns]
print(labels)

print(spot_mixtures)
print('Spot mixture order correct?: ', np.all(spot_mixtures.index.values == data.obs_names.values))  # Check if in correct order

# NOTE: using the same key in data.obs & data.uns
data.obs['cell_type'] = labels  # Adding the dominant cell type labels per spot
data.obs['cell_type'] = data.obs['cell_type'].astype('category')
data.uns['cell_type'] = spot_mixtures  # Adding the cell type scores

# Cluster plot
st.pl.cluster_plot(data, use_label='cell_type')

# Ligand-receptor analysis
lrs = st.tl.cci.load_lrs(['connectomeDB2020_lit'], species='human')
print(len(lrs))

# Run the ligand-receptor interaction analysis
st.tl.cci.run(
    data, 
    lrs,
    min_spots=20,  # Filter out any LR pairs with no scores for less than min_spots
    distance=None,  # None defaults to spot+immediate neighbours; distance=0 for within-spot mode
    n_pairs=100,  # Number of random pairs to generate; low as example, recommend ~10,000
    n_cpus=4,  # Number of CPUs for parallel. If None, detects & use all available.
)

# View the LR summary
lr_info = data.uns['lr_summary']  # A dataframe detailing the LR pairs ranked by number of significant spots
print('\n', lr_info)

# Adjust p-values
st.tl.cci.adj_pvals(data, correct_axis='spot', pval_adj_cutoff=0.05, adj_method='fdr_bh')

# Ligand-receptor summary plots
st.pl.lr_summary(data, n_top=500)
st.pl.lr_summary(data, n_top=50, figsize=(10, 3))

# Diagnostics for ligand-receptor interactions
st.pl.lr_diagnostics(data, figsize=(10, 2.5))

# Ligand-receptor interactions for spots
st.pl.lr_n_spots(data, n_top=50, figsize=(11, 3), max_text=100)
st.pl.lr_n_spots(data, n_top=500, figsize=(11, 3), max_text=100)

# Run GO enrichment analysis for ligand-receptor interactions
r_path = "D:/R-4.3.1"
st.tl.cci.run_lr_go(data, r_path)

# Ligand-receptor GO enrichment plots
st.pl.lr_go(data, lr_text_fp={'weight': 'bold', 'size': 10}, rot=15, figsize=(12, 3.65), n_top=15, show=False)

# Plot results for the best ligand-receptor interaction
best_lr = data.uns['lr_summary'].index.values[0]
stats = ['lr_scores', 'p_vals', 'p_adjs', '-log10(p_adjs)']
fig, axes = plt.subplots(ncols=len(stats), figsize=(16, 6))

for i, stat in enumerate(stats):
    st.pl.lr_result_plot(data, use_result=stat, use_lr=best_lr, show_color_bar=False, ax=axes[i])
    axes[i].set_title(f'{best_lr} {stat}')

# Further visualizations for ligand-receptor results
fig, axes = plt.subplots(ncols=2, figsize=(8, 6))
st.pl.lr_result_plot(data, use_result='-log10(p_adjs)', use_lr=best_lr, show_color_bar=False, ax=axes[0])
st.pl.lr_result_plot(data, use_result='lr_sig_scores', use_lr=best_lr, show_color_bar=False, ax=axes[1])
axes[0].set_title(f'{best_lr} -log10(p_adjs)')
axes[1].set_title(f'{best_lr} lr_sig_scores')

# Plot ligand-receptor interactions with arrows
st.pl.lr_plot(data, best_lr, inner_size_prop=0.1, outer_mode='binary', pt_scale=5, use_label=None, show_image=True, sig_spots=False)

st.pl.lr_plot(data, best_lr, outer_size_prop=1, outer_mode='binary', pt_scale=20, use_label=None, show_image=True, sig_spots=True)

# Plot different modes of ligand-receptor interactions
st.pl.lr_plot(data, best_lr, inner_size_prop=0.04, middle_size_prop=.07, outer_size_prop=.4, outer_mode='continuous', pt_scale=60, use_label=None, show_image=True, sig_spots=False)

st.pl.lr_plot(data, best_lr, inner_size_prop=0.04, middle_size_prop=.07, outer_size_prop=.4, outer_mode='continuous', pt_scale=60, use_label=None, show_image=True, sig_spots=True)

# Plot ligand-receptor interactions with arrows and labels
st.pl.lr_plot(data, best_lr, inner_size_prop=0.08, middle_size_prop=.3, outer_size_prop=.5, outer_mode='binary', pt_scale=50, show_image=True, arrow_width=10, arrow_head_width=10, sig_spots=True, show_arrows=True)

st.pl.lr_plot(data, best_lr, inner_size_prop=0.08, middle_size_prop=.3, outer_size_prop=.5, outer_mode='binary', pt_scale=150, use_label='cell_type', show_image=True, sig_spots=True)

# Cell-cell interaction analysis
st.tl.cci.run_cci(
    data, 
    'cell_type',  # Spot cell information either in data.obs or data.uns
    min_spots=3,  # Minimum number of spots for LR to be tested
    spot_mixtures=True,  # If True will use the label transfer scores
    cell_prop_cutoff=0.2,  # Spot considered to have cell type if score > 0.2
    sig_spots=True,  # Only consider neighbourhoods of spots which had significant LR scores
    n_perms=100  # Permutations of cell information to get background, recommend ~1000
)

# Check the results of cell-cell interaction
st.pl.cci_check(data, 'cell_type')

# Visualize the results for selected ligand-receptor pairs
pos_1 = st.pl.ccinet_plot(data, 'cell_type', return_pos=True)
lrs = data.uns['lr_summary'].index.values[0:3]

for best_lr in lrs[0:3]:
    st.pl.ccinet_plot(data, 'cell_type', best_lr, min_counts=2, figsize=(10, 7.5), pos=pos_1)

# Chord plot for ligand-receptor interactions
st.pl.lr_chord_plot(data, 'cell_type')


for lr in lrs:
    st.pl.lr_chord_plot(data, 'cell_type', lr)

# Map the ligand-receptor interactions on the spatial data
st.pl.lr_cci_map(data, 'cell_type', lrs=None, min_total=100, figsize=(20, 4))

# Map the cell-cell interactions on the spatial data
st.pl.cci_map(data, 'cell_type')
lrs = data.uns['lr_summary'].index.values[0:3]
for lr in lrs[0:3]:
    st.pl.cci_map(data, 'cell_type', lr)

# Final visualization with arrows and ligand-receptor scores
best_lr = lrs[0]

st.pl.lr_plot(data, best_lr, outer_size_prop=1, outer_mode=None, pt_scale=40, use_label='cell_type', show_arrows=True, show_image=True, sig_spots=False, sig_cci=True, arrow_head_width=4, arrow_width=1, cell_alpha=.8)

st.pl.lr_plot(data, best_lr, outer_size_prop=1, outer_mode=None, pt_scale=10, use_label='cell_type', show_arrows=True, show_image=True, sig_spots=False, sig_cci=True, arrow_head_width=4, arrow_width=2, arrow_cmap='YlOrRd', arrow_vmax=1.5)

```



# 