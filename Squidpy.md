# Squidpy



```bash

https://squidpy.readthedocs.io/en/stable/notebooks/tutorials/tutorial_read_spatial.html


pip install squidpy

conda install -c conda-forge squidpy


import squidpy as sq

```



```bash

from numpy.random import default_rng
import matplotlib.pyplot as plt
import scanpy as sc
import squidpy as sq
from anndata import AnnData
sc.logging.print_header()
print(f"squidpy=={sq.__version__}")

rng = default_rng(42)
counts = rng.integers(0, 15, size=(10, 100))  # feature matrix
coordinates = rng.uniform(0, 10, size=(10, 2))  # spatial coordinates
image = rng.uniform(0, 1, size=(10, 10, 3))  # image


adata = AnnData(counts, obsm={"spatial": coordinates})

sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata)
adata

sq.pl.spatial_scatter(adata, shape=None, color="leiden", size=50)



sq.gr.spatial_neighbors(adata, radius=3.0)
sq.pl.spatial_scatter(adata,color="leiden",
    connectivity_key="spatial_connectivities",
    edges_color="black",shape=None,
    edges_width=1,size=3000,)

plt.imshow(image)


spatial_key = "spatial"
library_id = "tissue42"
adata.uns[spatial_key] = {library_id: {}}
adata.uns[spatial_key][library_id]["images"] = {}
adata.uns[spatial_key][library_id]["images"] = {"hires": image}
adata.uns[spatial_key][library_id]["scalefactors"] = {
    "tissue_hires_scalef": 1,
    "spot_diameter_fullres": 0.5,}


sq.pl.spatial_scatter(adata, color="leiden")

adata.uns[spatial_key][library_id]["scalefactors"] = {
    "tissue_hires_scalef": 0.5,
    "spot_diameter_fullres": 0.5,}
sq.pl.spatial_scatter(adata, color="leiden", size=2)


img = sq.im.ImageContainer(image)
img.show()


非常重要----举个例子：
import squidpy as sq
import scanpy as sc
import anndata as ad
import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt

os.chdir('F://spe.lesson//k9')
## 更改完再检查一下有没有更改成功
print(os.getcwd())

#adata = sq.datasets.mibitof()
#adata.write("mibitof.h5ad")
adata=sc.read('mibitof.h5ad')
adata


df=adata.obs

sc.pl.spatial(adata[adata.obs["library_id"] == 'point16'], color="Cluster",
        library_id='point16', title='point16', frameon=False,show=False)
        

fig, ax = plt.subplots(ncols=3, nrows=1, figsize=(14, 3))
i = 0
for library_id in adata.uns["spatial"].keys():
    sc.pl.spatial(
        adata[adata.obs["library_id"] == library_id], color="Cluster",
        library_id=library_id, title=library_id, frameon=False, ax=ax[i], show=False
    )
    i += 1


sq.pl.spatial_scatter(adata,library_key= 'point16' ,color="Cluster" )


adata.uns.keys()

adata.uns['spatial'].keys()

adata.uns['spatial']['point23'].keys()
```



```bash

https://squidpy.readthedocs.io/en/stable/notebooks/tutorials/tutorial_image_container.html

import numpy as np
import squidpy as sq
arr = np.ones((100, 100, 3))
arr[40:60, 40:60] = [0, 0.7, 1]

print(arr.shape)
img = sq.im.ImageContainer(arr, layer="img1")
img

arr_seg = np.zeros((100, 100))
arr_seg[40:60, 40:60] = 1
img.add_img(arr_seg, layer="seg1")
img
img["seg2"] = arr_seg
img
print(list(img))
img["img1"]
img.rename("seg2", "new-name")


img.show(layer="img1")



crop1 = img.crop_corner(30, 40, size=(30, 30), scale=1)
crop1.show(layer="img1")

crop2 = crop1.crop_corner(0, 0, size=(40, 40), scale=0.5)
crop2.show(layer="img1")


print(crop1.data.attrs)
print(crop2.data.attrs)

sq.im.ImageContainer.uncrop([crop1], shape=img.shape).show(layer="img1")
sq.im.ImageContainer.uncrop([crop2], shape=(50, 50)).show(layer="img1")

img.data

img_on_disk = sq.datasets.visium_hne_image()
print(type(img_on_disk["image"].data))

img_on_disk.compute()
print(type(img_on_disk["image"].data))
```



```bash

import squidpy as sq
import scanpy as sc
import anndata as ad
import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt

## 更改完再检查一下有没有更改成功
print(os.getcwd())

#adata = sq.datasets.mibitof()
#adata.write("mibitof.h5ad")
adata=sc.read('mibitof.h5ad')
adata

##########################################################################
df=adata.obs

sq.pl.spatial_scatter(adata,library_key= 'point16' ,color="Cluster" )
sq.pl.spatial_scatter(adata,library_key= 'point23' ,color="Cluster" )


sq.pl.spatial_scatter(adata,  color="Cluster" )

sc.pl.spatial(adata[adata.obs["library_id"] == 'point16'], color="Cluster",
        library_id='point16', title='point16', frameon=False,show=False)

fig, ax = plt.subplots(ncols=3, nrows=1, figsize=(14, 3))
i = 0
for library_id in adata.uns["spatial"].keys():
    sc.pl.spatial(
        adata[adata.obs["library_id"] == library_id], color="Cluster",
        library_id=library_id, title=library_id, frameon=False, ax=ax[i], show=False
    )
    i += 1


adata.uns.keys()
adata.uns['spatial'].keys()
adata.uns['spatial']['point23'].keys()
adata.uns['spatial']['point23']['images'].keys()


adata.var_names
adata.uns.keys()
adata.uns['spatial']['point8']['images'].keys()

df1=adata.obs
adata.obs.library_id
fig, ax = plt.subplots(ncols=3, nrows=1, figsize=(14, 3))
i = 0
for library_id in adata.uns["spatial"].keys():
    sc.pl.spatial(
        adata[adata.obs["library_id"] == library_id], color="Cluster",
        library_id=library_id, title=library_id, frameon=False, ax=ax[i], show=False
    )
    i += 1


plt.figure(figsize=(8, 8))
sc.pl.spatial(
        adata[adata.obs["library_id"] == 'point8'], library_id='point8', frameon=False, cmap='Spectral_r', size=0, ax=plt.gca(), img_key='segmentation'
    )
##########################################################################




img =sq.im.ImageContainer(adata.uns["spatial"]['point16']["images"]["hires"], 
                          library_id='point16', layer='image')


img.add_img(adata.uns["spatial"]['point16']["images"]["segmentation"],
            library_id='point16', layer="segmentation")
img["segmentation"].attrs["segmentation"] = True


img.show(layer='image')
## 细胞分割之后的数据
img.show(layer='segmentation')
## 颜色区分
# 145_CD45 - a immune cell marker (cyan).
# 174_CK - a tumor marker (magenta).
# 113_vimentin - a mesenchymal cell marker (yellow).

img.show("image", segmentation_layer="segmentation", segmentation_alpha=0.5)
img.show("image", segmentation_layer="segmentation", segmentation_alpha=0.8)

for i, cmap in zip(range(3), ['Reds', 'Greens', 'Blues']):
    img.show(layer='image', channel=i, cmap=cmap)


img
img.crop_center(y=0.5, x=0.5, radius=200).show(layer='image')
img.crop_center(y=0.5, x=0.5, radius=100).show(layer='image')

img.crop_corner(y=0.2, x=0.3, size=400).show(layer='image')
img.crop_corner(y=0.2, x=0.3, size=200).show(layer='image')

subset_img= img.crop_corner(y=0.2, x=0.3, size=200)
img, subset_img

adata_crop = subset_img.subset(adata)
adata_crop
adata


sc.pl.spatial(
        adata_crop[adata_crop.obs["library_id"] == 'point16'], color=['CD45', 'CD3'],
        library_id='point16', frameon=False, cmap='Spectral_r', size=1.3)
        
        

import squidpy as sq
import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

## 更改完再检查一下有没有更改成功
print(os.getcwd())
## 加载数据
#adata = sq.datasets.mibitof()
#adata.write("mibitof.h5ad")
adata=sc.read('mibitof.h5ad')
adata

imgs = []
for library_id in adata.uns["spatial"].keys():
    img = sq.im.ImageContainer(adata.uns["spatial"][library_id]["images"]["hires"], library_id=library_id, layer='image')
    img.add_img(adata.uns["spatial"][library_id]["images"]["segmentation"], library_id=library_id, layer="segmentation")
    img["segmentation"].attrs["segmentation"] = True
    imgs.append(img)
img = sq.im.ImageContainer.concat(imgs)

```



```bash
https://squidpy.readthedocs.io/en/latest/notebooks/examples/image/compute_features.html


import squidpy as sq
import scanpy as sc
import anndata as ad
import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt

print(os.getcwd())

img = sq.datasets.visium_hne_image_crop()
adata = sq.datasets.visium_hne_adata_crop()
img
img['image']
fig, axes = plt.subplots(1, 3, figsize=(5, 5))
for i in range(3):
    img.show(layer='image', ax=axes[i], channel=i, cmap='gray')

import seaborn as sns

#  - `'smooth'` - :func:`skimage.filters.gaussian`
#  - `'gray'` - :func:`skimage.color.rgb2gray`
crop = img.crop_corner(0, 0, size=1000)
sq.im.process(crop, layer="image", method="smooth", sigma=4)  # smooth, gray


fig, axes = plt.subplots(1, 3, figsize=(15, 4))
crop.show("image_smooth", cmap="gray", ax=axes[0])
axes[1].imshow(crop["image_smooth"][:, :, 0, 0] < 90)
_ = sns.histplot(np.array(crop["image_smooth"]).flatten(), bins=50, ax=axes[2])
plt.tight_layout()

sq.im.calculate_image_features(adata, img, features="summary", key_added="summary_features",
                               show_progress_bar=True, n_jobs=3)

adata.obsm["summary_features"].head()

sq.im.calculate_image_features(adata,img,
    features="summary",features_kwargs={
        "summary": {
            "quantiles": [0.2],
            "channels": [0, 1],}
    },key_added="summary_features_2",
    mask_circle=True,show_progress_bar=False,)

adata.obsm["summary_features_2"].head()
sc.pl.spatial(
    sq.pl.extract(adata, "summary_features_2"),
    color=["summary_ch-0_quantile-0.2", "summary_ch-1_quantile-0.2"],)

sq.im.calculate_image_features(adata,img,
    features="histogram",
    features_kwargs={"histogram": {"bins": 3, "channels": [0, 1]}},
    key_added="histogram_features",)

adata.obsm['histogram_features']

sc.pl.spatial(sq.pl.extract(adata, "histogram_features"),
    color=[None, "histogram_ch-0_bin-0", "histogram_ch-0_bin-1", "histogram_ch-0_bin-2"],
    bw=True)


features_kwargs={"texture": {
            "props": ["contrast"],
            "channels": [0],
            "distances": 1,
            "angles": (0, np.pi / 4) }
}

sq.im.calculate_image_features(adata,img,
    features="texture",key_added="texture_features",
    spot_scale=2,show_progress_bar=False,
    features_kwargs = features_kwargs)

adata.obsm['texture_features']
sc.pl.spatial(
    sq.pl.extract(adata, "texture_features"),
    color=[None, "texture_ch-0_contrast_dist-1_angle-0.00", "texture_ch-0_contrast_dist-1_angle-0.79"],
    bw=True,)

crop = img.crop_corner(0, 0, size=1000)
adata_crop = crop.subset(adata)
# smooth image
sq.im.process(crop, layer="image", method="smooth", sigma=4)

# plot the result
fig, axes = plt.subplots(1, 2)
for layer, ax in zip(["image", "image_smooth"], axes):
    crop.show(layer, ax=ax)
    ax.set_title(layer)

fig, axes = plt.subplots(1, 3, figsize=(15, 4))
crop.show("image_smooth", cmap="gray", ax=axes[0])
axes[1].imshow(crop["image_smooth"][:, :, 0, 0] < 90)
_ = sns.histplot(np.array(crop["image_smooth"]).flatten(), bins=50, ax=axes[2])
plt.tight_layout()

sq.im.segment(img=crop, layer="image_smooth", method="watershed", thresh=90, geq=False, channel=0)

print(crop)
print(f"Number of segments in crop: {len(np.unique(crop['segmented_watershed']))}")
fig, axes = plt.subplots(1, 2)
crop.show("image", channel=0, ax=axes[0])
_ = axes[0].set_title("H&E")
crop.show("segmented_watershed", cmap="jet", interpolation="none", ax=axes[1])
_ = axes[1].set_title("segmentation")


sq.im.calculate_image_features(adata_crop,crop,
    layer="image",features="segmentation",
    key_added="segmentation_features",
    features_kwargs={
        "segmentation": {
            "label_layer": "segmented_watershed",
            "props": ["label", "area", "centroid"],
            "channels": [1, 2],}
    }, mask_circle=True)

adata_crop.obsm["segmentation_features"].head()

adata_sml = adata[:50].copy()

# calculate default features
sq.im.calculate_image_features(
    adata_sml, img, features="summary", key_added="features", show_progress_bar=False)
# calculate features with masking
sq.im.calculate_image_features(adata_sml,img,
    features="summary",key_added="features_masked",
    mask_circle=True,show_progress_bar=False,)
# calculate features with scaling and larger context
sq.im.calculate_image_features(adata_sml,img,
    features="summary",key_added="features_scaled",
    mask_circle=True,spot_scale=2,
    scale=0.5,show_progress_bar=False,)

# plot distribution of median for different cropping options
sns.displot(
    {
        "features": adata_sml.obsm["features"]["summary_ch-0_quantile-0.5"],
        "features_masked": adata_sml.obsm["features_masked"]["summary_ch-0_quantile-0.5"],
        "features_scaled": adata_sml.obsm["features_scaled"]["summary_ch-0_quantile-0.5"],
    },
    kind="kde",)


```



```bash

## 可以启用一个新的环境
conda install -c conda-forge napari pyqt
pip install deprecated


import squidpy as sq
print(f"squidpy=={sq.__version__}")
import os
import scanpy as sc

adata = sq.datasets.visium_hne_adata()
# adata.write('visium_hne_adata.h5ad')
adata=sc.read('visium_hne_adata.h5ad')
img = sq.datasets.visium_hne_image()


viewer = img.interactive(adata)
##运行 viewer



```



```bash


from numpy.random import default_rng
import matplotlib.pyplot as plt
import scanpy as sc
import squidpy as sq
from anndata import AnnData
sc.logging.print_header()
print(f"squidpy=={sq.__version__}")

rng = default_rng(42)
counts = rng.integers(0, 15, size=(10, 100))  # feature matrix
coordinates = rng.uniform(0, 10, size=(10, 2))  # spatial coordinates
image = rng.uniform(0, 1, size=(10, 10, 3))  # image


adata = AnnData(counts, obsm={"spatial": coordinates})

sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata)
adata

sq.pl.spatial_scatter(adata, shape=None, color="leiden", size=50)



sq.gr.spatial_neighbors(adata, radius=3.0)
sq.pl.spatial_scatter(adata,color="leiden",
    connectivity_key="spatial_connectivities",
    edges_color="black",shape=None,
    edges_width=1,size=3000,)

plt.imshow(image)


spatial_key = "spatial"
library_id = "tissue42"
adata.uns[spatial_key] = {library_id: {}}
adata.uns[spatial_key][library_id]["images"] = {}
adata.uns[spatial_key][library_id]["images"] = {"hires": image}
adata.uns[spatial_key][library_id]["scalefactors"] = {
    "tissue_hires_scalef": 1,
    "spot_diameter_fullres": 0.5,}


sq.pl.spatial_scatter(adata, color="leiden")

adata.uns[spatial_key][library_id]["scalefactors"] = {
    "tissue_hires_scalef": 0.5,
    "spot_diameter_fullres": 0.5,}
sq.pl.spatial_scatter(adata, color="leiden", size=2)


img = sq.im.ImageContainer(image)
img.show()
```



```bash

import numpy as np
import squidpy as sq
arr = np.ones((100, 100, 3))
arr[40:60, 40:60] = [0, 0.7, 1]

print(arr.shape)
img = sq.im.ImageContainer(arr, layer="img1")
img

arr_seg = np.zeros((100, 100))
arr_seg[40:60, 40:60] = 1
img.add_img(arr_seg, layer="seg1")
img
img["seg2"] = arr_seg
img
print(list(img))
img["img1"]
img.rename("seg2", "new-name")


img.show(layer="img1")



crop1 = img.crop_corner(30, 40, size=(30, 30), scale=1)
crop1.show(layer="img1")

crop2 = crop1.crop_corner(0, 0, size=(40, 40), scale=0.5)
crop2.show(layer="img1")


print(crop1.data.attrs)
print(crop2.data.attrs)

sq.im.ImageContainer.uncrop([crop1], shape=img.shape).show(layer="img1")
sq.im.ImageContainer.uncrop([crop2], shape=(50, 50)).show(layer="img1")

img.data

img_on_disk = sq.datasets.visium_hne_image()
print(type(img_on_disk["image"].data))

img_on_disk.compute()
print(type(img_on_disk["image"].data))
```



```bash

import squidpy as sq
import scanpy as sc
import anndata as ad
import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt

## 更改完再检查一下有没有更改成功
print(os.getcwd())

#adata = sq.datasets.mibitof()
#adata.write("mibitof.h5ad")
adata=sc.read('mibitof.h5ad')
adata
df=adata.obs

sq.pl.spatial_scatter(adata,library_key= 'point16' ,color="Cluster" )
sq.pl.spatial_scatter(adata,library_key= 'point23' ,color="Cluster" )


sq.pl.spatial_scatter(adata,  color="Cluster" )

sc.pl.spatial(adata[adata.obs["library_id"] == 'point16'], color="Cluster",
        library_id='point16', title='point16', frameon=False,show=False)

fig, ax = plt.subplots(ncols=3, nrows=1, figsize=(14, 3))
i = 0
for library_id in adata.uns["spatial"].keys():
    sc.pl.spatial(
        adata[adata.obs["library_id"] == library_id], color="Cluster",
        library_id=library_id, title=library_id, frameon=False, ax=ax[i], show=False
    )
    i += 1


adata.uns.keys()
adata.uns['spatial'].keys()
adata.uns['spatial']['point23'].keys()
adata.uns['spatial']['point23']['images'].keys()


adata.var_names
adata.uns.keys()
adata.uns['spatial']['point8']['images'].keys()

df1=adata.obs
adata.obs.library_id
fig, ax = plt.subplots(ncols=3, nrows=1, figsize=(14, 3))
i = 0
for library_id in adata.uns["spatial"].keys():
    sc.pl.spatial(
        adata[adata.obs["library_id"] == library_id], color="Cluster",
        library_id=library_id, title=library_id, frameon=False, ax=ax[i], show=False
    )
    i += 1


plt.figure(figsize=(8, 8))
sc.pl.spatial(
        adata[adata.obs["library_id"] == 'point8'], library_id='point8', frameon=False, cmap='Spectral_r', size=0, ax=plt.gca(), img_key='segmentation'
    )




img =sq.im.ImageContainer(adata.uns["spatial"]['point16']["images"]["hires"], 
                          library_id='point16', layer='image')


img.add_img(adata.uns["spatial"]['point16']["images"]["segmentation"],
            library_id='point16', layer="segmentation")
img["segmentation"].attrs["segmentation"] = True

img.show(layer='image')
img.show(layer='segmentation')

img.show("image", segmentation_layer="segmentation", segmentation_alpha=0.5)
img.show("image", segmentation_layer="segmentation", segmentation_alpha=0.8)

for i, cmap in zip(range(3), ['Reds', 'Greens', 'Blues']):
    img.show(layer='image', channel=i, cmap=cmap)


img
img.crop_center(y=0.5, x=0.5, radius=200).show(layer='image')
img.crop_center(y=0.5, x=0.5, radius=100).show(layer='image')

img.crop_corner(y=0.2, x=0.3, size=400).show(layer='image')
img.crop_corner(y=0.2, x=0.3, size=200).show(layer='image')

subset_img= img.crop_corner(y=0.2, x=0.3, size=200)
img, subset_img

adata_crop = subset_img.subset(adata)
adata_crop
adata


sc.pl.spatial(
        adata_crop[adata_crop.obs["library_id"] == 'point16'], color=['CD45', 'CD3'],
        library_id='point16', frameon=False, cmap='Spectral_r', size=1.3)

```



```bash
06 部分代码
import squidpy as sq
print(f"squidpy=={sq.__version__}")
import os
import scanpy as sc
adata = sq.datasets.visium_hne_adata()
# adata.write('visium_hne_adata.h5ad')
adata=sc.read('visium_hne_adata.h5ad')
img = sq.datasets.visium_hne_image()
viewer = img.interactive(adata)
df=viewer.adata.obs
viewer.screenshot(canvas_only=False)
viewer.screenshot(canvas_only=False)
```





