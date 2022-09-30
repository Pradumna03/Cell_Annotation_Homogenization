# Task1_Ten_60_BIO
Scanpy is a python library which is used to analyze single cell Gene Expression Data.

Anndata is a type of matrix data structure organization class that has a specific form of data matrix for a specific type of data 
![Anndata_str](https://user-images.githubusercontent.com/99180702/193336621-2d53ca6c-e8bc-4682-a0a8-df23859d6621.png)

```
import scvi
import scanpy as sc
import numpy as np
import sys
import matplotlib
import leidenalg
```

Importing the libraries needed

```
adata = sc.read_h5ad('/content/drive/MyDrive/train.h5ad')
```

Reading the data from a mounted google drive in Google Colab.

```
sc.pp.filter_genes(adata, min_counts=3)
```
Filters genes that are expressed in less than 3 cells.

```
#sc.pp.normalize_total(adata, target_sum=1e4)
#sc.pp.log1p(adata)
```
These steps are skipped as the data provided already contains the adata.layers["Counts"] layer which is created after doing these steps and can then be stored as a layer. A layer can be understood as freezed state of observation vs. variable expression matrix X. 

```
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=1200,
    subset=True,
    layer="counts",
    flavor="seurat_v3",
    batch_key="patient"
)
```
This step serves to identify the highly variable genes and the various attributes are:

n_top_genes = To identify the top n variable genes.

subset = True means that adata would be converted to a subset of highly variable genes. If False it will not convert

layer = Layer of Data to be used for analysis

flavor = seurat_v3 used

batch_key = patient is used as a batch key here to identify the batches of data that have been extracted from different sources(batches) and then put together for analysis.

We now create PC's or Principal Components which serve to reduce dimensionality of our matrix by serving as efficient markers of distinction. These PC's can be created by:

```
sc.pp.pca(adata)
```
All information about this function can be accessed through:
```
?sc.pl.pca_overview
```
We can also plot the variance ratio by:
```
sc.pl.pca_variance_ratio(adata, log=True)
```
The plot looks like:![variance_ratio](https://user-images.githubusercontent.com/99180702/193353671-e08bef1d-6acd-485d-8b2e-ecb48a1160ed.png)

If the variance ratio difference will be large for PC's we can simply use normalization by CPM approach to remove the highly expressed PC's and make the distribution better and more even for the rest of the PC's. Here that is not needed.

The steps to creating a set of marker genes involves computing neighbourhood graphs and then embedding and clustering the neighbourhood graph. I do not have an absolutely clear idea about the maths involved in these steps and the working basis of these steps.
```
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)
sc.tl.leiden(adata)
sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
```
The ranked genes can now be used as marker genes for sell type prediction.

A model can also be trained for cell type annotation using scVI.

```
model = scvi.model.SCVI(adata)
model.train()
```
