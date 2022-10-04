# Task1_Ten_60_BIO
Scanpy is a python library which is used to analyze single cell Gene Expression Data.

Anndata is a type of matrix data structure organization class that has a specific form of data matrix for a specific type of data 
![Anndata_str](https://user-images.githubusercontent.com/99180702/193336621-2d53ca6c-e8bc-4682-a0a8-df23859d6621.png)

```
import scvi
import scanpy as sc
import numpy as np
import sys
```

Importing the libraries needed

```
adata = sc.read_h5ad('File Path') # Training Dataset
bdata = sc.read_h5ad('File Path') # Testing Dataset
```

```
sc.pp.filter_genes(adata, min_counts=3)
sc.pp.filter_genes(bdata, min_counts=3)
```
Filters genes that are expressed in less than 3 cells.

```
#sc.pp.normalize_total(adata, target_sum=1e4)
#sc.pp.log1p(adata)
```
These steps are skipped as the data provided already contains the adata.layers["Counts"] layer which is created after doing these steps and can then be stored as a layer. A layer can be understood as frozen state of observation vs. variable expression matrix X after some mathematical computation on the dataframe. 

```
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=1200,
    subset=True,
    layer="counts",
    flavor="seurat_v3",
    batch_key="patient"
)
sc.pp.highly_variable_genes(
    bdata,
    n_top_genes=1200,
    subset=True,
    layer="counts",
    flavor="seurat_v3",
    batch_key="patient"
)
```

```
scvi.model.SCANVI.setup_anndata(adata, batch_key="batch", labels_key="labels")
model = scvi.model.SCANVI(adata, "Unknown")
model.train() # model being trained on adata 
model.save("my_model/") # To save the model
latent = model.get_latent_representation(bdata) # Dimensionality Reduction and storing of values as a matrix in .obsm in bdata
bdata.obsm["X_scVI"] = latent
bdata.obs["pred_label"] = model.predict(bdata) # model prediction for bdata test data set
col1=bdata.obs["pred_label"] # predicted cell type values stored in bdata column pred_label
col2=bdata.obs["cell_types"] # accurate cell type values stored in bdata column cell_types
Arr1=np.array(col1) 
Arr2=np.array(col2)
Total_trues=np.sum(Arr1==Arr2) # Individual matching values between predicted and accurate cell type counted
efficiency=Total_trues/Arr1.size() # Efficiency calculated
print(efficiency)
```
