import os
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import gc
import time
import matplotlib.pyplot as plt
import scipy.sparse as sparse
from mygene import MyGeneInfo
from functions_python import compute_expression_margins, cell_specific_gene_expression, compute_positivity

#### Lung cancer post-EGFR-TKI
# Overall structure of .h5ad file
### Lung cancer post-EGFR-TKI single cells
adata_lung_primary_all=sc.read_h5ad('Lung_primary_posttki.h5ad')

#### if the count matrix applied log-normalization, transform to the original counts (only CP10K)
if sparse.issparse(adata_lung_primary_all.X):
    X = adata_lung_primary_all.X.copy()
    X.data = np.expm1(X.data)
    adata_lung_primary_all.X = X
else:
    adata_lung_primary_all.X = np.expm1(adata_lung_primary_all.X)

adata_lung_primary_tumor = adata_lung_primary_all[(adata_lung_primary_all.obs['cell_type_reannotation'].isin(['Cancer cell']))].copy()
adata_lung_normal =  adata_lung_primary_all[(adata_lung_primary_all.obs['cell_type_reannotation'].isin(['Epithelial cell']))].copy()
gc.collect()
