import scanpy as sc
import pandas as pd
import numpy as np
import os

import anndata
import anndata2ri
import logging

import rpy2.rinterface_lib.callbacks as rcb
import rpy2.robjects as ro

rcb.logger.setLevel(logging.ERROR)
ro.pandas2ri.activate()
anndata2ri.activate()
from rpy2.robjects.packages import importr
decontX = importr('decontX')

data_dir = ''
for f in os.listdir(data_dir):
    adata = sc.read_h5ad(os.path.join(data_dir, f))

    batch_adatas = []

    for batch in adata.obs['batch'].cat.categories:
        adata_temp = adata[adata.obs['batch'] == batch].copy()
        adata_temp.layers['counts'] = adata_temp.X
        adata_temp_adj = adata_temp.copy()
        adata_temp_adj = decontX.decontX(adata_temp_adj)
        adata_temp.layers['counts_adj'] = np.round(adata_temp_adj.layers['decontXcounts'])
        batch_adatas += [adata_temp]

    adata_adj = anndata.concat(batch_adatas)

    adata_adj.X = adata_adj.layers['counts_adj']

    adata_adj.write(f)
