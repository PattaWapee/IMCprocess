import os, sys
sys.path.append(os.path.join('../IMCprocess/'))

from IMCprocess.clustifier import *

#from IMCprocess.clustifier import *
# PYTHONPATH=src python -m pytest tests/test_cluster.py
import anndata as ad
import numpy as np
from scipy.sparse import csr_matrix

# Create anndata object for test
counts = csr_matrix(np.random.poisson(1, size=(1000, 20)))
adata = ad.AnnData(counts)
adata.obs_names = [f"Cell_{i+1}" for i in range(adata.n_obs)]
adata.var_names = [f"marker_{i+1}" for i in range(adata.n_vars)]
adata.obsm["spatial"] = np.array(list(zip(np.random.rand(1000), np.random.rand(1000))))

print(adata.to_df().head())