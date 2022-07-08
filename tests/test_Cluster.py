from IMCprocess.cluster import *
# PYTHONPATH=src python -m pytest tests/test_cluster.py
import anndata as ad
import numpy as np
from scipy.sparse import csr_matrix

# Create anndata object for test
counts = csr_matrix(np.random.poisson(1, size=(1000, 20)))
adata = ad.AnnData(counts)
adata.obs_names = [f"Cell_{i}" for i in range(adata.n_obs)]
adata.var_names = [f"marker_{i}" for i in range(adata.n_vars)]
