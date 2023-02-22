import os, sys
sys.path.append(os.path.join('../IMCprocess/'))

from IMCprocess.clustifier import clustifier as cl
import unittest

# PYTHONPATH=src python -m pytest tests/test_cluster.py
import anndata as ad
import scanpy as sc
import numpy as np
from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt



class TestClustifier(unittest.TestCase):
    def create_random_anndata(self):
        # Create anndata object for test
        counts = csr_matrix(np.random.poisson(1, size=(1000, 20)))
        adata = ad.AnnData(counts)
        adata.obs_names = [f"Cell_{i+1}" for i in range(adata.n_obs)]
        adata.var_names = [f"marker_{i+1}" for i in range(adata.n_vars)]
        adata.obsm["spatial"] = np.array(list(zip(np.random.rand(1000), np.random.rand(1000))))
        return(adata)

    def test_plot_cumulative_var(self):
        # Check if it returns a figure object
        adata = self.create_random_anndata()
        adata = sc.tl.pca(adata, n_comps=len(adata.var_names)-1, copy=True)
        cml_var_explained = np.cumsum(adata.uns['pca']['variance_ratio'])
        fig = cl.plot_cumulative_var(adata, cml_var_explained)
        self.assertIsInstance(fig, type(plt.figure()))

if __name__ == '__main__':
    unittest.main()
