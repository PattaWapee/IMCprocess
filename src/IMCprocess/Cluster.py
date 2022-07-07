import numpy as np
import scanpy as sc
import pandas as pd
from matplotlib import pyplot as plt
import phenograph
import copy
from scipy import sparse
from sklearn.metrics import adjusted_rand_score
import seaborn as sns


def plot_cumulative_var(adata, cml_var_explained):
    """
    The function to plot the cumulative explained variance ratio 
    as a function of the number of components


    :function: plot the cumulative explained variance ratio as a function of the number of components
    :returns: TODO

    """
    x = range(len(adata.uns['pca']['variance_ratio']))
    y = cml_var_explained
    fig = plt.figure()
    plt.scatter(x, y, s=4)
    plt.xlabel('Number of principle components')
    plt.ylabel('Cumulative explained variance')
    #plt.title('Cumulative explained variance')
    plt.show()
    return(fig)


def runPCA(adata):
    adata = sc.tl.pca(adata, n_comps=len(adata.var_names)-1, copy=True)
    cml_var_explained = np.cumsum(adata.uns['pca']['variance_ratio'])
    plot_cumulative_var(adata, cml_var_explained)

    # set the minimum cumulative fraction of variance explained
    min_cml_frac = 0.8
    # find the number of PCs that together explain the minimal cum. fraction of variance
    n_pcs = next(idx for idx, cml_frac in enumerate(
        cml_var_explained) if cml_frac > min_cml_frac)
    return adata, n_pcs


def run_phenograph(adata):
    k = 30  # choose k
    sc.settings.verbose = 0
    communities, graph, Q = phenograph.cluster(
        pd.DataFrame(adata.obsm['X_pca']), k=k)  # run PhenoGraph
    # store the results in adata_level1:
    adata.obs['PhenoGraph_clusters'] = pd.Categorical(communities)
    adata.uns['PhenoGraph_Q'] = Q
    adata.uns['PhenoGraph_k'] = k

    return adata
