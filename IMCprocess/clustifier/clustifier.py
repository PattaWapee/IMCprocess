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
    :returns: figure of cumulative explain variance plot

    """
    x = range(len(adata.uns['pca']['variance_ratio']))
    y = cml_var_explained
    fig = plt.figure()
    plt.scatter(x, y, s=4)
    plt.xlabel('Number of principle components')
    plt.ylabel('Cumulative explained variance')
    # plt.title('Cumulative explained variance')
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


def runPhenograph(adata, k=30):
    """
    Parameters:
    ___________
    adata = anndata object
    k = number of nearest neighbors to use in the first step of graph constructions
        default setting is 30
    return:
    _______
    adata with PhenoGraph_clusters data


    """

    sc.settings.verbose = 0
    communities, graph, Q = phenograph.cluster(
        pd.DataFrame(adata.obsm['X_pca']), k=k)  # run PhenoGraph
    # store the results in adata_level1:
    adata.obs['PhenoGraph_clusters'] = pd.Categorical(communities)
    adata.uns['PhenoGraph_Q'] = Q
    adata.uns['PhenoGraph_k'] = k

    return adata


def copy_property(adata1, adata2):
    """

    To transfer property from adata1 to adata2
    adata1 is clustered anndata with filtered markers
    adata2 is anndata with complete markers
    Parameters:
    ___________
    adata1 = anndata object with property to transfer to adata2
    adata2 = anndata object to get property from adata1
    return:
    _______
    adata2 with PhenoGraph_clusters data

    """
    adata2.obs = adata1.obs.copy()
    adata2.uns = adata1.uns.copy()
    adata2.obsm = adata1.obsm.copy()
    adata2.obsp = adata1.obsp.copy()
    # transfer PCs to adata2
    # adata1 and adata2 has different number of vars so we save in uns
    adata2.uns['PCs'] = adata1.varm['PCs']
    return adata2


def main_clustering(adata, markers):
    """
    function to run phenograph clustering for specific markers selections
    Parameters:
    ___________
    adata: anndata object
    markers: list of markers for clustering
    return:
    _______
    adata: full markers adata with PhenoGraph_clusters property

    """
    # 1. Filter ann data with markers
    adata_markers = adata[:, markers]

    # 2. run PCA
    adata_markers, n_pcs = runPCA(adata_markers)

    # 3. run Phenograph
    adata_markers = runPhenograph(adata_markers, k=30)

    # 4.Computing the neighborhood graph
    sc.pp.neighbors(adata_markers, n_pcs=n_pcs)

    # 5. Embedding the neighborhood graph with UMAP
    sc.tl.umap(adata_markers)

    # 6. copy phenograph property from fileter markers anndata to full markers anndata
    copy_property(adata_markers, adata)

    # 7. add dendogram
    sc.tl.dendrogram(adata, groupby='PhenoGraph_clusters')

    return adata
