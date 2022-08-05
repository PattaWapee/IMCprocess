import scanpy as sc


def filter_adata(adata, obs_var, obs_value_ls):
    """
    filter anndata by list of obs value
    Parameters:
    ___________
    adata = anndata object
    k = number of nearest neighbors to use in the first step of graph constructions
        default setting is 30
    return:
    _______
    adata with PhenoGraph_clusters data


    """
    pass


def heatmap_dendrogram(adata, groupby, figsize):
    """
    plot heatmap with dendrogram from anndata
    Parameters:
    ___________
    adata = anndata object
    groupby = The obs group label for plotting in heatmap ex. 'PhenoGraph_clusters' (colname from anndata.obs)
    return:
    _______
    list of Axes
    Note:
    _______
    To save plot from this object.
    >> plot['heatmap_ax'].figure.savefig('output_plot/Fib_level2_AB.png', dpi = 300)

    """

    fig = sc.pl.heatmap(adata,
                        adata.var_names,
                        groupby=groupby,
                        swap_axes=True,
                        cmap='viridis',
                        dendrogram=True,
                        show=False,
                        figsize=figsize)
    return fig
