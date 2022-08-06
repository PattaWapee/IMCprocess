import scanpy as sc
import pickle
from upsetplot import plot
from upsetplot import UpSet


def read_pickle_obj(pickle_file):
    open_file = open(pickle_file, "rb")
    data = pickle.load(open_file)
    open_file.close()
    return(data)


def save_pickle_obj(file_name, obj_file):
    open_file = open(file_name, "wb")
    pickle.dump(obj_file, open_file)
    open_file.close()


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


######################
### PLOT FUNCTIONS ###
######################


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


def plt_upset(binary_df):
    df = binary_df.astype(bool)
    df['index'] = df.index
    df = df.set_index(binary_df.columns.to_list())
    upset_res = UpSet(df, subset_size='count',
                      show_counts=True, sort_categories_by=None)
    plot = UpSet(df, subset_size='count',
                 show_counts=True, sort_categories_by=None).plot()

    return plot
