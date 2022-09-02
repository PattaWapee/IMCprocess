import numpy as np
import scanpy as sc
import pandas as pd
from matplotlib import pyplot as plt


def transfer_spatial(adata_with_xy, adata_get_xy):
    """
    transfer spatial data by matching index cells 
    adata_with_xy = anndata with spatial data
    adata_get_xy = anndata to be transfered spatial data from adata_with_xy
    """

    adata_filter = adata_with_xy[adata_get_xy.obs_names]
    adata_get_xy.obsm['spatial'] = adata_filter.obsm['spatial']

    return adata_get_xy


def plt_intensity(adata, marker, path, name):
    fig = adata.to_df()[marker].plot.density()
    plt.title(name + marker)
    if path != None:
        plt.savefig(path+'/' + marker + '_intenisty.png', dpi=300)
    plt.show()
    return fig
