import numpy as np
import scanpy as sc
import squidpy as sq
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


def plt_spatial(adata, obs_col, output_path, name):
    """
    function to plot 
    """
    fig = sc.pl.spatial(adata, color=obs_col,
                        title=[name+' '+obs_col],
                        return_fig=True,
                        #save = output_path + '/' + marker+ '_' + name +'.png',
                        spot_size=10, show=False)
    fig2 = fig[0].get_figure()
    fig2.tight_layout()
    fig2.savefig(output_path+obs_col+'_'+name+'_spatial.png', dpi=300)


def run_spatial_nhood(adata, obs_col, radius):
     # add None to obs_col if there is nan in obs_col
    adata.obs[obs_col] = pd.Categorical(adata.obs[obs_col])
    if adata.obs[obs_col].isna().any():
        adata.obs[obs_col] = adata.obs[obs_col].cat.add_categories(['None'])
        adata.obs[obs_col].fillna('None', inplace=True)

    sq.gr.spatial_neighbors(adata, radius=radius,coord_type='generic')
    sq.gr.nhood_enrichment(adata, cluster_key=obs_col)
    # calculate pvalue from zscore and save to uns
    adata.uns[obs_col+'_nhood_enrichment']['pval'] = z_to_pvalue(
        adata.uns[obs_col+'_nhood_enrichment']['zscore'])

def plt_spatial_nhood(adata, obs_col, output_path, name, img_size=(8, 5)):
    plt.figure(figsize=img_size)
    sq.pl.nhood_enrichment(adata, cluster_key=obs_col,
                           annotate=True,
                           save=output_path + obs_col + '_' +
                           name + '_NH_enrich.png',
                           dpi=300
                           )


def plt_interaction_mat(adata, obs_col, output_path, name):

    sq.gr.interaction_matrix(adata, cluster_key=obs_col)
    sq.pl.interaction_matrix(adata, cluster_key=obs_col,
                             annotate=True,
                             save=output_path + obs_col + '_' +
                             name + '_interaction_mat.png',
                             dpi=300)


import scipy.stats
def z_to_pvalue(z_array):
    df = pd.DataFrame(z_array)
    pval = scipy.stats.norm.sf(abs(df))
    return(pval)