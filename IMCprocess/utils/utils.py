import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad
import pickle
from upsetplot import plot
from upsetplot import UpSet
import json
from skimage import io
from skimage import measure


def transfer_obs(adata_with_info, adata_get_info, info_name):
    """
    transfer obs level info from adata_with_info to adata_get_info
    by matching index
    """
    info = adata_with_info.obs[info_name]
    adata_get_info.obs[info_name] = info[info.index.isin(
        adata_get_info.obs_names)]
    return adata_get_info


def merge_adata_obs(adata1, adata2, obs_ls1, obs_ls2):
    '''
    merge two anndata and obs by create new anndata with merge obs columns
    Parameters:
    ___________
    adata1, adata2 = anndata objects for merging
    obs_ls1 = list of obs variable from adata1
    obs_ls2 = list of obs variable from adata2
    return:
    _______
    merge_adata with new obs columns as obs1+'_'+obs2
    '''
    mergeobs = pd.concat([pd.concat([adata1.obs[obs1].rename(obs1+'_'+obs2),
                                     adata2.obs[obs2].rename(obs1+'_'+obs2)])
                          if obs1 != obs2 else
                          pd.concat([adata1.obs[obs1], adata2.obs[obs2]])
                          for obs1, obs2 in zip(obs_ls1, obs_ls2)], axis=1)
    # change all columns to categories
    mergeobs[mergeobs.columns] = mergeobs[mergeobs.columns].astype('category')

    new_adata1 = ad.AnnData(adata1.to_df())
    new_adata2 = ad.AnnData(adata2.to_df())

    merge_adata = ad.concat([new_adata1, new_adata2])
    merge_adata.obs = mergeobs
    return merge_adata

def add_level2_to_level1_obs(level1_adata, level2_adata_list):
    '''
    Add level2 annotation to level1 adata object
    '''
    level1_adata.obs['level2'] = level1_adata.obs['level1_annotated'].astype(str)
    for lev2_adata in level2_adata_list:
        lev2_col = lev2_adata.obs.columns[lev2_adata.obs.columns.str.contains('level2')][0]
        lev2_adata.obs['level2'] = lev2_adata.obs[
            'level1_annotated'].astype('str') + '_' + lev2_adata.obs[lev2_col].astype('str')
        level1_adata.obs.loc[lev2_adata.obs.level2.index, 'level2'] = lev2_adata.obs.level2
    return level1_adata

def get_map_dict(annot_cluster_dict):
    map_dict = {}
    for key, value in annot_cluster_dict.items():
        for i in value:
            map_dict[i] = key
    return(map_dict)

def read_pickle_obj(pickle_file):
    open_file = open(pickle_file, "rb")
    data = pickle.load(open_file)
    open_file.close()
    return(data)


def save_pickle_obj(file_name, obj_file):
    open_file = open(file_name, "wb")
    pickle.dump(obj_file, open_file)
    open_file.close()


def filter_adata_obs(adata, obs_var, obs_value_ls):
    """
    filter anndata by list of obs value
    Parameters:
    ___________
    adata = anndata object
    obs_var = obs variable for filtering
    obs_value_ls = list of obs variable categories for filtering
    return:
    _______
    filtered adata

    """
    adata_filtered = adata[adata.obs[obs_var].isin(obs_value_ls)]
    return adata_filtered



######################
### PLOT FUNCTIONS ###
######################

def heatmap_dendrogram(adata, groupby, figsize, markers = 'all'):
    """
    plot heatmap with dendrogram from anndata
    Parameters:
    ___________
    adata = anndata object
    groupby = The obs group label for plotting in heatmap 
            ex. 'PhenoGraph_clusters' (colname from anndata.obs)
    figsize = tuple of figure size
    markers = marker list to plot in heatmap (default = 'all')
    return:
    _______
    list of Axes
    Note:
    _______
    To save plot from this object.
    >> plot['heatmap_ax'].figure.savefig('output_plot/Fib_level2_AB.png', dpi = 300)

    """
    if markers != 'all':
        adata = adata[:, markers]
    fig = sc.pl.heatmap(adata,
                        adata.var_names,
                        groupby=groupby,
                        swap_axes=True,
                        cmap='viridis',
                        dendrogram=True,
                        show=False,
                        figsize=figsize)
    return fig

def plt_umap(adata, groupby):
    fig = sc.pl.umap(adata,
                     color=[groupby],
                     add_outline=True,
                     legend_loc='on data',
                     legend_fontsize=8,
                     palette='tab20',
                     frameon=True,
                     show=False)
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

#################
### CELL MASK ###
#################


def label_class(json_all, label_index, name_class_cell, color, mask):
    """
    function to create geoJson object for cell mask with color label

    Parameters:
    ___________
    json_all = json object for append new cell mask
    label_index = cell index to create json object with color label
    name_class_cell = name of the cell 
    color = color to label cell
    mask = mask file
    return:
    _______
    json object

    """
    for label_i in label_index:
        contour = measure.find_contours(mask == label_i, 0.5)
        evalCloseness = contour[0][-1] == contour[0][0]
        isClosed = False not in evalCloseness
        if isClosed:
            y, x = contour[0].T
            coordinates_inv = np.dstack((x, y))
            coordinates = coordinates_inv[0].tolist()
            json_object = {'type': 'Feature',  # 'id': 'PathAnnotationObject',
                           'id': 'PathDetectionObject',
                           'geometry': {'type': 'Polygon', 'coordinates': [coordinates]}, 'properties': {"classification": {"name": name_class_cell, "colorRGB": color}, 'isLocked': False, 'measurements': []}}
            json_all.append(json_object)
        else:
            coords = contour[0]
            coords_closed = np.append(coords, [coords[0]], axis=0)
            contour_closed = [coords_closed]
            y, x = contour_closed[0].T
            coordinates_inv = np.dstack((x, y))
            coordinates = coordinates_inv[0].tolist()
            json_object = {'type': 'Feature',
                           # 'id': 'PathAnnotationObject',
                           'id': 'PathDetectionObject',
                           'geometry': {'type': 'Polygon', 'coordinates': [coordinates]}, 'properties': {"classification": {"name": name_class_cell, "colorRGB": color}, 'isLocked': False, 'measurements': []}}
            json_all.append(json_object)

    return(json_all)


def export_labelmask_geoJson(mask_file, output_file, label_dict, color_label):
    """
    function to create geoJson cell mask with multiple cell label
    Parameters:
    ___________
        mask = probability mask file from cellprofiler
        output_file = geoJson output path file
        label_dict = dictionary of label cell
        color_label = list of color for cell
            (if there are only two type of cells: positive & negative,
            we use color [-1,-8245601 (purple)]

            )
    return:
    _______
    json object
    """

    # read the image as np-array
    mask = io.imread(mask_file)
    # regionprops_table below requires inte
    mask = mask.astype(np.uint32)
    # measure properties of all cells
    results = measure.regionprops_table(mask, properties=['label', 'coords'])
    # convert results into pandas dataframe
    results_df = pd.DataFrame.from_dict(results)

    # label_class
    name_class_list = list(label_dict.keys())
    #color_label = [-1, -8245601]
    json_export = []
    for i in range(len(name_class_list)):
        json_export = label_class(
            json_export, label_dict[name_class_list[i]], name_class_list[i], color_label[i], mask)

    with open(output_file, 'w') as outfile:
        json.dump(json_export, outfile)

    return json_export
