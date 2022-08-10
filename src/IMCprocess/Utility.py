import scanpy as sc
import numpy as np
import pandas as pd
import pickle
from upsetplot import plot
from upsetplot import UpSet
import json
from skimage import io
from skimage import measure


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
