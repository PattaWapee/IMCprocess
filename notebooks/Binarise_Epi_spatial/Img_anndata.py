import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd
import sys, os
#import clustifier as cl
#import utils as ut

other_dir_path = os.path.abspath(os.path.join(os.getcwd(), '..', '..'))
sys.path.append(other_dir_path)
import IMCprocess.utils.utils as ut
from IMCprocess.clustifier import clustifier as cl




class Img_anndata:
    def __init__(self, dfs, img_ids):
        self.dfs = dfs
        self.img_ids = img_ids
        self.adata = self.create_anndata()
        self.adata_dict = {}

    def get_img_id(self):
        return(self.adata.obs['img_id'].unique())

    def create_img_adata(self, img_id):
        return(self.adata[self.adata.obs['img_id'] == img_id])

    def create_spatial(self, Xseries, Yseries):
        self.adata.obsm["spatial"] = np.array(list(zip(Xseries, Yseries)))
        return self.adata

    def create_anndata(self):
        """
        Description
        -----------
            Create anndata from pandas dataframe for single or multiple images
            with meta data of X,Y and img_id
        Parameters
        ----------

        Returns
        -------
            adata : anndata
        """
        adata_ls = []
        for df, img_id in zip(self.dfs, self.img_ids):
            if sum(df.columns.isin(['X_position', 'Y_position'])) == 2:
                df_sc = df.dropna(axis= 'columns').drop(
                    ['Image_id','Mask_id','X_position','Y_position'],
                    axis = 1, errors='ignore')
            ada = sc.AnnData(df_sc)
            ada.var_names = df_sc.columns
            ada.obs['img_id'] = [img_id] * len(df)
            if sum(df.columns.isin(['X_position', 'Y_position'])) == 2:
                ada.obsm["spatial"] = np.array(
                    list(zip(df['X_position'], df['Y_position'])))
            adata_ls.append(ada)
        merge_adata = sc.concat(adata_ls, index_unique="_")

        return(merge_adata)

    #######################
    ##### Clustering ######
    #######################

def cluster_phenograph(img_anndata, adata_input, markers, name_level):
    """
    Description
    -----------
        Cluster cells using phenograph algorithm
    Parameters
    ----------

    Returns
    -------
    """
    adata = adata_input.copy()
    cl.main_clustering(adata, markers)
    add_adata_dict(img_anndata, adata, name_level)
        

def add_adata_dict(img_anndata, adata_new, name_adata):
    img_anndata.adata_dict[name_adata] = adata_new

def annotate_cluster(img_anndata, name_level, annote_cluster_dict):
    map_cluster_dict = ut.get_map_dict(annote_cluster_dict)
    img_anndata.adata_dict[name_level].obs[name_level+'_annotated'] = (img_anndata.adata_dict[
        name_level].obs['PhenoGraph_clusters'].map(
            map_cluster_dict).astype('category'))
