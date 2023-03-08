import numpy as np
import scanpy as sc
import pandas as pd
from matplotlib import pyplot as plt
import os, sys
import anndata as ad

import itertools
import random
import string
other_dir_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.append(other_dir_path)
from IMCprocess.img_anndata.img_anndata import Img_anndata
import IMCprocess.spatial_analysis.spatial_analysis as sa

######################################################
# Function to find enrichment of markers combination #
######################################################


## Create a random img data for testing
def rand_marker_XY(ncol = 9, nrow=1000):
    columns_list=[string.ascii_uppercase[i] for i in range(ncol)]
    df = pd.DataFrame(np.random.rand(nrow, ncol), columns=columns_list)
    df['X_position'] = [random.uniform(0, 256) for _ in range(nrow)]
    df['Y_position'] = [random.uniform(0, 256) for _ in range(nrow)]
    return(df)

def binary_dataframe(df, threshold=0.5):
    '''
    Return a DataFrame with binary columns based on a threshold
    '''
    binary_df = pd.DataFrame()
    for col in df.columns:
        binary_col = np.where(df[col] > threshold, 1, 0)
        binary_df[col] = binary_col
    return binary_df

def label_marker(binary_df, marker_list):
    df = binary_df[marker_list]
    for col in df.columns:
        df.loc[df[col] == 1, col] = col+'+'
        df.loc[df[col] == 0, col] = col+'-'
    df_label = pd.DataFrame(df.apply(lambda row: ''.join(row.astype(str)), axis=1), columns=['label'])
    df_label['label'] = df_label['label'].astype('category')
    return df_label

def add_obs_label(adata, marker_list):
    binary_df = binary_dataframe(adata.to_df())
    label_df = label_marker(binary_df, marker_list)
    label_df.index = adata.obs.index
    adata.obs = pd.concat([adata.obs, label_df], axis=1)
    return adata
def get_zscore(merged_adata):
    '''
    Calculate Neighborhood enrichment zscore for each label
        Parameters
        ----------
        merged_adata: anndata object of two cell types
        Returns
        -------
        zscore: zscore dataframe of each label
    '''
    sa.run_spatial_nhood(merged_adata, 'label', 15)
    zscore = merged_adata.uns['label_nhood_enrichment']['zscore'][-1]
    label_index = merged_adata.obs['label'].cat.categories
    zscore = pd.DataFrame(zscore[:-1], index=label_index[:-1], columns=['ZZcell1']).T
    return zscore

def label_merge_adata(cell1_adata, cell2_adata, marker_list):
    # 1. label cell1 with zzcell1 and cell2 with combination markers + and -
    cell2_adata = add_obs_label(cell2_adata,marker_list)
    cell1_adata.obs['label'] = 'ZZcell1'
    # 2. Merge two cell types
    merged_adata = ad.concat([cell1_adata, cell2_adata], axis=0)
    merged_adata.obs['label'] = merged_adata.obs['label'].astype('category')
    return merged_adata

def col_combinations(adata, n=2):
    '''
    Return all combinations of columns in a dataframe
    '''
    combinations = list(itertools.combinations(adata.to_df().columns, n))
    combinations_as_lists = [list(c) for c in combinations]
    return(combinations_as_lists)

def zscore_combinations(cell1_adata, cell2_adata, n=2):
    '''
    Get zscore list of all combinations of columns
    '''
    zscore_ls = []
    for set in col_combinations(cell2_adata, n = 2):
        print('start testing combination: ')
        print(set)
        adata1 = cell1_adata.copy()
        adata2 = cell2_adata.copy()
        merged_adata = label_merge_adata(adata1, adata2, set)
        # 3. calculate zscore of neighborhood enrichment
        zscore = get_zscore(merged_adata)
        zscore_ls.append(zscore)
    return zscore_ls

## main function to test a function of enrichment of markers combination
def main():
    print('Start testing enrichment of markers combination')
    ## 1. Create random anndata for Two cell types
    df1 = rand_marker_XY()
    df2 = rand_marker_XY()
    df2.index = list(range(1001, 2001))
    cell1_adata = Img_anndata(dfs=[df1], img_ids=['R1']).adata
    cell2_adata = Img_anndata(dfs=[df2], img_ids=['R1']).adata

    zscore_ls = []
    # 2. label cell1 with zzcell1 and cell2 with combination markers + and -
    #  and merge two cell types and calculate zscore of neighborhood enrichment

    zscore_ls = zscore_combinations(cell1_adata, cell2_adata, n=2)
    '''
    for set in col_combinations(cell2_adata, n = 2):
        print(set)
        adata1 = cell1_adata.copy()
        adata2 = cell2_adata.copy()
        merged_adata = label_merge_adata(adata1, adata2, set)
        # 3. calculate zscore of neighborhood enrichment
        zscore = get_zscore(merged_adata)
        zscore_ls.append(zscore)

    '''
    print(len(zscore_ls))
                             


if __name__ == '__main__':
    main()


