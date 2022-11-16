import anndata as ad
import scanpy as sc
import numpy as np

def create_anndata(df, X, Y, img_id):
    """
        Description
        -----------
            Create anndata from pandas dataframe 
            with meta data of X,Y and img_id
        Parameters
        ----------
            df : pandas dataframe
            X: list or pandas series for x position
            Y: list or pandas series for y position
            img: image id
            
        Returns
        -------
            adata : anndata
    """
    
    df.index = [f"Cell_{i:d}" for i in range(1,len(df)+1)]
    adata = sc.AnnData(df)
    adata.obs_names = [f"Cell_{i:d}" for i in range(1,len(adata)+1)]
    adata.var_names = df.columns
    adata.obsm["spatial"] = np.array(list(zip(X,Y)))
    adata.obs['img_id'] = [img_id]* len(df)
    return(adata)