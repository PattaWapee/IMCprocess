import numpy as np
import scanpy as sc
import pandas as pd
from matplotlib import pyplot as plt


def plt_intensity(adata, marker, path, name):
    fig = adata.to_df()[marker].plot.density()
    plt.title(name + marker)
    if path != None:
        plt.savefig(path+'/' + marker + '_intenisty.png', dpi=300)
    plt.show()
    return fig
