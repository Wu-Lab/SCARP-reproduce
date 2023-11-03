import pandas as pd
import scanpy as sc
import numpy as np
import re
import matplotlib.pyplot as plt
import networkx as nx
from scipy.sparse import csr_matrix
import seaborn as sns
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import precision_recall_curve


# [a,b], [c,d]
def is_pchic_validated(a, b, c, d, thre1=0.5, thre2=1000):
    if a > c:
        peak1 = [c, d]
        peak2 = [a, b]
    else:
        peak1 = [a, b]
        peak2 = [c, d]

    if peak2[0] > peak1[1]:
        return 0
    else:
        judge1 = max((peak1[1] - peak2[0]) / (b - a), (peak1[1] - peak2[0]) / (d - c))
        judge2 = abs((a+b)/2-(c+d)/2)
        if (judge1 > thre1) or (judge2 < thre2):
            return 1
        else:
            return 0


def unpivot(frame):
    N, K = frame.shape
    data_unpivot = {
        "Regulatory score": frame.to_numpy().ravel("F"),
        "Peaks": np.asarray(frame.columns).repeat(N),
        "Promoters": np.tile(np.asarray(frame.index), K),
    }
    return pd.DataFrame(data_unpivot, columns=["Promoters", "Peaks", "Regulatory score"])


def filter_cells(adata, keep_cells):
    adata_df = pd.DataFrame(adata.X.todense(),
                            index=adata.obs.index,
                            columns=adata.var.index)
    adata_new = sc.AnnData(adata_df.loc[keep_cells])
    adata_new.obs = adata.obs.loc[keep_cells]
    adata_new.var = adata.var
    adata_new.X = csr_matrix(adata_new.X)
    return adata_new


def filter_peaks(adata, keep_peaks):
    adata_df = pd.DataFrame(adata.X.todense(),
                            index=adata.obs.index,
                            columns=adata.var.index)
    adata_new = sc.AnnData(adata_df[keep_peaks])
    adata_new.obs = adata.obs
    adata_new.var = adata.var.loc[keep_peaks]
    adata_new.X = csr_matrix(adata_new.X)
    return adata_new


def getNClusters_Louvain(adata, n_cluster, range_min=0, range_max=3, max_steps=50, random_seed=1):
    np.random.seed(random_seed)
    temp_step = 0
    temp_min = float(range_min)
    temp_max = float(range_max)
    while temp_step < max_steps:
        temp_resolution = temp_min + ((temp_max - temp_min) / 2)
        sc.tl.louvain(adata, resolution=temp_resolution)
        temp_clusters = adata.obs['louvain'].nunique()

        if temp_clusters > n_cluster:
            temp_max = temp_resolution
        elif temp_clusters < n_cluster:
            temp_min = temp_resolution
        else:
            return [1, adata]
        temp_step += 1

    print('Cannot find the number of clusters')
    print('Clustering solution from last iteration is used:' +
          str(temp_clusters) + ' at resolution ' + str(temp_resolution))
    sc.tl.louvain(adata, resolution=temp_resolution)
    adata.obs['louvain'].nunique()
    return [0, adata]