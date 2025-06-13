import itertools
import networkx as nx
from networkx.algorithms.bipartite import biadjacency_matrix
import faiss
import numpy as np
import pandas as pd
import scipy.sparse
import seaborn as sns
from anndata import AnnData
import sklearn

import scanpy as sc
import anndata as ad

from tqdm.auto import tqdm


def make_dist_bins(dist, bins):
    r"""
    ``bins`` are in KB
    """
    labels = [f"{bins[i]}-{bins[i+1]} kb" for i in range(len(bins) - 1)]
    bins = np.asarray(bins) * 1e3
    return pd.cut(dist, bins, labels=labels, include_lowest=True)

# def boxplot(x=None, y=None, hue=None, data=None, bins_factor=2):
#     """
#     Box plot with marginal frequency distribution histogram
#     """
#     # Ensure the specified columns exist in the DataFrame
#     assert x in data and y in data and hue in data

#     # Make a copy of the DataFrame to avoid modifying the original data
#     plot_data = data.copy(deep=False)

#     # Convert columns to appropriate data types
#     if not pd.api.types.is_categorical_dtype(plot_data[x]):
#         plot_data[x] = plot_data[x].astype("category")
#     if not pd.api.types.is_categorical_dtype(plot_data[hue]):
#         plot_data[hue] = plot_data[hue].astype("category")
#     plot_data[y] = plot_data[y].astype(float)

#     # Create the JointGrid instance
#     g = sns.JointGrid(x=x, y=y, data=plot_data, height=5)
#     # Hide the right y-axis
#     # g.ax_joint.yaxis.set_visible(False)
#     g.ax_marg_y.set_visible(False)
#     print(g.ax_joint)
#     # Boxplot on joint
#     sns.boxplot(
#         x=x, y=y, hue=hue, data=plot_data,
#         saturation=1.0, showmeans=True,
#         meanprops=dict(marker="^", markerfacecolor="white", markeredgecolor="black"),
#         boxprops=dict(edgecolor="black"), medianprops=dict(color="black"),
#         whiskerprops=dict(color="black"), capprops=dict(color="black"),
#         flierprops=dict(marker=".", markerfacecolor="black", markeredgecolor="none", markersize=3),
#         ax=g.ax_joint
#     )


#     # Group data by 'x' and 'hue' columns and calculate fractions
#     # grouped_data = plot_data.groupby([x, hue])[y].size().groupby(level=0).apply(lambda x: x / float(x.sum())).reset_index(name='frac')
#     grouped_data = plot_data.groupby([x, hue])[y].size().groupby(level=0).apply(lambda x: x / float(x.sum()))

#     # Initialize bottom array for the stacked bar chart
#     bottom = np.zeros(len(plot_data[x].cat.categories))

#     # Loop through each level of the hue and plot the bars
#     for name, group in grouped_data.groupby(hue):
#         g.ax_marg_x.bar(group[x], group['frac'], bottom=bottom, label=name, width=0.7, edgecolor="black")
#         bottom += group['frac'].values

#     # Add a legend if there are hue levels
#     if len(plot_data[hue].cat.categories) > 1:
#         g.ax_marg_x.legend(title=hue)

#     return g
def boxplot(x=None, y=None, hue=None, data=None):
    r"""
    Box plot with marginal distributions
    """
    assert x in data and y in data and hue in data
    data = data.copy(deep=False)
    if not pd.api.types.is_categorical_dtype(data[x]):
        data[x] = data[x].astype("category")
    if not pd.api.types.is_categorical_dtype(data[hue]):
        data[hue] = data[hue].astype("category")
    data[y] = data[y].astype(float)

    g = sns.JointGrid(x=x, y=y, data=data, height=5)
    sns.boxplot(
        x=x, y=y, hue=hue, data=data,
        saturation=1.0, showmeans=True,
        meanprops=dict(marker="^", markerfacecolor="white", markeredgecolor="black"),
        boxprops=dict(edgecolor="black"), medianprops=dict(color="black"),
        whiskerprops=dict(color="black"), capprops=dict(color="black"),
        flierprops=dict(marker=".", markerfacecolor="black", markeredgecolor="none", markersize=3),
        ax=g.ax_joint
    )
    sns.kdeplot(
        y=y, hue=hue, data=data,
        common_norm=False, shade=True, legend=False, ax=g.ax_marg_y
    )
    data = data.groupby(x)[hue].value_counts(normalize=True).rename("frac").reset_index()
    bottom = np.zeros(data[x].cat.categories.size)
    for _, d in data.groupby(hue):
        g.ax_marg_x.bar(d[x], d["frac"], bottom=bottom, width=0.7, edgecolor="black")
        bottom += d["frac"]
    return g


def metacell_corr(rna, atac, use_rep, n_meta=200, skeleton=None, method="spr"):
    print("Clustering metacells...")
    rna_agg, atac_agg = get_metacells_paired(rna, atac, use_rep, n_meta=n_meta)
    print("Computing correlation...")
    return _metacell_corr(rna_agg, atac_agg, skeleton=skeleton, method=method)

def get_metacells_paired(rna, atac, use_rep, n_meta=150):
    X_pca = np.ascontiguousarray(rna.obsm[use_rep])
    kmeans = faiss.Kmeans(X_pca.shape[1], n_meta, gpu=False, seed=0)

    kmeans.train(X_pca)
    _, rna.obs["metacell"] = kmeans.index.search(X_pca, 1)
    atac.obs["metacell"] = rna.obs["metacell"].to_numpy()
    # # 使用 sklearn 的 KMeans
    # kmeans = KMeans(n_clusters=n_meta, random_state=0)
    # kmeans.fit(X_pca)

    # # 获取每个样本的簇标签
    # rna.obs["metacell"] = kmeans.labels_
    # atac.obs["metacell"] = rna.obs["metacell"].to_numpy()


    rna_agg = aggregate_obs(rna, "metacell")
    atac_agg = aggregate_obs(atac, "metacell")
    common_metacells = np.intersect1d(rna_agg.obs_names, atac_agg.obs_names)
    rna_agg = rna_agg[common_metacells, :].copy()
    atac_agg = atac_agg[common_metacells, :].copy()
    return rna_agg, atac_agg



def metacell_regr(rna, atac, use_rep, n_meta=200, skeleton=None, model="Lasso", **kwargs):
    print("Clustering metacells...")
    rna_agg, atac_agg = get_metacells_paired(rna, atac, use_rep, n_meta=n_meta)
    print("Computing regression...")
    return _metacell_regr(rna_agg, atac_agg, skeleton=skeleton, model=model, **kwargs)


def aggregate_obs(
        adata, by: str, X_agg = "sum",
        obs_agg= None,
        obsm_agg = None,
        layers_agg = None
) -> AnnData:
    r"""
    Aggregate obs in a given dataset by certain categories

    Parameters
    ----------
    adata
        Dataset to be aggregated
    by
        Specify a column in ``adata.obs`` used for aggregation,
        must be discrete.
    X_agg
        Aggregation function for ``adata.X``, must be one of
        ``{"sum", "mean", ``None``}``. Setting to ``None`` discards
        the ``adata.X`` matrix.
    obs_agg
        Aggregation methods for ``adata.obs``, indexed by obs columns,
        must be one of ``{"sum", "mean", "majority"}``, where ``"sum"``
        and ``"mean"`` are for continuous data, and ``"majority"`` is for
        discrete data. Fields not specified will be discarded.
    obsm_agg
        Aggregation methods for ``adata.obsm``, indexed by obsm keys,
        must be one of ``{"sum", "mean"}``. Fields not specified will be
        discarded.
    layers_agg
        Aggregation methods for ``adata.layers``, indexed by layer keys,
        must be one of ``{"sum", "mean"}``. Fields not specified will be
        discarded.

    Returns
    -------
    aggregated
        Aggregated dataset
    """
    obs_agg = obs_agg or {}
    obsm_agg = obsm_agg or {}
    layers_agg = layers_agg or {}

    by = adata.obs[by]
    agg_idx = pd.Index(by.cat.categories) \
        if pd.api.types.is_categorical_dtype(by) \
        else pd.Index(np.unique(by))
    agg_sum = scipy.sparse.coo_matrix((
        np.ones(adata.shape[0]), (
            agg_idx.get_indexer(by),
            np.arange(adata.shape[0])
        )
    )).tocsr()
    agg_mean = agg_sum.multiply(1 / agg_sum.sum(axis=1))

    agg_method = {
        "sum": lambda x: agg_sum @ x,
        "mean": lambda x: agg_mean @ x,
        "majority": lambda x: pd.crosstab(by, x).idxmax(axis=1).loc[agg_idx].to_numpy()
    }

    X = agg_method[X_agg](adata.X) if X_agg and adata.X is not None else None
    obs = pd.DataFrame({
        k: agg_method[v](adata.obs[k])
        for k, v in obs_agg.items()
    }, index=agg_idx.astype(str))
    obsm = {
        k: agg_method[v](adata.obsm[k])
        for k, v in obsm_agg.items()
    }
    layers = {
        k: agg_method[v](adata.layers[k])
        for k, v in layers_agg.items()
    }
    for c in obs:
        if pd.api.types.is_categorical_dtype(adata.obs[c]):
            obs[c] = pd.Categorical(obs[c], categories=adata.obs[c].cat.categories)
    return AnnData(
        X=X, obs=obs, var=adata.var,
        obsm=obsm, varm=adata.varm, layers=layers,
        dtype=None if X is None else X.dtype
    )


def densify(arr) -> np.ndarray:
    r"""
    Convert a matrix to dense regardless of original type.

    Parameters
    ----------
    arr
        Input array (either sparse or dense)

    Returns
    -------
    densified
        Densified array
    """
    if scipy.sparse.issparse(arr):
        return arr.toarray()
    if isinstance(arr, np.ndarray):
        return arr
    return np.asarray(arr)


def _metacell_corr(
        *adatas, skeleton = None, method = "spr",
        prep_fns = None
) -> nx.Graph:
    if skeleton is None:
        raise ValueError("Missing required argument `skeleton`!")
    if set.intersection(*(set(adata.var_names) for adata in adatas)):
        raise ValueError("Overlapping features are currently not supported!")
    prep_fns = prep_fns or [None] * len(adatas)
    if not len(prep_fns) == len(adatas):
        raise ValueError("Length of `prep_fns` must match the number of datasets!")
    for adata, prep_fn in zip(adatas, prep_fns):
        if prep_fn:
            prep_fn(adata)
    adata = ad.concat(adatas, axis=1)
    edgelist = nx.to_pandas_edgelist(skeleton)
    source = adata.var_names.get_indexer(edgelist["source"])
    target = adata.var_names.get_indexer(edgelist["target"])
    X = densify(adata.X.T)
    if method == "spr":
        X = np.array([scipy.stats.rankdata(x) for x in X])
    elif method != "pcc":
        raise ValueError(f"Unrecognized method: {method}!")
    mean = X.mean(axis=1)
    meansq = np.square(X).mean(axis=1)
    std = np.sqrt(meansq - np.square(mean))
    edgelist["corr"] = np.array([
        ((X[s] * X[t]).mean() - mean[s] * mean[t]) / (std[s] * std[t])
        for s, t in zip(source, target)
    ])
    return nx.from_pandas_edgelist(edgelist, edge_attr=True, create_using=type(skeleton))


def _metacell_regr(
        *adatas: AnnData, skeleton: nx.DiGraph = None,
        model: str = "Lasso", **kwargs
) -> nx.DiGraph:
    if skeleton is None:
        raise ValueError("Missing required argument `skeleton`!")
    for adata in adatas:
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
    if set.intersection(*(set(adata.var_names) for adata in adatas)):
        raise ValueError("Overlapping features are currently not supported!")
    adata = ad.concat(adatas, axis=1)

    targets = [node for node, in_degree in skeleton.in_degree() if in_degree]
    biadj = biadjacency_matrix(
        skeleton, adata.var_names, targets, weight=None
    ).astype(bool).T.tocsr()
    X = densify(adata.X)
    Y = densify(adata[:, targets].X.T)
    coef = []
    model = getattr(sklearn.linear_model, model)
    for target, y, mask in tqdm(zip(targets, Y, biadj), total=len(targets), desc="metacell_regr"):
        # print(f'target:{target},y:{y.shape},mask:{mask.shape}')
        # print(f'mask_indices:{mask.indices}')
        X_ = X[:, mask.indices]
        lm = model(**kwargs).fit(X_, y)
        # print(f'X_: {X_.shape}, lm.coef_: {lm.coef_.shape}')
        coef.append(pd.DataFrame({
            "source": adata.var_names[mask.indices],
            "target": target,
            "regr": lm.coef_
        }))
    coef = pd.concat(coef)
    return nx.from_pandas_edgelist(coef, edge_attr=True, create_using=type(skeleton)),coef

