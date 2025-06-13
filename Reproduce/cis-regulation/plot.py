r"""
Plotting functions
"""

from typing import Callable, List, Optional, Union

import matplotlib.axes as ma
import numpy as np
import pandas as pd
import seaborn as sns
import scanpy as sc
import sklearn.metrics
from matplotlib import rcParams



def roc(
        true: np.ndarray, pred: np.ndarray, max_points: int = 500,
        ax: Optional[ma.Axes] = None, **kwargs
) -> ma.Axes:
    r"""
    Plot an ROC curve

    Parameters
    ----------
    true
        True labels
    pred
        Prediction values
    max_points
        Maximal number of points on the ROC curve, beyond which the points
        are equidistantly subsampled.
    ax
        Existing axes to plot on
    **kwargs
        Additional keyword arguments passed to :func:`seaborn.lineplot`

    Returns
    -------
    ax
        Plot axes
    """
    fpr, tpr, _ = sklearn.metrics.roc_curve(true, pred)
    idx = np.linspace(
        0, fpr.size, min(fpr.size, max_points), endpoint=False
    ).round().astype(int)
    idx[-1] = fpr.size - 1  # Always keep the last point
    data = pd.DataFrame({"FPR": fpr[idx], "TPR": tpr[idx]})
    ax = sns.lineplot(x="FPR", y="TPR", data=data, ax=ax, **kwargs)
    return ax


def prc(
        true: np.ndarray, pred: np.ndarray, max_points: int = 500,
        ax: Optional[ma.Axes] = None, **kwargs
) -> ma.Axes:
    r"""
    Plot a precision-recall curve

    Parameters
    ----------
    true
        True labels
    pred
        Prediction values
    max_points
        Maximal number of points on the precision-recall curve, beyond which
        the points are equidistantly subsampled.
    ax
        Existing axes to plot on
    **kwargs
        Additional keyword arguments passed to :func:`seaborn.lineplot`

    Returns
    -------
    ax
        Plot axes
    """
    prec, rec, _ = sklearn.metrics.precision_recall_curve(true, pred)
    idx = np.linspace(
        0, prec.size, min(prec.size, max_points), endpoint=False
    ).round().astype(int)
    idx[-1] = prec.size - 1  # Always keep the last point
    data = pd.DataFrame({"Precision": prec[idx], "Recall": rec[idx]})
    ax = sns.lineplot(x="Recall", y="Precision", data=data, ax=ax, **kwargs)
    return ax
