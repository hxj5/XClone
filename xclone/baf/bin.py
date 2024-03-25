# bin.py - processing bin_xdata

import logging
import matplotlib.pylab as plt
import numpy as np
import scipy as sp


def cal_cell_BAF(xdata, AD_key = "AD", DP_key = "DP", BAF_key = "BAF"):
    """
    xdata: anndata.
    AD_key: Layer name for AD used.
    DP_key: Layer name for DP used.
    BAF_key: Layer name for the calculated BAF.
    """
    if sp.sparse.issparse(xdata.layers[AD_key]):
        AD = xdata.layers[AD_key].A
    else:
        AD = xdata.layers[AD_key]
    if sp.sparse.issparse(xdata.layers[DP_key]):
        DP = xdata.layers[DP_key].A
    else:
        DP = xdata.layers[DP_key]
    
    xdata.layers[BAF_key] = AD / DP
    return xdata


def cap_extreme_count(
    xdata,
    quantile = 0.99,
    backup = True,
    verbose = False
):
    """
    remove extreme counts influence.
    """
    ad_counts = xdata.layers["ad_bin"].A.sum(axis = 0)
    dp_counts = xdata.layers["dp_bin"].A.sum(axis = 0)

    # TODO: the figure will be overwritten by next step ("after capping").
    if verbose:
        logging.info("distribution of ad_counts and dp_counts before capping:")
        plt.plot(ad_counts)
        plt.plot(dp_counts)

    cutoff = np.quantile(dp_counts, quantile)
    flag = dp_counts > cutoff
    ratio_replace = dp_counts[flag] / cutoff

    if backup:
        xdata.layers["ad_bin_backup"] = xdata.layers["ad_bin"].copy()
        xdata.layers["dp_bin_backup"] = xdata.layers["dp_bin"].copy()
    
    # mask to, and insert your replacement values:
    xdata.layers["ad_bin"][:, flag] = np.ceil(xdata.layers["ad_bin"][:, flag] / ratio_replace)
    xdata.layers["dp_bin"][:, flag] = np.ceil(xdata.layers["dp_bin"][:, flag] / ratio_replace)

    if verbose:
        ad_counts = xdata.layers["ad_bin"].A.sum(axis = 0)
        dp_counts = xdata.layers["dp_bin"].A.sum(axis = 0)
        logging.info("distribution of ad_counts and dp_counts after capping:")
        plt.plot(ad_counts)
        plt.plot(dp_counts)

    return xdata
