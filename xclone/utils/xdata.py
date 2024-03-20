# xdata.py

import logging
import numpy as np
import scipy as sp


def check_sanity_layer(xdata, layer):
    state = 0

    mtx = None
    if sp.sparse.issparse(xdata.layers[layer]):
        mtx = xdata.layers[layer].copy().A
    else:
        mtx = xdata.layers[layer].copy()
    
    # detect nan Value
    nan_count = np.isnan(mtx).sum()
    if nan_count > 0:
        logging.warning("NaN values in layer '%s'!" % layer)
        state |= (1<<0)
    
    # detect negative Value
    if np.any(mtx < 0):
        logging.warning("negative values in layer '%s'!" % layer)
        state |= (1<<1)
    
    return(state)


def check_unanno_cells(xdata, remove_unanno = True, alt_cell_type = "unannotated", verbose = True):
    cell_anno_key = "cell_type"
    if remove_unanno:
        valid_cells = xdata.obs[cell_anno_key] == xdata.obs[cell_anno_key]
        new_xdata = xdata[valid_cells, :].copy()
        if verbose:
            n_cells = xdata.shape[0]
            logging.info("filter out %d (out of %d) cells." %
                (n_cells - valid_cells.sum(), n_cells))
    else:
        new_xdata = xdata.copy()
        new_xdata.obs[cell_anno_key].fillna(alt_cell_type, inplace = True)
    return new_xdata


def remove_XY(xdata):
    flag = ~(xdata.var["chrom"].isin(["X", "Y"]))
    return xdata[:, flag].copy()
