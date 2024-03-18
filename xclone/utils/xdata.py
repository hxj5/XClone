# adata.py

import numpy as np
import scipy as sp
from sys import stdout, stderr


def check_sanity_layer(adata, layer):
    func = "check_sanity_layer"
    state = 0

    mtx = None
    if sp.sparse.issparse(adata.layers[layer]):
        mtx = adata.layers[layer].copy().A
    else:
        mtx = adata.layers[layer].copy()
    
    # detect nan Value
    nan_count = np.isnan(mtx).sum()
    if nan_count > 0:
        stderr.write("[W::%s] NaN values in layer '%s'!\n" % (func, layer))
        state |= (1<<0)
    
    # detect negative Value
    if np.any(mtx < 0):
        stderr.write("[W::%s] negative values in layer '%s'!\n" % (func, layer))
        state |= (1<<1)
    
    return(state)


def check_unanno_cells(adata, remove_unanno = True, alt_cell_type = "unannotated", verbose = True):
    func = "check_unanno_cells"
    cell_anno_key = "cell_type"
    if remove_unanno:
        valid_cells = adata.obs[cell_anno_key] == adata.obs[cell_anno_key]
        new_adata = adata[valid_cells, :].copy()
        if verbose:
            n_cells = adata.shape[0]
            stdout.write("[I::%s] filter out %d (out of %d) cells." %
                (func, n_cells - valid_cells.sum(), n_cells))
    else:
        new_adata = adata.copy()
        new_adata.obs[cell_anno_key].fillna(alt_cell_type, inplace = True)
    return new_adata


def remove_XY(adata):
    flag = ~(adata.var["chrom"].isin(["X", "Y"]))
    return adata[:, flag].copy()


