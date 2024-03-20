# baf.py - main functions

import logging
from ..utils.xdata import check_sanity_layer, check_unanno_cells


def check_sanity(xdata, verbose = True):
    state = 0
    
    if verbose:
        logging.info("begin ...")

    st = check_sanity_layer(xdata, "AD")
    state |= st

    if verbose:
        logging.info("check sanity of layer 'AD', state %d." % st)
        
    st = check_sanity_layer(xdata, "DP")
    state |= st

    if verbose:
        logging.info("check sanity of layer 'DP', state %d." % st)
        
    xdata = check_unanno_cells(
        xdata, remove_unanno = True, verbose = verbose
    )

    return((state, xdata))
