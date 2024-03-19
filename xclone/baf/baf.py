# baf.py - main functions

from sys import stdout, stderr
from ..utils.xdata import check_sanity_layer, check_unanno_cells


def check_sanity(xdata, verbose = True):
    func = "check_sanity"
    state = 0
    
    if verbose:
        stdout.write("[I::%s] begin...\n" % func)

    st = check_sanity_layer(xdata, "AD")
    state |= st

    if verbose:
        stdout.write("[I::%s] check sanity of layer AD, state %d.\n" %
            (func, st))
        
    st = check_sanity_layer(xdata, "DP")
    state |= st

    if verbose:
        stdout.write("[I::%s] check sanity of layer DP, state %d.\n" %
            (func, st))
        
    xdata = check_unanno_cells(
        xdata, remove_unanno = True, verbose = verbose
    )

    return((state, xdata))
