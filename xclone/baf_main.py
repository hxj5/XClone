# baf_main.py

from sys import stdout, stderr
from .baf.baf import check_sanity

def run_baf(xdata, verbose = True):
    func = "run_baf"
    ret = -1

    if verbose:
        stdout.write("[I::%s] begin...\n" % func)

    if verbose:
        stdout.write("[I::%s] check sanity ...\n" % func)

    ret, xdata = check_sanity(xdata, verbose = verbose)

    return((ret, xdata))
