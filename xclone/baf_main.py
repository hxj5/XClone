# baf_main.py

import logging
from sys import stdout, stderr
from .baf.baf import check_sanity

def run_baf(xdata, verbose = True):
    ret = -1

    if verbose:
        logging.info("begin ...")

    if verbose:
        logging.info("check sanity ...")

    ret, xdata = check_sanity(xdata, verbose = verbose)

    return((ret, xdata))
