# baf_main.py

import logging
from sys import stdout, stderr
from .baf.baf import check_sanity, do_local_phasing

def run_baf(xdata, verbose = True):
    ret = -1

    if verbose:
        logging.info("begin ...")

    if verbose:
        logging.info("check sanity ...")

    ret, xdata = check_sanity(xdata, verbose = verbose)
    if ret < 0:
        logging.error("check sanity retcode: %d." % ret)
        return((ret, xdata))
    elif ret != 0:
        logging.warning("check sanity retcode: %d." % ret)

    if verbose:
        logging.info("do local phasing ...")

    xdata = do_local_phasing(
        xdata, 
        chrom_lst = ["8", "10", "18"],
        #chrom_lst = ["10"],
        phasing_len = 100,
        reg_n_proc = 1,
        bin_n_proc = 1,
        feature_mode = "gene",   # or "block"
        var_add = None,
        verbose = verbose)

    return((ret, xdata))
