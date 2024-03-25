# baf_main.py

import logging
from sys import stdout, stderr
from .baf.baf import check_sanity, do_local_phasing, feature2bin

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
        logging.debug("xdata after check_sanity:")
        logging.debug(str(xdata))

    if verbose:
        logging.info("do local phasing ...")

    xdata, local_stat = do_local_phasing(
        xdata, 
        chrom_lst = ["8", "10", "18"],
        #chrom_lst = ["10"],
        phasing_len = 100,
        reg_n_proc = 1,
        bin_n_proc = 1,
        verbose = verbose)
    
    if verbose:
        logging.debug("xdata after local_phasing:")
        logging.debug(str(xdata))
    
    if verbose:
        logging.info("merge features into bin ...")
    
    merge_xdata = feature2bin(
        xdata = xdata,
        stat = local_stat,
        feature_mode = "gene",
        var_add = None,
        verbose = verbose
    )

    if verbose:
        logging.debug("xdata after feature2bin:")
        logging.debug(str(xdata))
        logging.debug("merge_xdata after feature2bin:")
        logging.debug(str(merge_xdata))

    return((ret, xdata))
