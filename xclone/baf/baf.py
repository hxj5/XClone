# baf.py - main functions

import logging
import multiprocessing
import numpy as np
from scipy import sparse

from .phasing import Local_Phasing, Global_Phasing
from ..utils.xdata import check_sanity_layer, check_unanno_cells, remove_XY


def check_sanity(xdata, verbose = True):
    if verbose:
        logging.info("begin ...")

    state = 0

    if verbose:
        logging.info("check unannotated cells ...")
    xdata = check_unanno_cells(
        xdata, remove_unanno = True, verbose = verbose
    )
    if verbose:
        logging.info("remove chromosome X and Y ...")
    xdata = remove_XY(xdata)

    if verbose:
        logging.info("check sanity of layer 'AD' ...")
    st = check_sanity_layer(xdata, "AD")
    state |= st
    if verbose:
        logging.info("layer 'AD' state %d." % st)

    if verbose:
        logging.info("check sanity of layer 'DP' ...")
    st = check_sanity_layer(xdata, "DP")
    state |= st
    if verbose:
        logging.info("layer 'DP' state %d." % st)

    return((state, xdata))


def do_local_phasing(
    xdata, 
    chrom_lst = None, 
    phasing_len = 100, 
    reg_n_proc = 1,
    bin_n_proc = 1,
    feature_mode = "gene",  # or "block"
    var_add = None,
    verbose = False
):
    """
    Func:
    phasing_len: default for 100 genes.
    """
    if verbose:
        logging.info("begin ...")

    reg_key = "chrom"
    chrom_lst_raw = chrom_lst
    if chrom_lst is None:
        chrom_lst = xdata.var[reg_key].drop_duplicates(keep = "first")

    if reg_n_proc > 1:
        if verbose:
            logging.info("using multi-processing: n_proc = %d ..." % reg_n_proc)
        result = []
        pool = multiprocessing.Pool(processes = reg_n_proc)
        for chrom in chrom_lst:
            reg_xdata = xdata[:, xdata.var[reg_key] == chrom]
            AD_reg = reg_xdata.layers["AD"].T
            DP_reg = reg_xdata.layers["DP"].T
            result.append(pool.apply_async(
                do_local_phasing_reg, 
                (chrom, AD_reg, DP_reg, phasing_len, bin_n_proc, verbose),
                callback = None)
            )
        pool.close()
        pool.join()
        result = [res.get() for res in result]
    else:
        if verbose:
            logging.info("using single process ...")
        result = []
        for chrom in chrom_lst:
            reg_xdata = xdata[:, xdata.var[reg_key] == chrom]
            AD_reg = reg_xdata.layers["AD"].T
            DP_reg = reg_xdata.layers["DP"].T
            RV_reg = do_local_phasing_reg(
                chrom, AD_reg, DP_reg,
                phasing_len = phasing_len, 
                n_proc = bin_n_proc,
                verbose = verbose)
            result.append(RV_reg)

    ## process the data from the list to a dataset
    for i, RV_reg in zip(range(len(result)), result):
        if i == 0:
            AD_phased = RV_reg["AD_phased"]
            theta_bin = RV_reg["theta_bin"]
            bin_idx = RV_reg["bin_idx"]
            bin_idx_lst = RV_reg["bin_idx_lst"]
            ad_bin_softcnt = RV_reg["ad_bin_softcnt"]
            ad_bin = RV_reg["ad_bin"]
            dp_bin = RV_reg["dp_bin"]
            allele_flip_local = RV_reg["allele_flip_local"]
        else:
            AD_phased = np.vstack((AD_phased, RV_reg["AD_phased"]))
            theta_bin = np.vstack((theta_bin, RV_reg["theta_bin"]))
            bin_idx = np.append(bin_idx, RV_reg["bin_idx"])
            bin_idx_lst = np.append(bin_idx_lst, RV_reg["bin_idx_lst"])
            ad_bin_softcnt = np.vstack((ad_bin_softcnt, RV_reg["ad_bin_softcnt"]))
            ad_bin = np.vstack((ad_bin, RV_reg["ad_bin"]))
            dp_bin = np.vstack((dp_bin, RV_reg["dp_bin"]))
            allele_flip_local = np.append(allele_flip_local, RV_reg["allele_flip_local"])
 
    update_xdata = None
    if chrom_lst_raw is None:
        update_xdata = xdata.copy()
    else:
        update_xdata = xdata[:, xdata.var[reg_key].isin(chrom_lst_raw)].copy()
    update_xdata.layers["AD_phased"] = AD_phased.T        # to be cell x region
    update_xdata.var["bin_idx"] = bin_idx_lst
    update_xdata.var["allele_flip_local"] = allele_flip_local
    try:
        update_xdata.var["allele_flip_local"].replace({0: False, 1: True}, inplace = True)
    except Exception as e:
        logging.warning(str(e))
    else:
        logging.info("get allele flip status from local phasing.")  
    update_xdata.obsm["theta_bin"] = theta_bin.T

    return update_xdata


def do_local_phasing_reg(
    reg_idx,
    AD,     # feature x cell
    DP, 
    phasing_len = 100,
    n_proc = 1, 
    verbose = False
):
    """
    Func:
    region: default chr_based. 
    """
    if verbose:
        logging.info("begin to process region '%s' ..." % reg_idx)

    n_bins = int(AD.shape[0] / phasing_len)
    last_bin_len = AD.shape[0] % phasing_len
    
    ## do phasing
    if n_proc > 1:
        if verbose:
            logging.debug("using multi-processing, n_proc = %d ..." % n_proc)
        result = []
        pool = multiprocessing.Pool(processes = n_proc)
        if last_bin_len == 0:
            for ib in range(n_bins):
                idx = range(ib*phasing_len, (ib+1)*phasing_len)
                result.append(pool.apply_async(
                    do_local_phasing_bin,
                    (AD[idx, :], DP[idx, :], reg_idx, ib, verbose),
                    callback = None)
                )
        else:
            for ib in range(n_bins + 1):
                if ib == n_bins:
                    idx = range(-last_bin_len, 0)
                else:
                    idx = range(ib*phasing_len, (ib+1)*phasing_len)
                result.append(pool.apply_async(
                    do_local_phasing_bin,
                    (AD[idx, :], DP[idx, :], reg_idx, ib, verbose),
                    callback = None)
                )
        pool.close()
        pool.join()
        result = [res.get() for res in result]
    else:
        if verbose:
            logging.debug("using single process ...")
        result = []
        if last_bin_len == 0:
            for ib in range(n_bins):
                idx = range(ib*phasing_len, (ib+1)*phasing_len)
                RV_bin = do_local_phasing_bin(AD[idx, :], DP[idx, :], reg_idx, ib, verbose)
                result.append(RV_bin)
        else:
            for ib in range(n_bins + 1):
                if ib == n_bins:
                    idx = range(-last_bin_len, 0)
                else:
                    idx = range(ib*phasing_len, (ib+1)*phasing_len)
                RV_bin = do_local_phasing_bin(AD[idx, :], DP[idx, :], reg_idx, ib, verbose)
                result.append(RV_bin)

    ## resolve result
    for i, RV_bin in zip(range(len(result)), result):
        if i == 0:
            AD_phased = RV_bin["AD_phased"]
            ad_bin_softcnt = RV_bin["ad_bin_softcnt"]
            ad_bin = RV_bin["ad_bin"]
            dp_bin = RV_bin["dp_bin"]
            theta_bin = RV_bin["theta_bin"]
            bin_idx = RV_bin["bin_idx"]
            bin_idx_lst = RV_bin["bin_idx_lst"]
            allele_flip_local = RV_bin["flip"]
        else:
            AD_phased = np.vstack((AD_phased, RV_bin["AD_phased"]))
            ad_bin_softcnt = np.vstack((ad_bin_softcnt, RV_bin["ad_bin_softcnt"]))
            ad_bin = np.vstack((ad_bin, RV_bin["ad_bin"]))
            dp_bin = np.vstack((dp_bin, RV_bin["dp_bin"]))
            theta_bin = np.vstack((theta_bin, RV_bin["theta_bin"]))
            bin_idx = np.append(bin_idx, RV_bin["bin_idx"])
            bin_idx_lst = np.append(bin_idx_lst, RV_bin["bin_idx_lst"])
            allele_flip_local = np.append(allele_flip_local, RV_bin["flip"])

    ## resolve results for global phasing input
    RV_reg = {}
    RV_reg["AD_phased"] = AD_phased
    RV_reg["bin_idx"] = bin_idx
    RV_reg["ad_bin_softcnt"] = ad_bin_softcnt
    RV_reg["ad_bin"] = ad_bin
    RV_reg["dp_bin"] = dp_bin
    RV_reg["theta_bin"] = theta_bin
    RV_reg["allele_flip_local"] = allele_flip_local
    RV_reg["bin_idx_lst"] = bin_idx_lst    # for global phasing record

    return RV_reg


def do_local_phasing_bin(AD, DP, reg_idx, bin_idx, verbose = False):
    """
    Func:
    bin: default for 100 genes
    """
    if verbose:
        logging.debug("processing region '%s' bin '%d' ..." % (reg_idx, bin_idx))

    ## Z: allele flipping probability
    ad_sum, ad_sum1, dp_sum, Z, thetas, _logLik_new = Local_Phasing(AD, DP, verbose = verbose)
    if np.isnan(_logLik_new):
        if verbose:
            logging.warning("logLik is NaN, region '%s' bin '%d' ..." % (reg_idx, bin_idx))

    is_flip = np.argmax(Z, axis = 1)
    
    ## important
    BD = DP - AD
    AD_phased = AD * np.expand_dims(1 - is_flip, axis = 1) + \
        BD * np.expand_dims(is_flip, axis = 1)

    RV = {}
    RV["bin_idx"] = bin_idx
    RV['bin_idx_lst'] = np.repeat(bin_idx, AD_phased.shape[0])
    RV['Z'] = Z
    RV['flip'] = is_flip
    RV["AD_phased"] = AD_phased
    RV['ad_bin_softcnt'] = ad_sum[0, :]
    RV['ad_bin'] = ad_sum1[0, :]
    RV['dp_bin'] = dp_sum[0, :]
    RV['theta_bin'] = np.array(thetas)[:, 0]
    RV['logLik'] = _logLik_new
    return RV
