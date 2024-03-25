# baf.py - main functions

import anndata as ad
import logging
import multiprocessing
import numpy as np
import pandas as pd
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


def do_global_phasing(xdata, bin_xdata, verbose = False):
    """
    TODO
    - test the function.
    - get filp status for bins (genes).
    - add plots for global phasing part (Xlayer = "BAF_globbal_phased" or global_phased = True).
    """
    if verbose:
        logging.info("begin ...")

    ## Global_Phasing
    is_flips, distances, theta_new = Global_Phasing(xdata.obsm["theta_bin"].T)

    ## apply global phasing method on original AD_phased
    bin_value_counts = xdata.var["bin_idx_cum"].value_counts()
    bin_items = xdata.var["bin_idx_cum"].drop_duplicates(keep = "first")
    
    ## check is_flips the same length with bin_item
    if is_flips.shape[0] != bin_items.shape[0]:
        raise ValueError("length of is_flips (%d) is different from bin_items (%d)." %
            (is_flips.shape[0], bin_items.shape[0]))

    gene_flip = None
    for i, bin_id in enumerate(bin_items):
        if i == 0:
            gene_flip = np.repeat(is_flips[i], bin_value_counts[bin_id])
        else:
            gene_flip = np.append(gene_flip, np.repeat(is_flips[i], bin_value_counts[bin_id]))

    # get gene_flip for AD_Phased
    AD_phased = xdata.layers["AD_phased"]
    BD_phased = xdata.layers["DP"] - xdata.layers["AD_phased"]
    AD_global_phased = AD_phased + 0
    AD_global_phased[:, gene_flip] = BD_phased[: , gene_flip] + 0

    new_xdata = xdata.copy()
    new_xdata.layers["AD_global_phased"] = AD_global_phased
    new_xdata.var["allele_flip_global"] = gene_flip

    if verbose:
        logging.info("get allele flip status from global phasing.")

    if {'allele_flip_local', 'allele_flip_global'}.issubset(new_xdata.var.columns):
        new_xdata.var["allele_flip"] =  new_xdata.var["allele_flip_local"] ^ new_xdata.var["allele_flip_global"]
        logging.info("get final allele flip status.")

    ## apply global phasing method on AD_phased bins
    ad_bin_softcnt = bin_xdata.layers["ad_bin_softcnt"]
    ad_bin = bin_xdata.layers["ad_bin"]
    dp_bin = bin_xdata.layers["dp_bin"]
    
    bd_bin_softcnt = dp_bin - ad_bin_softcnt
    ad_bin_softcnt_phased = ad_bin_softcnt + 0
    ad_bin_softcnt_phased[: , is_flips] = bd_bin_softcnt[:, is_flips] + 0

    bd_bin = dp_bin - ad_bin
    ad_bin_phased = ad_bin + 0
    ad_bin_phased[: , is_flips] = bd_bin[:, is_flips] + 0

    # the `ad_bin_softcnt_phased` is weighted counts(bin counts is weighted)(soft phasing).
    # `ad_bin_phased`` is hard phasing counts.
    bin_xdata.layers["ad_bin_softcnt_phased"] = ad_bin_softcnt_phased
    bin_xdata.layers["ad_bin_phased"] = ad_bin_phased
    
    return new_xdata, bin_xdata


def do_local_phasing(
    xdata, 
    chrom_lst = None, 
    phasing_len = 100, 
    reg_n_proc = 1,
    bin_n_proc = 1,
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
    AD_phased = theta_bin = bin_idx = bin_idx_lst = None
    ad_bin_softcnt = ad_bin = dp_bin = allele_flip_local = None
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

    RV = {}
    RV["reg_key"] = reg_key
    RV["phasing_len"] = phasing_len
    #RV["AD_phased"] = AD_phased
    #RV["theta_bin"] = theta_bin
    #RV["bin_idx"] = bin_idx
    RV["bin_idx_lst"] = bin_idx_lst
    RV["ad_bin_softcnt"] = ad_bin_softcnt
    RV["ad_bin"] = ad_bin
    RV["dp_bin"] = dp_bin
    #RV["allele_flip_local"] = allele_flip_local

    new_xdata = None
    if chrom_lst_raw is None:
        new_xdata = xdata.copy()
    else:
        new_xdata = xdata[:, xdata.var[reg_key].isin(chrom_lst_raw)].copy()
    new_xdata.layers["AD_phased"] = AD_phased.T        # to be cell x feature
    new_xdata.var["bin_idx"] = bin_idx_lst
    new_xdata.var["allele_flip_local"] = allele_flip_local
    try:
        new_xdata.var["allele_flip_local"].replace({0: False, 1: True}, inplace = True)
    except Exception as e:
        logging.warning(str(e))
    else:
        logging.info("get allele flip status from local phasing.")
    new_xdata.obsm["theta_bin"] = theta_bin.T

    return new_xdata, RV


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
    RV = {}
    RV["AD_phased"] = AD_phased
    RV["bin_idx"] = bin_idx
    RV["ad_bin_softcnt"] = ad_bin_softcnt
    RV["ad_bin"] = ad_bin
    RV["dp_bin"] = dp_bin
    RV["theta_bin"] = theta_bin
    RV["allele_flip_local"] = allele_flip_local
    RV["bin_idx_lst"] = bin_idx_lst    # for global phasing record

    return RV


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


def feature2bin(
    xdata,
    stat,
    feature_mode = "gene",  # or "block"
    var_add = None,
    verbose = False
):
    # 
    def get_bin_features(bin_var, group_key = "bin_idx_cum"):
        feature_lst = []
        feature_dict = {}
        groups = bin_var.groupby(group_key).groups
        for key, idx in groups.items():
            fet_lst = bin_var.loc[idx]["feature"].copy().tolist()
            feature_lst.append(fet_lst)
            feature_dict[key] = fet_lst
        return feature_lst, feature_dict

    if not feature_mode:
        raise ValueError
    if feature_mode.lower() not in ("gene", "block"):
        raise ValueError

    # calculate and assign the bin_idx_cum for bins across all regions (chromosomes).
    # originally implemented in the `process_bin_id` function.
    xv = xdata.var.drop_duplicates(["chrom", "bin_idx"], keep = "first")[["chrom", "bin_idx"]]
    xv["bin_idx_cum"] = range(xv.shape[0])
    xdata.var = pd.merge(xdata.var, xv, on = ["chrom", "bin_idx"], how = "left")

    ##check theta_bin reuslts first and there are nan value
    ## save for visualization
    ad_bin_softcnt = sparse.csr_matrix(stat["ad_bin_softcnt"])
    ad_bin = sparse.csr_matrix(stat["ad_bin"])
    dp_bin = sparse.csr_matrix(stat["dp_bin"])

    ## generate bin_xdata var
    bin_var = xdata.var.copy()
    bin_var.drop_duplicates("bin_idx_cum", keep = "first", inplace = True)
    bin_var.rename(columns={"end": "feature1_end"}, inplace = True)

    xv = xdata.var.drop_duplicates("bin_idx_cum", keep = "last")
    bin_var["end"] = xv["end"].copy().tolist()
    bin_var["bin_end_arm"] = xv["arm"].copy().tolist()
    bin_var["bin_end_band"] = xv["band"].copy().tolist()

    xv =  xdata.var.copy()
    feature_lst, feature_dict = get_bin_features(xv, group_key = "bin_idx_cum")
    bin_var["feature_lst"] = [",".join(x) for x in feature_lst]

    var_keep_lst = [
        "chrom", "start", "end", "feature", "arm", "band", 
        "feature1_end", "bin_end_arm", "bin_end_band", 
        "bin_idx", "bin_idx_cum", 
        "feature_lst"]

    if var_add:
        for v in var_add:
            if v not in var_keep_lst and v in xdata.var.columns:
                var_keep_lst.append(v)

    bin_var = bin_var[var_keep_lst]
    bin_var["bin_features_cnt"] = bin_var["feature_lst"].str.len()

    bin_xdata = ad.AnnData(
        ad_bin.T,
        obs = xdata.obs.copy(),
        var = bin_var)

    ## soft phasing
    bin_xdata.layers["ad_bin_softcnt"] = ad_bin_softcnt.T

    ## hard phasing
    bin_xdata.layers["ad_bin"] = ad_bin.T
    bin_xdata.layers["dp_bin"] = dp_bin.T

    bin_xdata.uns["local_phasing_key"] = stat["reg_key"]
    bin_xdata.uns["local_phasing_len"] = stat["phasing_len"]
    return bin_xdata
