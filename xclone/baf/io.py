# io.py - BAF data input/output.

import anndata as ad
import os
import pandas as pd
import scipy as sp
from scipy import io
from scipy import sparse

from ..blib.region import format_chrom


def load_xcltk_data(data_dir, cell_anno_fn = None, feature_anno_fn = None, ref_cell_types = None):
    file_prefix = "xcltk"

    if ref_cell_types is None:
        if cell_anno_fn is not None:
            raise ValueError
    elif isinstance(ref_cell_types, list) or isinstance(ref_cell_types, tuple):
        ref_cell_types = list(ref_cell_types)
    else:
        ref_cell_types = [ref_cell_types]

    AD_mtx = load_matrix(os.path.join(data_dir, file_prefix + ".AD.mtx"))
    DP_mtx = load_matrix(os.path.join(data_dir, file_prefix + ".DP.mtx"))
    OTH_mtx = load_matrix(os.path.join(data_dir, file_prefix + ".OTH.mtx"))

    cells = features = None
    if cell_anno_fn is None:
        cells = load_samples(os.path.join(data_dir, file_prefix + ".samples.tsv"))
    else:
        cells = load_cell_anno(cell_anno_fn)
    
    if feature_anno_fn is None:
        features = load_regions(os.path.join(data_dir, file_prefix + ".region.tsv"))
    else:
        features = load_feature_anno(feature_anno_fn)
    
    adata = ad.AnnData(
        X = AD_mtx, 
        obs = features,
        var = cells)
    adata.layers["AD"] = AD_mtx
    adata.layers["DP"] = DP_mtx
    adata.layers["OTH"] = OTH_mtx
    adata.uns["ref_cell_types"] = [ref_cell_types]

    adata = adata.transpose()      # cell x feature
    adata.obs.index = adata.obs.index.astype(str)      # otherwise, anndata will report error of integer index

    return(adata)


def save_xcltk_data(adata, out_dir):
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    file_prefix = "xcltk"
    
    save_samples(adata.obs[["cell"]],
        fn = os.path.join(out_dir, file_prefix + ".samples.tsv"))
    save_regions(adata.var[["chrom", "start", "end", "feature"]],
        fn = os.path.join(out_dir, file_prefix + ".region.tsv"))
    save_matrix(adata.layers["AD"].transpose(), 
        fn = os.path.join(out_dir, file_prefix + ".AD.mtx"))
    save_matrix(adata.layers["DP"].transpose(), 
        fn = os.path.join(out_dir, file_prefix + ".DP.mtx"))
    save_matrix(adata.layers["OTH"].transpose(), 
        fn = os.path.join(out_dir, file_prefix + ".OTH.mtx"))
    save_cell_anno(adata.obs, 
        fn = os.path.join(out_dir, file_prefix + ".cell_anno.tsv"))
    save_feature_anno(adata.var,
        fn = os.path.join(out_dir, file_prefix + ".feature_anno.tsv"))


def load_cell_anno(fn):
    df = pd.read_csv(fn, header = None, sep = "\t")
    df.columns = df.columns.astype(str)
    df.columns.values[:2] = ["cell", "cell_type"]
    return(df)


def save_cell_anno(df, fn):
    df.to_csv(fn, sep = "\t", header = False, index = False)


def load_feature_anno(fn):
    df = pd.read_csv(fn, header = None, sep = "\t")
    df.columns = df.columns.astype(str)
    df.columns.values[:6] = ["chrom", "start", "end", "feature", "arm", "band"]
    return(df)


def save_feature_anno(df, fn):
    df.to_csv(fn, sep = "\t", header = False, index = False)


def load_matrix(fn):
    mtx = None
    try:
        mtx = sp.io.mmread(fn)
    except:
        mtx = io.mmread(fn)
    mtx = mtx.toarray()    # convert from sparse matrix to ndarray to support slicing.
    return(mtx)


def save_matrix(mtx, fn):
    mtx = sparse.csr_matrix(mtx)   # convert from ndarray to sparse matrix to be fully compatible with .mtx format.
    io.mmwrite(fn, mtx)


def load_regions(fn):
    df = pd.read_csv(fn, header = None, sep = "\t")
    df.columns = ["chrom", "start", "end", "feature"]
    df["chrom"] = df["chrom"].map(format_chrom)
    return(df)


def save_regions(df, fn):
    df.to_csv(fn, sep = "\t", header = False, index = False)


def load_samples(fn):
    df = pd.read_csv(fn, header = None)
    df.columns = ["cell"]
    return(df)


def save_samples(df, fn):
    df.to_csv(fn, sep = "\t", header = False, index = False)
