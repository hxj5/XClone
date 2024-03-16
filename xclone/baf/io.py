# io.py - BAF data input/output.

import anndata as ad
import os
import pandas as pd
import scipy as sp
from scipy import io
from scipy import sparse


def load_data(data_dir, file_prefix = "xcltk"):
    samples = load_samples(os.path.join(data_dir, file_prefix + ".samples.tsv"))
    regions = load_regions(os.path.join(data_dir, file_prefix + ".region.tsv"))
    AD_mtx = load_matrix(os.path.join(data_dir, file_prefix + ".AD.mtx"))
    DP_mtx = load_matrix(os.path.join(data_dir, file_prefix + ".DP.mtx"))
    OTH_mtx = load_matrix(os.path.join(data_dir, file_prefix + ".OTH.mtx"))
    
    adata = ad.AnnData(
        X = AD_mtx, 
        obs = regions,
        var = samples)
    adata.layers["AD"] = AD_mtx
    adata.layers["DP"] = DP_mtx
    adata.layers["OTH"] = OTH_mtx

    adata = adata.transpose()      # cell x region
    return(adata)


def save_data(adata, out_dir, file_prefix = "xcltk"):
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    
    save_samples(adata.obs, 
        fn = os.path.join(out_dir, file_prefix + ".samples.tsv"))
    save_regions(adata.var, 
        fn = os.path.join(out_dir, file_prefix + ".region.tsv"))
    save_matrix(adata.layers["AD"], os.path.join(out_dir, file_prefix + ".AD.mtx"))
    save_matrix(adata.layers["DP"], os.path.join(out_dir, file_prefix + ".DP.mtx"))
    save_matrix(adata.layers["OTH"], os.path.join(out_dir, file_prefix + ".OTH.mtx"))


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
    df.columns = ["chrom", "start", "end", "region"]
    return(df)


def save_regions(df, fn):
    df.to_csv(fn, sep = "\t", header = False, index = False)


def load_samples(fn):
    df = pd.read_csv(fn, header = None)
    df.columns = ["cell"]
    return(df)


def save_samples(df, fn):
    df.to_csv(fn, sep = "\t", header = False, index = False)

