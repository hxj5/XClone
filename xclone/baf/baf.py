# baf.py - main functions

def exclude_XY(adata):
    flag = ~(adata.var["chrom"].isin(["X", "Y"]))
    return adata[:, flag].copy()
