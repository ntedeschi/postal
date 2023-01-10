import anndata as ad
from pathlib import Path
import scanpy as sc


def leiden(adata, neighborhood, resolution, use_rep, key_added):
    sc.pp.neighbors(adata, use_rep=use_rep, n_neighbors=neighborhood)
    sc.tl.leiden(adata, key_added=key_added, resolution=resolution)
