import re
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd


def read_counts(f: Path) -> pd.DataFrame:
    counts = pd.read_csv(f, index_col=0)
    return counts


def read_celldata(f: Path) -> pd.DataFrame:
    celldata = pd.read_csv(f, index_col=0)
    celldata.index = celldata.index.astype(str)
    celldata.index.name = ""
    celldata["center_x_clockwise90"] = celldata["center_y"]
    celldata["center_y_clockwise90"] = celldata["center_x"].max() - celldata["center_x"]
    return celldata


def read_transcripts(f: Path) -> pd.DataFrame:
    trx = pd.read_csv(f)
    trx = trx.rename(columns={"name": "gene_name"})
    return trx


def mkad(
    counts: pd.DataFrame, celldata: pd.DataFrame, trxdata: pd.DataFrame
) -> ad.AnnData:
    X = counts.to_numpy()
    var = pd.DataFrame(index=counts.columns)
    obs = pd.DataFrame(index=counts.index)
    obs.index = obs.index.astype(str)

    adata = ad.AnnData(X, dtype=np.int32, var=var, obs=obs)
    adata.obs = pd.merge(adata.obs, celldata, left_index=True, right_index=True)
    adata.obsm["spatial"] = adata.obs[["center_y", "center_x"]].to_numpy()
    adata.uns["transcripts"] = trxdata
    return adata


def _barcode_names(f: Path) -> pd.DataFrame:
    h = h5py.File(f, "r")
    p = re.compile("B_\\d+")
    barcode_keys = [i for i in h.keys() if p.match(i)]
    names = []
    for b in barcode_keys:
        names = names + [i.decode("UTF-8") for i in h[b]["names"]]
    b = pd.DataFrame({"gene_name": names, "serial": False})
    return b


def _serial_names(f) -> pd.DataFrame:
    h = h5py.File(f, "r")
    names = [i.decode("UTF-8") for i in h["S"]["names"]]
    s = pd.DataFrame({"gene_name": names, "serial": True})
    return s


def add_probe_type(adata: ad.AnnData, f: Path) -> None:
    b = _barcode_names(f)
    s = _serial_names(f)
    adata.var["barcode"] = adata.var.index.isin(b["gene_name"])
    adata.var["serial"] = adata.var.index.isin(s["gene_name"])

