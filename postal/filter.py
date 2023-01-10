import anndata as ad  # type: ignore
import numpy as np
import scanpy as sc  # type: ignore
from scipy import stats  # type: ignore
from typing import Tuple


def _mads(x: np.ndarray, mads_factor: np.float64) -> Tuple[np.int64, np.int64]:
    """Calculate the min and max of x using median absolution deviation."""

    med = np.median(x)
    mad = stats.median_abs_deviation(x)
    max_x = med + mads_factor * mad
    min_x = med - mads_factor * mad
    min_x = max(np.float64(1), min_x)

    min_x = min_x.astype(np.int64)
    max_x = max_x.astype(np.int64)

    return min_x, max_x


def _filter_cells_by_transcript_counts(
    adata: ad.AnnData, min_counts: np.int64, max_counts: np.int64
):
    sc.pp.filter_cells(adata, min_counts=min_counts)
    sc.pp.filter_cells(adata, max_counts=max_counts)


def _filter_cells_by_detected_genes(
    adata: ad.AnnData, min_genes: np.int64, max_genes: np.int64
):
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_cells(adata, max_genes=max_genes)


def filter_with_mads(adata: ad.AnnData, mads_factor: np.float32):
    counts = adata.obs.detected_transcripts.to_numpy()
    min_x, max_x = _mads(counts, mads_factor)
    _filter_cells_by_transcript_counts(adata, min_x, max_x)

    genes = adata.obs.detected_genes.to_numpy()
    min_x, max_x = _mads(genes, mads_factor)
    _filter_cells_by_detected_genes(adata, min_x, max_x)


def filter_manually(
    adata: ad.AnnData,
    min_counts: np.int64,
    max_counts: np.int64,
    min_genes: np.int64,
    max_genes: np.int64,
):
    counts = adata.obs.detected_transcripts.to_numpy()
    _filter_cells_by_transcript_counts(adata, min_counts, max_counts)

    genes = adata.obs.detected_genes.to_numpy()
    _filter_cells_by_detected_genes(adata, min_genes, max_genes)
