from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from pandas.testing import assert_frame_equal
import pytest
import scanpy as sc
from scipy import stats

from postal.filter import (
    filter_with_mads,
    _filter_cells_by_detected_genes,
    _filter_cells_by_transcript_counts,
    filter_manually,
)


@pytest.fixture()
def paths():
    tests = Path("tests")
    data = tests / "data"
    return {"tests": tests, "data": data}


@pytest.fixture()
def get_adata(paths):
    adata = ad.read(paths["data"] / "test_after_qc.h5ad")
    return adata


def test_mads():
    x = np.array(range(100))

    mads_factor = np.float64(0.5)
    med = np.median(x)
    mad = stats.median_abs_deviation(x)
    max_x = med + mads_factor * mad
    min_x = med - mads_factor * mad
    min_x = max(np.float64(1), min_x)

    min_x = min_x.astype(np.int64)
    max_x = max_x.astype(np.int64)

    assert (min_x, max_x) == (37, 62)

    mads_factor = np.float64(1000)
    med = np.median(x)
    mad = stats.median_abs_deviation(x)
    max_x = med + mads_factor * mad
    min_x = med - mads_factor * mad
    min_x = max(np.float64(1), min_x)

    min_x = min_x.astype(np.int64)
    max_x = max_x.astype(np.int64)

    assert (min_x, max_x) == (1, 25049)


def test_filter_cells_by_transcript_counts(get_adata):
    adata = get_adata.copy()
    assert adata.shape == (40954, 715)

    _filter_cells_by_transcript_counts(adata, min_counts=100, max_counts=500)
    sc.pp.filter_cells(adata, max_counts=500)
    sc.pp.filter_cells(adata, min_counts=100)

    assert adata.shape == (15823, 715)


def test_filter_cells_by_detected_genes(get_adata):
    adata = get_adata
    assert adata.shape == (40954, 715)

    _filter_cells_by_detected_genes(adata, min_genes=20, max_genes=100)
    #    sc.pp.filter_cells(adata, max_genes=100)
    #    sc.pp.filter_cells(adata, min_genes=20)

    assert adata.shape == (14426, 715)


def test_filter_with_mads(get_adata):
    adata = get_adata
    assert adata.shape == (40954, 715)

    mads_factor = 3
    filter_with_mads(adata, mads_factor)

    assert adata.shape == (38943, 715)


def test_filter_manually(get_adata):
    adata = get_adata
    assert adata.shape == (40954, 715)

    filter_manually(
        adata,
        min_counts=100,
        max_counts=500,
        min_genes=20,
        max_genes=100,
    )

    assert adata.shape == (11698, 715)
