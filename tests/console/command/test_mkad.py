from pathlib import Path

import anndata as ad
import pandas as pd
import pytest
from cleo.application import Application
from cleo.testers.command_tester import CommandTester
from pandas.testing import assert_frame_equal

from postal.console.command.mkad import MKADCommand
from postal.mkad import (add_probe_type, read_celldata, read_counts,
                         read_transcripts)


@pytest.fixture()
def tester() -> CommandTester:
    application = Application()
    application.add(MKADCommand())
    command = application.find("mkad")
    return CommandTester(command)

@pytest.fixture()
def paths():
    tests = Path("tests")
    data = tests / "data"
    return {"tests": tests, "data": data}


def test_read_celldata(tester, paths):
    paths = paths
    path = paths['data']
    file = path / "cell_data.csv"
    celldata = read_celldata(file)

    celldata_ref = pd.read_csv(path / "cell_data_ref.csv", index_col=0)
    celldata_ref.index = celldata_ref.index.astype(str)
    celldata_ref.index.name = ''

    assert_frame_equal(celldata_ref, celldata)


def test_read_counts(tester, paths):
    paths = paths
    path = paths['data']
    file = path / "counts.csv"
    counts = read_counts(file)
    counts_ref = pd.read_csv(path / "counts.csv", index_col=0)

    assert_frame_equal(counts_ref, counts)


def test_read_transcripts(tester, paths):
    paths = paths
    path = paths['data']
    file = path / "cell_data.csv"
    file = path / "transcript_data.csv"
    transcripts = read_transcripts(file)
    transcripts_ref = pd.read_csv(path / "transcript_data.csv")
    transcripts_ref = transcripts_ref.rename(columns={"name": "gene_name"})

    assert_frame_equal(transcripts_ref, transcripts)


def test_anndata(tester, paths):
    paths = paths
    path = paths['data']
    counts = path / "counts.csv"
    celldata = path / "cell_data.csv"
    trxdata = path / "transcript_data.csv"
    decode = path / "decode.hdf5"
    output = path / "outs" / "test.h5ad"
    config_file = path / "config.yaml"

    # tester.execute(args=f"{counts} {celldata} {trxdata} {decode} {output}")
    tester.execute(args=f"{config_file}")

    adata_ref = ad.read(path / "test.h5ad")
    adata = ad.read(path / "outs" / "test.h5ad")

    assert_frame_equal(adata_ref.obs, adata.obs)
    output.unlink()

def test_add_probe_type(paths):
    path = paths['data']
    adata = ad.read(path / "test.h5ad")
    hfile = path / "decode.hdf5"
    add_probe_type(adata, hfile)
