from pathlib import Path
import pytest

from postal.config import read_config


@pytest.fixture()
def linux_paths():
    base = Path("/workspace/postal")
    tests = base / "tests"
    data = tests / "data"
    return {"base": base, "tests": tests, "data": data}


def test_read_mkad_config(linux_paths):
    path = linux_paths["data"]
    file = path / "config.yaml"

    pc = read_config("mkad", file)
    mkad = pc.mkad
    assert mkad.counts_file == path / "counts.csv"
    assert mkad.cell_data_file == path / "cell_data.csv"
    assert mkad.transcript_data_file == path / "transcript_data.csv"
    assert mkad.decode_file == path / "decode.hdf5"
    assert mkad.outs == path / "outs"
    assert mkad.anndata_file == "test.h5ad"
