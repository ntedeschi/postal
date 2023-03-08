from pathlib import Path

import anndata as ad
import pandas as pd
import pytest
from cleo.application import Application
from cleo.testers.command_tester import CommandTester
from pandas.testing import assert_frame_equal

from postal.console.command.qc import QCCommand
from postal.qc import QC


@pytest.fixture()
def tester() -> CommandTester:
    application = Application()
    application.add(QCCommand())
    command = application.find("qc")
    return CommandTester(command)

@pytest.fixture()
def paths():
    tests = Path("tests")
    data = tests / "data"
    return {"tests": tests, "data": data}


@pytest.mark.skip()
def test_qc(tester, paths):
    paths = paths
    path = paths['data']
    config_file = path / "config.yaml"
    tester.execute(args=f"{config_file}")
