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
def linux_paths():
    base = Path("/workspace/postal")
    tests = base / "tests"
    data = tests / "data"
    return {"base": base, "tests": tests, "data": data}


#@pytest.mark.skip()
def test_qc(tester, linux_paths):
    paths = linux_paths
    path = paths['data']
    config_file = path / "config.yaml"
    tester.execute(args=f"{config_file}")
