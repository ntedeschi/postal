from pathlib import Path

import anndata as ad
import pandas as pd
import pytest
from cleo.application import Application
from cleo.testers.command_tester import CommandTester

from postal.console.command.cluster import ClusterCommand
import postal.cluster

@pytest.fixture()
def tester() -> CommandTester:
    application = Application()
    application.add(ClusterCommand())
    command = application.find("cluster")
    return CommandTester(command)


@pytest.fixture()
def paths():
    tests = Path("tests")
    data = tests / "data"
    return {"tests": tests, "data": data}


def test_leiden(tester, paths):
    path = paths['data']
    config_file = path / "config.yaml"

    tester.execute(args=f"{config_file}")
