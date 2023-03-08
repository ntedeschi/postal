from pathlib import Path

import anndata as ad
import pytest
from cleo.application import Application
from cleo.testers.command_tester import CommandTester

from postal.console.command.filter import FilterCommand
import postal.latent

@pytest.fixture()
def tester() -> CommandTester:
    application = Application()
    application.add(FilterCommand())
    command = application.find("filter")
    return CommandTester(command)


@pytest.fixture()
def paths():
    tests = Path("tests")
    data = tests / "data"
    return {"tests": tests, "data": data}


def test_use_mads(tester, paths):
    path = paths['data']
    config_file = path / "config.yaml"
    tester.execute(args=f"{config_file}")
    
def test_filter_manually(tester, paths):
    path = paths['data']
    config_file = path / "manual_filter.yaml"
    tester.execute(args=f"{config_file}")
