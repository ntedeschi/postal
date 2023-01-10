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
def linux_paths():
    base = Path("/workspace/postal")
    tests = base / "tests"
    data = tests / "data"
    return {"base": base, "tests": tests, "data": data}


def test_use_mads(tester, linux_paths):
    path = linux_paths['data']
    config_file = path / "config.yaml"
    tester.execute(args=f"{config_file}")
    
def test_filter_manually(tester, linux_paths):
    path = linux_paths['data']
    config_file = path / "manual_filter.yaml"
    tester.execute(args=f"{config_file}")
