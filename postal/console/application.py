from cleo.application import Application

from postal.console.command.mkad import MKADCommand
from postal.console.command.qc import QCCommand
from postal.console.command.filter import FilterCommand
from postal.console.command.cluster import ClusterCommand
from postal.console.command.latent import LatentCommand

application = Application()
application.add(MKADCommand())
application.add(QCCommand())
application.add(FilterCommand())
application.add(ClusterCommand())
application.add(LatentCommand())


def run() -> None:
    application.run()
