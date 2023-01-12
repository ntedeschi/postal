import anndata as ad
from pathlib import Path

from cleo.commands.command import Command
from cleo.helpers import argument, option

from postal import cluster
from postal.config import read_config


class ClusterCommand(Command):
    name = "cluster"
    description = "Cluster data."

    arguments = [
        argument(
            "config_file",
            description="YAML file with cluster parameters.",
            optional=False,
        ),
    ]

    def handle(self):
        pc = read_config("cluster", self.argument("config_file"))
        config = pc.cluster
        adata = ad.read(config.adata_in)
        if config.method == "leiden":
            cluster.leiden(
                adata,
                config.neighborhood,
                config.resolution,
                config.use_rep,
                config.key_added,
            )
            adata.write(config.adata_out)
