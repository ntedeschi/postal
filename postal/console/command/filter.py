import anndata as ad
from pathlib import Path

from cleo.commands.command import Command
from cleo.helpers import argument, option

from postal import filter
from postal.filter import filter_with_mads, filter_manually
from postal.config import read_config


class FilterCommand(Command):
    name = "filter"
    description = "Filter AnnData cells."

    arguments = [
        argument(
            "config_file",
            description="YAML file with filter parameters.",
            optional=False,
        ),
    ]

    def handle(self):
        pc = read_config(self.argument("config_file"))
        config = pc.filterc
        adata = ad.read(config.adata_in)
        if config.use_mads:
            filter_with_mads(adata, config.mads_factor)
        else:
            filter_manually(
                adata,
                config.min_counts,
                config.max_counts,
                config.min_genes,
                config.max_genes,
            )

        adata.write(config.adata_out)
