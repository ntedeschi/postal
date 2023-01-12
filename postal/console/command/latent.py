import anndata as ad
from pathlib import Path

from cleo.commands.command import Command
from cleo.helpers import argument, option

from postal import latent
from postal.config import read_config


class LatentCommand(Command):
    name = "latent"
    description = "Compute latent representation of data."

    arguments = [
        argument(
            "config_file",
            description="YAML file with latent parameters.",
            optional=False,
        ),
    ]

    def handle(self):
        pc = read_config("latent", self.argument("config_file"))
        config = pc.latent
        adata = ad.read(config.adata_in)
        if config.method == "scvi":
            latent.method_scvi(
                adata,
                config.n_layers,
                config.n_latent,
                config.model_path,
                config.obsm_col,
                config.use_gpu,
            )
            adata.write(config.adata_out)
