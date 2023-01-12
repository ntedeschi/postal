from pathlib import Path

from cleo.commands.command import Command
from cleo.helpers import argument, option

from postal.mkad import (
    add_probe_type,
    mkad,
    read_celldata,
    read_counts,
    read_transcripts,
)
from postal.config import read_config


class MKADCommand(Command):
    name = "mkad"
    description = "Create an AnnData object from a count matrix and metadata."

    arguments = [
        argument(
            "config_file",
            description="YAML file with mkad parameters.",
            optional=False,
        ),
    ]

    def handle(self):
        pc = read_config("mkad", self.argument("config_file"))
        config = pc.mkad
        counts = read_counts(config.counts_file)
        celldata = read_celldata(config.cell_data_file)
        transcripts = read_transcripts(config.transcript_data_file)
        adata = mkad(counts, celldata, transcripts)
        add_probe_type(adata, config.decode_file)

        outs = config.outs
        adata.write(outs / config.anndata_file)

