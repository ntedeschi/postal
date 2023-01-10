import anndata as ad
from pathlib import Path

from cleo.commands.command import Command
from cleo.helpers import argument, option

from postal.qc import QC
from postal.config import read_config


class QCCommand(Command):
    name = "qc"
    description = "Calculate standard QC metrics."

    arguments = [
        argument(
            "config_file",
            description="YAML file with qc module parameters.",
            optional=False,
        ),
    ]

    def handle(self):
        pc = read_config(self.argument("config_file"))
        config = pc.qc
        adata = ad.read(config.adata_in)
        outs = config.outs
        experiment = config.experiment
        refgene = config.refgene
        qc = QC(adata, experiment, outs, refgene)
        qc.cell_area_distribution()
        qc.calculate_qc_metrics()
        summary_stats = qc.summary_stats()
        qc.gene_stats()
        qc.entities_per_cell("detected_genes", "Genes per cell")
        qc.entities_per_cell("detected_transcripts", "Transcripts per cell")
        qc.cell_area_vs_detected_genes()
#        qc.segmentation_stats()
        
        qc.adata.write(outs / f"{experiment}_qc.h5ad")
