from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional
import yaml


@dataclass
class MKADConfig:
    """Parameters for MKAD module"""

    counts_file: Path
    cell_data_file: Path
    transcript_data_file: Path
    decode_file: Path
    outs: Path
    anndata_file: str


def read_mkad(d: Dict):
    a = d["arguments"]
    mkad = MKADConfig(
        Path(a["counts_file"]),
        Path(a["cell_data_file"]),
        Path(a["transcript_data_file"]),
        Path(a["decode_file"]),
        Path(a["outs"]),
        a["anndata_file"],
    )
    return mkad


@dataclass
class QCConfig:
    adata_in: Path
    experiment: str
    refgene: str
    outs: Path


def read_qc(d: Dict):
    a = d["arguments"]
    qc = QCConfig(
        Path(a["adata_in"]),
        a["experiment"],
        a["refgene"],
        Path(a["outs"]),
    )
    return qc



@dataclass
class FilterConfig:
    pass


@dataclass
class LatentConfig:
    method: str
    adata_in: Path
    adata_out: Path
    n_layers: int
    n_latent: int
    model_path: Path
    obsm_col: str
    use_gpu: bool | str


def read_latent(d: Dict):
    a = d["arguments"]
    latent = LatentConfig(
        a["method"],
        Path(a["adata_in"]),
        Path(a["adata_out"]),
        a["n_layers"],
        a["n_latent"],
        a["model_path"],
        a["obsm_col"],
        a["use_gpu"],
    )
    return latent


@dataclass
class ClusterConfig:
    method: str
    adata_in: Path
    adata_out: Path
    neighborhood: int
    resolution: float
    use_rep: str
    key_added: str



def read_cluster(d: Dict):
    a = d["arguments"]
    if a["method"] == "leiden":
        cluster = ClusterConfig(
            a["method"],
            Path(a["adata_in"]),
            Path(a["adata_out"]),
            a["neighborhood"],
            a["resolution"],
            a["use_rep"],
            a["key_added"],
        )
    return cluster


@dataclass
class NormalizationConfig:
    pass


@dataclass
class SCVIConfig:
    pass


@dataclass
class PostalConfig:
    mkad: Optional[MKADConfig] = None
    qc: Optional[QCConfig] = None
    filtering: Optional[FilterConfig] = None
    latent: Optional[LatentConfig] = None
    cluster: Optional[ClusterConfig] = None
    normalization: Optional[NormalizationConfig] = None
    scvi: Optional[SCVIConfig] = None


def read_config(path: Path):
    with open(path, mode="rt", encoding="utf-8") as f:
        ds = yaml.safe_load(f)

    pc = PostalConfig()

    for d in ds:
        m = d["module"]
        if m == "mkad":
            pc.mkad = read_mkad(d)
        if m == "qc":
            pc.qc = read_qc(d)
        if m == "latent":
            pc.latent = read_latent(d)
        if m == "cluster":
            pc.cluster = read_cluster(d)

    return pc
