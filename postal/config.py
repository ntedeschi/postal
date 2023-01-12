from dataclasses import dataclass
import numpy as np
from pathlib import Path
from typing import Dict, Optional
import yaml  # type: ignore


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

    if not "counts_file" in a:
        raise ValueError("mkad: missing 'counts_file' argument.  Check yaml file.")

    if not "cell_data_file" in a:
        raise ValueError("mkad: missing 'cell_data_file' argument.  Check yaml file.")

    if not "transcript_data_file" in a:
        raise ValueError("mkad: missing 'transcript_data_file' argument.  Check yaml file.")

    if not "decode_file" in a:
        raise ValueError("mkad: missing 'decode_file' argument.  Check yaml file.")

    if not "outs" in a:
        raise ValueError("mkad: missing 'outs' argument.  Check yaml file.")

    counts_file = Path(a["counts_file"])
    if not counts_file.exists():
        raise FileNotFoundError(f"{counts_file} does not exist")

    cell_data_file = Path(a["cell_data_file"])
    if not cell_data_file.exists():
        raise FileNotFoundError(f"{cell_data_file} does not exist")

    transcript_data_file = Path(a["transcript_data_file"])
    if not transcript_data_file.exists():
        raise FileNotFoundError(f"{transcript_data_file} does not exist")

    decode_file = Path(a["decode_file"])
    if not decode_file.exists():
        raise FileNotFoundError(f"{decode_file} does not exist")

    outs = Path(a["outs"])
    if not outs.exists():
        outs.mkdir(parents=True)

    mkad = MKADConfig(
        counts_file,
        cell_data_file,
        transcript_data_file,
        decode_file,
        outs,
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
    if not "experiment" in a:
        raise ValueError("qc: missing 'experiment' argument.  Check yaml file.")

    if not "adata_in" in a:
        raise ValueError("qc: missing 'adata_in' argument.  Check yaml file.")

    if not "outs" in a:
        raise ValueError("qc: missing 'outs' argument.  Check yaml file.")

    if not "refgene" in a:
        raise ValueError("qc: missing 'refgene' argument.  Check yaml file.")

    adata_in = Path(a["adata_in"])
    if not adata_in.exists():
        raise FileNotFoundError(f"qc: 'adata_in = {adata_in} does not exist")

    outs = Path(a["outs"])
    if not outs.exists():
        outs.mkdir()

    qc = QCConfig(
        adata_in,
        a["experiment"],
        a["refgene"],
        outs,
    )

    return qc


@dataclass
class FilterConfig:
    adata_in: Path
    adata_out: Path
    use_mads: bool = True
    mads_factor: Optional[np.float64] = None
    min_counts: Optional[int] = None
    max_counts: Optional[int] = None
    min_genes: Optional[int] = None
    max_genes: Optional[int] = None


def read_filter(d: Dict):
    a = d["arguments"]
    use_mads = a["use_mads"]
    if use_mads:
        mads_factor = a["mads_factor"]
        min_counts = None
        max_counts = None
        min_genes = None
        max_genes = None
    else:
        mads_factor = None
        min_counts = a["min_counts"]
        max_counts = a["max_counts"]
        min_genes = a["min_genes"]
        max_genes = a["max_genes"]

    ft = FilterConfig(
        Path(a["adata_in"]),
        Path(a["adata_out"]),
        use_mads,
        mads_factor,
        min_counts,
        max_counts,
        min_genes,
        max_genes,
    )

    return ft


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
    filterc: Optional[FilterConfig] = None
    latent: Optional[LatentConfig] = None
    cluster: Optional[ClusterConfig] = None
    normalization: Optional[NormalizationConfig] = None
    scvi: Optional[SCVIConfig] = None


def read_config(module, path: Path):
    with open(path, mode="rt", encoding="utf-8") as f:
        ds = yaml.safe_load(f)

    for d in ds:
        if d["module"] == module:
            arguments = d
            break

    pc = PostalConfig()

    if module == "mkad":
        pc.mkad = read_mkad(arguments)
    if module == "qc":
        pc.qc = read_qc(arguments)
    if module == "filter":
        pc.filterc = read_filter(arguments)
    if module == "latent":
        pc.latent = read_latent(arguments)
    if module == "cluster":
        pc.cluster = read_cluster(arguments)

    return pc
