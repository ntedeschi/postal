import anndata as ad
from pathlib import Path
import scvi


def method_scvi(
    adata: ad.AnnData,
    n_layers: int,
    n_latent: int,
    model_path: Path,
    obsm_col: str,
    use_gpu: bool | str,
):
    scvi.model.SCVI.setup_anndata(
        adata,
        layer=None,
    )

    model = scvi.model.SCVI(
        adata,
        n_layers=n_layers,
        n_latent=n_latent,
    )

    model.train(use_gpu=use_gpu)
    model.save(model_path)
    latent = model.get_latent_representation()
    adata.obsm[obsm_col] = latent
    adata.layers["scvi_normalized"] = model.get_normalized_expression(library_size=1e4)
