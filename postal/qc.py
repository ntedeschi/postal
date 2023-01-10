import anndata as ad
import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch
from matplotlib.cbook import boxplot_stats
import numpy as np
import pandas as pd
from pathlib import Path
import scanpy as sc
from scipy.stats import gaussian_kde
import seaborn as sns
from typing import List, Optional


class QC:
    """
    Quality control calculations on seqFish data.

    Parameters
    ----------
    adata
        AnnData object created by `mkad` module.
    experiment
        Name of experiment.
    refgene
        Gene used as high expressor reference.
    """

    def __init__(self, adata: ad.AnnData, experiment: str, outs: Path, refgene: str = "EEF2"):
        self.adata = adata.copy()
        self.experiment = experiment
        self.experiment = self.experiment.replace(" ", "_")
        self.experiment = self.experiment.replace("-", "_")
        self.outs = outs
        self.figdir = self.outs / "figures"
        self.refgene = refgene

        if not self.outs.exists():
            self.outs.mkdir()
        if not self.figdir.exists():
            self.figdir.mkdir()

        self._remove_unused()
        self.barcode = self.adata[:, self.adata.var["barcode"]].copy()
        self.serial = self.adata[:, self.adata.var["serial"]].copy()

    def cell_area_distribution(self):
        """
        Calculate histogram and spatial plot of cell area.
        """

        cell_metadata = self.adata.obs.copy()
        cp = sns.color_palette("colorblind")

        q = cell_metadata["area"].quantile([0.025, 0.975])
        cell_metadata["area_extremes"] = "_"
        cell_metadata.loc[
            cell_metadata["area"] < q.values[0], "area_extremes"
        ] = f"< {q.index[0]}"

        cell_metadata.loc[
            cell_metadata["area"] >= q.values[1], "area_extremes"
        ] = f">= {q.index[1]}"

        fig, axs = plt.subplots(
            1, 3, figsize=(15, 5), gridspec_kw={"wspace": 0.3, "hspace": 0.6}
        )
        sns.histplot(
            x=cell_metadata["area"], ax=axs[0], hue=cell_metadata["area_extremes"]
        ).set(xlabel="Area (pixels)", title="Histogram")

        sns.scatterplot(
            x=cell_metadata["center_x_clockwise90"],
            y=cell_metadata["center_y_clockwise90"],
            ax=axs[1],
            s=10,
            hue=cell_metadata["area_extremes"],
        )

        cell_metadata.loc[cell_metadata["area_extremes"] == "_", "area_extremes"] = None

        sns.scatterplot(
            x=cell_metadata["center_x_clockwise90"],
            y=cell_metadata["center_y_clockwise90"],
            ax=axs[2],
            s=10,
            hue=cell_metadata["area_extremes"],
            palette=[cp[i] for i in [0, 2]],
        )

        for i in [0, 1, 2]:
            axs[i].get_legend().remove()

        axs[1].title.set_text("Spatial")
        axs[2].title.set_text("Outliers")
        fig.subplots_adjust(top=0.85)
        extremes = (
            f"Lower 2.5% = {int(q.values[0])}    -   Upper 97.5% = {int(q.values[1])}"
        )

        fig.text(0.13, 0, extremes)
        fig.suptitle(f"{self.experiment} — Cell Area")
        figfile = self.figdir / f"{self.experiment}_cellarea.png"
        fig.savefig(figfile)

    def _remove_unused(self):
        used = [name for name in self.adata.var_names if not name.startswith("UNUSED")]
        self.adata = self.adata[:, used].copy()

    def calculate_qc_metrics(self):
        sc.pp.calculate_qc_metrics(
            self.adata, percent_top=None, inplace=True, log1p=False
        )
        sc.pp.calculate_qc_metrics(
            self.barcode, percent_top=None, inplace=True, log1p=False
        )
        sc.pp.calculate_qc_metrics(
            self.serial, percent_top=None, inplace=True, log1p=False
        )

        def _rename(adata):
            adata.obs = adata.obs.rename(
                columns={
                    "n_genes_by_counts": "detected_genes",
                    "total_counts": "detected_transcripts",
                }
            )
            adata.var = adata.var.rename(
                columns={
                    "n_cells_by_counts": "cells",
                    "pct_dropout_by_counts": "pct_cells_absent",
                }
            )

        def _trx_per_gene(adata):
            adata.obs["transcripts_per_gene"] = (
                adata.obs["detected_transcripts"] / adata.obs["detected_genes"]
            )

        def _mean_counts_wp(adata):
            adata.var["mean_counts_wp"] = adata.var["total_counts"] / adata.var["cells"]

        _rename(self.adata)
        _rename(self.barcode)
        _rename(self.serial)

        _trx_per_gene(self.adata)
        _trx_per_gene(self.barcode)
        _trx_per_gene(self.serial)

        _mean_counts_wp(self.adata)
        _mean_counts_wp(self.barcode)
        _mean_counts_wp(self.serial)

    def summary_stats(self) -> pd.DataFrame:
        """Return dataframe summarizing gene and transcript statistics."""

        def _refgene_absence(adata):
            if self.refgene in adata.var_names:
                return adata.var.loc[self.refgene, "pct_cells_absent"]
            else:
                return np.nan

        def _summary(adata):
            n_cells = adata.shape[0]
            n_genes = adata.shape[1]
            mean_genes = np.mean(adata.obs["detected_genes"])
            mean_trx = np.mean(adata.obs["detected_transcripts"])
            mean_trx_per_gene = np.mean(adata.obs["transcripts_per_gene"])

            df = pd.DataFrame(
                {
                    "cells": n_cells,
                    "genes": n_genes,
                    "genes/cell": mean_genes,
                    "transcripts/cell": mean_trx,
                    "transcripts/gene/cell": mean_trx_per_gene,
                    f"{self.refgene} pct absent": _refgene_absence(adata),
                },
                index=["a"],
            )
            return df

        df = pd.concat(
            [_summary(self.adata), _summary(self.barcode), _summary(self.serial)]
        )
        df.index = ["All genes", "Barcode genes", "Serial genes"]
        df = df.astype(
            {
                "cells": int,
                "genes": int,
                "transcripts/cell": np.float64,
                "transcripts/gene/cell": np.float64,
                f"{self.refgene} pct absent": np.float64,
            }
        )
        df.style.format(precision=1)

        return df

    def gene_stats(self):
        """Calculate and save figure of gene level statistics."""

        fig = plt.figure(constrained_layout=True, figsize=(12, 10))
        subfigs = fig.subfigures(3, 3, wspace=0.07, hspace=0.07)

        col_titles = ["All genes", "Barcode genes", "Serial genes"]
        row_vars = ["cells", "total_counts", "mean_counts_wp"]
        row_titles = [
            "Cells found in",
            "Total transcript counts",
            "Mean counts/cell\n when present",
        ]

        flierprops = dict(
            marker=".", markerfacecolor="indianred", markeredgecolor="none"
        )

        adatas = [self.adata, self.barcode, self.serial]

        for col in [0, 1, 2]:
            for row in [0, 1, 2]:
                rv = row_vars[row]  # Variable plotted this row (row variable)
                sf = subfigs[row][col]
                ax0, ax1 = sf.subplots(
                    2, 1, gridspec_kw={"height_ratios": [1, 3]}, sharex=False
                )
                ax0.boxplot(
                    x=adatas[col].var[rv], notch=0, flierprops=flierprops, vert=0
                )
                bx = boxplot_stats(adatas[col].var[rv]).pop(0)
                ax0.annotate(format(bx["med"], "g"), [bx["med"], 1.15])
                x = adatas[col].var[rv]
                ax1.hist(
                    x=x[np.logical_and(x <= bx["whishi"], x >= bx["whislo"])],
                    bins="auto",
                )
                # Connecting lines
                bx_wl = bx["whislo"] / 1e6 if bx["whislo"] > 1e6 else bx["whislo"]
                bx_wh = bx["whishi"] / 1e6 if bx["whishi"] > 1e6 else bx["whishi"]

                conn1 = ConnectionPatch(
                    xyA=(bx_wl, 1),
                    coordsA="data",
                    axesA=ax0,
                    xyB=(0, 1),
                    coordsB="axes fraction",
                    axesB=ax1,
                    color="gray",
                    ls=":",
                )
                sf.add_artist(conn1)

                conn2 = ConnectionPatch(
                    xyA=(bx_wh, 1),
                    coordsA="data",
                    axesA=ax0,
                    xyB=(1, 1),
                    coordsB="axes fraction",
                    axesB=ax1,
                    color="gray",
                    ls=":",
                )
                sf.add_artist(conn2)

                # Annotate mean with vertical line
                ax1.axvline(x=bx["med"], ls="--", c="orange")
                ax1.axvline(x=bx["mean"], ls="--", c="purple")

                # Annotate genes
                ann_genes = []
                ann_genes.append(self.refgene) if any(
                    adatas[col].var.index == self.refgene
                ) else None
                ann_genes.append(
                    adatas[col]
                    .var[adatas[col].var[rv] == max(adatas[col].var[rv])]
                    .index[0]
                )

                for ag in ann_genes:
                    ag_rv = adatas[col].var.loc[ag, rv]
                    if np.logical_and(ag_rv <= bx["whishi"], ag_rv >= bx["whislo"]):
                        ax1.annotate(
                            ag,
                            [ag_rv, 0],
                            fontsize=6,
                            horizontalalignment="right",
                            verticalalignment="center",
                        )
                    else:
                        ax0.annotate(
                            ag,
                            [ag_rv, 1],
                            xytext=[ag_rv, 0.6],
                            arrowprops=dict(
                                facecolor="black", headwidth=0, width=1, shrink=0.5
                            ),
                            fontsize=6,
                            horizontalalignment="right",
                            verticalalignment="center",
                        )

                    ax0.set_title(col_titles[col], fontsize=14) if row == 0 else ""
                    ax1.set_xlabel(row_titles[row], fontsize=12)
                    ax1.set_ylabel(
                        "Genes", fontsize=12
                    ) if col == 0 else ax1.set_ylabel("")

            fig.suptitle(f"{self.experiment} Gene Stats")
            fig.savefig(self.figdir / f"{self.experiment}_genestats.png")

    def entities_per_cell(
        self,
        row_var: str,
        row_title: str,
        col_titles: List[str] = ["All genes", "Barcode genes", "Serial genes"],
    ):
        """
        Plot histograms and spatial plots of an entity per cell

        An entity can be genes or transcripts or some other thing that
        is measured in a cell.
        """

        fig = plt.figure(constrained_layout=True, figsize=(12, 6))
        fig.suptitle(f"{self.experiment} — {row_title}")

        subfigs = fig.subfigures(1, 3, wspace=0.07, hspace=0.2)

        flierprops = dict(
            marker=".", markerfacecolor="indianred", markeredgecolor="none"
        )

        adatas = [self.adata, self.barcode, self.serial]
        for col in [0, 1, 2]:
            # Remove any nan from the variable (nan messes up boxplot_stats)
            ad_col = adatas[col][~np.isnan(adatas[col].obs[row_var])]

            ax = subfigs[col].subplots(
                3, 1, gridspec_kw={"height_ratios": [1, 2, 6]}, sharex=False
            )

            ax[0].boxplot(x=ad_col.obs[row_var], notch=0, flierprops=flierprops, vert=0)
            bx = boxplot_stats(ad_col.obs[row_var][~np.isnan(ad_col.obs[row_var])]).pop(
                0
            )
            ax[0].annotate(format(bx["med"], "g"), [bx["med"], 1.15])
            x = ad_col.obs[row_var]
            ax[1].hist(
                x=x[np.logical_and(x <= bx["whishi"], x >= bx["whislo"])],
                bins="auto",
            )
            obs = ad_col.obs.sort_values(row_var)
            sns.scatterplot(
                x=obs["center_x_clockwise90"],
                y=obs["center_y_clockwise90"],
                ax=ax[2],
                s=10,
                hue=obs[row_var],
            )
            ax[2].get_legend().remove()
            conn1 = ConnectionPatch(
                xyA=(bx["whislo"], 1),
                coordsA="data",
                axesA=ax[0],
                xyB=(0, 1),
                coordsB="axes fraction",
                axesB=ax[1],
                color="gray",
                ls=":",
            )
            subfigs[col].add_artist(conn1)

            conn2 = ConnectionPatch(
                xyA=(bx["whishi"], 1),
                coordsA="data",
                axesA=ax[0],
                xyB=(1, 1),
                coordsB="axes fraction",
                axesB=ax[1],
                color="gray",
                ls=":",
            )
            subfigs[col].add_artist(conn2)

            # Annotate mean with vertical line
            ax[1].axvline(x=bx["med"], ls="--", c="orange")
            ax[1].axvline(x=bx["mean"], ls="--", c="purple")

            ax[0].set_title(col_titles[col], fontsize=14)
            ax[1].set_xlabel(row_title, fontsize=12)
            ax[1].set_ylabel("Cells", fontsize=12) if col == 0 else ax[1].set_ylabel("")
            ax[2].set_ylabel("y", fontsize=14) if col == 0 else ax[2].set_ylabel("")
            ax[2].set_xlabel("x", fontsize=14)

        fig.savefig(self.figdir / f"{self.experiment}_cell_{row_var}.png")

    def cell_area_vs_detected_genes(self):
        fig, axs = plt.subplots(
            1, 3, figsize=(12, 5), gridspec_kw={"wspace": 0.5, "hspace": 1}
        )

        titles = ["All genes", "Barcode genes", "Serial genes"]
        adatas = [self.adata, self.barcode, self.serial]

        for i in [0, 1, 2]:
            x, y = adatas[i].obs["area"], adatas[i].obs["detected_genes"]
            xy = np.vstack([x, y])
            z = gaussian_kde(xy)(xy)

            axs[i].scatter(x, y, c=z, alpha=0.5, cmap="viridis")
            axs[i].set_title(titles[i])
            axs[i].set_ylabel("Detected genes")
            axs[i].set_xlabel("Cell area (pixels)")

        fig.subplots_adjust(top=0.85)
        fig.suptitle(f"{self.experiment} — Cell area vs Genes")
        fig.savefig(self.figdir / f"{self.experiment}_area_vs_genes.png")

    def segmentation_stats(self) -> Optional[pd.DataFrame]:
        """
        Calculates the number of transcripts not in cells.

        Requires self.adata.uns["transcipts"] DataFrame to exist
        and have a `cell` column the cell number in which the transcript
        is detected.
        """

        if "transcripts" not in self.adata.uns:
            return None
        transcripts = self.adata.uns["transcripts"]
        total_transcripts = int(transcripts["cell"].count())
        not_in_cell = transcripts.loc[transcripts["cell"] == 0]
        not_in_cell = int(not_in_cell["cell"].count())
        percent = (not_in_cell / total_transcripts) * 100
        df = pd.DataFrame(
            {
                "number of transcripts": total_transcripts,
                "number not in a cell": not_in_cell,
                "percent not in cell": percent,
            },
            index=pd.Series(range(1)),
        )

        return df

