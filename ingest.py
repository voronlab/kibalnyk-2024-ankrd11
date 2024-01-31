from functools import partial
from math import atan2, cos, pi, sin
from os import makedirs
from os.path import isdir
from pathlib import Path
from time import time
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
from anndata import AnnData, concat
from matplotlib.patches import Polygon


def scp(*args, **kwargs):
    """
    scp: spatial cell plot

    wrapper around sq.pl.spatial_scatter to make sure that axes are not inverted
    """
    sq.pl.spatial_scatter(*args, **kwargs)
    x_lim = plt.xlim()
    y_lim = plt.ylim()
    plt.xlim(min(x_lim), max(x_lim))
    plt.ylim(min(y_lim), max(y_lim))


def vp(adata: AnnData):
    """
    vp: vizgen pipeline

    basic single cell data pipeline
    mutates given object in-place with additional data
    """
    sc.pp.calculate_qc_metrics(adata, percent_top=(50, 100))
    sc.pp.filter_cells(adata, min_counts=10)
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata)
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.tl.leiden(adata)


def vr(path: str, tags: dict[str, Any]) -> AnnData:
    """
    vr: vizgen read

    reads a vizgen output folder
    tags dict allows for additional of sample specific metadata (sample id, condition, etc.)
    """
    if isdir(Path(path) / "Cellpose"):
        adata = sq.read.vizgen(
            path,
            counts_file="Cellpose/cellpose_cell_by_gene.csv",
            meta_file="Cellpose/cellpose_cell_metadata.csv",
            transformation_file="micron_to_mosaic_pixel_transform.csv",
        )
    else:
        adata = sq.read.vizgen(
            path,
            counts_file="cell_by_gene.csv",
            meta_file="cell_metadata.csv",
            transformation_file="micron_to_mosaic_pixel_transform.csv",
        )

    for k, v in tags.items():
        adata.obs[k] = v

    return adata


def within_poly(polygon: list[tuple[float, float]], point: tuple[float, float]) -> bool:
    """
    checks if a point is within a polygon
    adapted from: https://web.archive.org/web/20130126163405/http://geomalgorithms.com/a03-_inclusion.html
    """
    winding_num: int = 0
    x: float = point[0]
    y: float = point[1]
    for (cur_x, cur_y), (prv_x, prv_y) in zip(
        polygon, (polygon[i] for i in range(-1, len(polygon) - 1))
    ):
        if cur_y <= y:
            if (prv_y > y) and (
                (prv_x - cur_x) * (y - cur_y) > (x - cur_x) * (prv_y - cur_y)
            ):
                winding_num += 1
        elif (prv_y <= y) and (
            (prv_x - cur_x) * (y - cur_y) < (x - cur_x) * (prv_y - cur_y)
        ):
            winding_num -= 1

    return winding_num != 0


def rot_mat(theta: float) -> np.ndarray:
    """
    creates a rotation matrix using the angle theta
    """
    return np.array([[cos(theta), -sin(theta)], [sin(theta), cos(theta)]])


def scr(adata: AnnData, col: str, show_after: bool = False):
    """
    scr: spatial cell rotation

    creates a spatial plot of cells in sample
    expects two mouse clicks before being closed:
    - first on point that will be rotated to the top of the image
    - second on point that will be rotated to the bottom of the image
    after the first plot is closed, will create a second plot
    """
    coords = []
    scp(adata, color=col, shape=None)
    plt.figure(1).canvas.mpl_connect(
        "button_press_event", lambda e: coords.append((e.xdata, e.ydata))
    )
    plt.show()

    (nose_x, nose_y), (tail_x, tail_y) = coords
    shift = np.array([tail_x + (nose_x - tail_x) / 2, tail_y + (nose_y - tail_y) / 2])
    r = rot_mat(((pi / 2) - atan2(nose_y - shift[1], nose_x - shift[0])))
    adata.obsm["spatial"] = np.row_stack(
        [np.matmul(r, p - shift) for p in adata.obsm["spatial"]]
    )

    if show_after:
        scp(adata, color=col, shape=None)
        plt.show()


def scs(
    adata: AnnData,
    uns_slot_name: str,
    col: str,
    num_boxes: int = 1,
    show_after: bool = False,
    only_show: None | str = None,
):
    """
    scs: spatial cell subset

    first creates a spatial plot of cells in sample
    expects multiple mouse before before being closed
    - ***ALL*** mouse clicks will be used to define a bounding box around the heart
    after the second plot is closed, will create two more plots
    - one will be the entire sample with a bounding box drawn around the heart
    - the other will be of just the subsetted cells
    """
    bb_list: list[list[tuple[float, float]]] = []
    for _ in range(0, num_boxes):
        bb = []
        scp(
            adata if only_show is None else adata[adata.obs[only_show]],
            color=col,
            shape=None,
        )
        for bb in bb_list:
            plt.gca().add_patch(
                Polygon(bb, linewidth=2, edgecolor="r", facecolor="none")
            )
        plt.figure(1).canvas.mpl_connect(
            "button_press_event",
            lambda e: bb.append((e.xdata, e.ydata)),
        )
        plt.show()
        bb_list.append(bb)

    uns_slot: str = uns_slot_name or "bb-{time()}"
    adata.obs[f"within.{uns_slot}"] = list(
        map(
            lambda p: any(within_poly(bb, p) for bb in bb_list),
            adata.obsm["spatial"],
        )
    )
    adata.uns[uns_slot] = bb_list

    if show_after:
        # zoomed out, with outline
        scp(
            adata if only_show is None else adata[adata.obs[only_show]],
            color=col,
            shape=None,
        )
        for bb in adata.uns[uns_slot]:
            plt.gca().add_patch(
                Polygon(bb, linewidth=2, edgecolor="r", facecolor="none")
            )

        # zoomed in
        scp(
            adata[adata.obs[f"within.{uns_slot}"]],
            color=col,
            shape=None,
        )

        plt.show()


def merge(adatas: list[AnnData], prefix_tag: bool = False) -> AnnData:
    """
    wrapper around anndata.concat to easily merge multiple anndata objects

    optional argument allows tagging cell idents in each object with a prefix
    to ensure idents are globally unique
    """
    if prefix_tag:
        for adata in adatas:
            adata.obs_names = [
                f"{adata.obs['source'].iloc[0]}-{name}" for name in adata.obs_names
            ]

    return concat(adatas, join="outer")


if __name__ == "__main__":
    adatas = [
        (
            "202309211049_MsEmbryo-VS147-AH59-2-5-A_VMSC18002/region_1",
            {"source": "AH59-5", "isKO": False},
        ),
        (
            "202309151103_MsEmbryo-VS147-AH59-4-6-C_VMSC17602_reanalysis/region_1",
            {"source": "AH59-6", "isKO": False},
        ),
        (
            "20231204_Ankrd11_OFT/Region1",
            {"source": "AI04-2", "isKO": False},
        ),
        (
            "202309181306_MsEmbryo-VS147-AH50-59-1-A_VMSC07101/region_0",
            {"source": "AH50-1", "isKO": True},
        ),
        (
            "202309151104_MsEmbryo-VS147-AH59-2-5-B_VMSC16102/region_1",
            {"source": "AH59-2", "isKO": True},
        ),
        (
            "20231204_Ankrd11_OFT/Region0",
            {"source": "AI13-3", "isKO": True},
        ),
    ]
    for path, meta in adatas:
        adata = vr(Path("data") / path, meta)

        # apply rotation
        scr(adata, "Postn")
        # create heart bounding box
        scs(adata, "bb.heart", "Postn", show_after=True)

        outdir = Path("out") / path
        makedirs(outdir, exist_ok=True)

        adata.to_df().T.to_csv(outdir / "counts.csv")
        adata.obs[
            [
                "source",
                "isKO",
                "within.bb.heart",
            ]
        ].to_csv(outdir / "meta.csv")
        pd.concat(
            [
                pd.DataFrame(adata.obs_names, columns=["cell"]),
                pd.DataFrame(adata.obsm["spatial"], columns=["x", "y"]),
            ],
            axis=1,
        ).to_csv(outdir / "centroids.csv")
        pd.DataFrame(
            pd.Series(
                [
                    adata.uns["bb.heart"],
                ],
                [
                    "bb.heart",
                ],
            ),
            columns=["coords"],
        ).to_csv(outdir / "bb.csv")
