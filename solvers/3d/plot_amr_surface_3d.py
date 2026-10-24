#!/usr/bin/env python3
"""Render AMR wall pressure on the exact three-dimensional CRM surface."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/quick-cfd-matplotlib")
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import TwoSlopeNorm
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np


def configure_axis(axis, bounds, view, title):
    axis.set_xlim(*bounds[0])
    axis.set_ylim(*bounds[1])
    axis.set_zlim(*bounds[2])
    axis.set_box_aspect(tuple(high - low for low, high in bounds), zoom=1.42)
    axis.set_proj_type("persp", focal_length=0.9)
    axis.view_init(elev=view[0], azim=view[1])
    axis.set_title(title, pad=4)
    axis.set_axis_off()


def add_surface(axis, vertices, colors, edges=False):
    collection = Poly3DCollection(
        vertices,
        facecolors=colors,
        edgecolors=(0.02, 0.03, 0.04, 0.13) if edges else "none",
        linewidths=0.08 if edges else 0.0,
        rasterized=True,
    )
    axis.add_collection3d(collection)
    return collection


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--geometry", required=True, type=Path)
    parser.add_argument("--run", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--cp-min", type=float, default=-0.8)
    parser.add_argument("--cp-max", type=float, default=1.2)
    args = parser.parse_args()

    with np.load(args.geometry) as source:
        points = np.asarray(source["points"], dtype=np.float64)
        triangles = np.asarray(source["triangles"], dtype=np.int64)
    table = np.genfromtxt(args.run / "surface.csv", delimiter=",", names=True)
    table = np.atleast_1d(table)
    triangle_ids = np.asarray(table["triangle"], dtype=np.int64)
    if not np.array_equal(triangle_ids, np.arange(len(triangles))):
        raise ValueError("surface.csv must contain every source triangle in order")
    cp = np.asarray(table["Cp"], dtype=np.float64)
    result = json.loads((args.run / "result.json").read_text())

    vertices = points[triangles]
    norm = TwoSlopeNorm(vmin=args.cp_min, vcenter=0.0, vmax=args.cp_max)
    cmap = plt.get_cmap("RdBu_r")
    colors = cmap(norm(np.clip(cp, args.cp_min, args.cp_max)))
    padding = np.array((0.18, 0.18, 0.10))
    lower = points.min(axis=0) - padding
    upper = points.max(axis=0) + padding
    bounds = tuple((float(lo), float(hi)) for lo, hi in zip(lower, upper))

    figure = plt.figure(figsize=(18, 7), layout="constrained")
    axes = [figure.add_subplot(1, 3, index + 1, projection="3d")
            for index in range(3)]
    views = (
        ((28, -118), "Upper-front perspective: AMR wall pressure"),
        ((-27, -62), "Lower-front perspective: AMR wall pressure"),
        ((15, 153), "Rear-quarter perspective: AMR wall pressure"),
    )
    for axis, (view, title) in zip(axes, views):
        add_surface(axis, vertices, colors)
        configure_axis(axis, bounds, view, title)

    scalar = ScalarMappable(norm=norm, cmap=cmap)
    scalar.set_array(cp)
    colorbar = figure.colorbar(scalar, ax=axes, shrink=0.80, pad=0.015)
    colorbar.set_label("pressure coefficient  $C_p$")
    status = "converged" if result.get("converged") else "not converged"
    figure.suptitle(
        "NASA CRM full-aircraft 3D surface — exact triangles colored by "
        "area-averaged Cartesian-AMR wall-fragment pressure\n"
        f"{len(triangles):,} triangles | {result.get('cells', 0):,} active "
        f"controls | M={result.get('Mach', float('nan')):g} | {status}",
        fontsize=16,
    )
    figure.text(
        0.5, 0.008,
        "No geometric interpolation: every color is the conservative area average of exact cut-wall fragments mapped to its source triangle.",
        ha="center", fontsize=10,
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.out, dpi=210, facecolor="white")
    plt.close(figure)

    metadata = {
        "figure": args.out.name,
        "geometry": args.geometry.name,
        "run": args.run.name,
        "surface_triangles": int(len(triangles)),
        "active_controls": int(result.get("cells", 0)),
        "Mach": result.get("Mach"),
        "converged": bool(result.get("converged", False)),
        "Cp_range_rendered": [args.cp_min, args.cp_max],
        "rendering": (
            "exact source triangles painted by conservative area-averaged "
            "AMR cut-wall-fragment Cp; no geometric interpolation"
        ),
    }
    args.out.with_suffix(".json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
