#!/usr/bin/env python3
"""Plot a continuous CRM AMR convergence record from restart segments."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def parse_segment(value: str) -> tuple[int, Path]:
    try:
        offset, path = value.split(":", 1)
        return int(offset), Path(path)
    except ValueError as error:
        raise argparse.ArgumentTypeError("segment must be OFFSET:HISTORY.csv") from error


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--segment", action="append", type=parse_segment, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--nasa-cl", type=float, default=0.54239)
    parser.add_argument("--nasa-cd", type=float, default=0.02995)
    parser.add_argument("--nasa-cm", type=float, default=-0.06249)
    args = parser.parse_args()

    records: list[np.ndarray] = []
    segment_metadata = []
    for segment_index, (offset, path) in enumerate(args.segment):
        data = np.atleast_1d(np.genfromtxt(path, delimiter=",", names=True))
        if data.size == 0:
            raise RuntimeError(f"empty history: {path}")
        iteration = offset + data["iteration"]
        if segment_index and iteration[0] == records[-1]["iteration"][-1]:
            data = data[1:]
            iteration = iteration[1:]
        record = np.empty(
            data.size,
            dtype=[
                ("iteration", "f8"),
                ("CD", "f8"),
                ("CL", "f8"),
                ("CM_y", "f8"),
                ("global", "f8"),
                ("near", "f8"),
            ],
        )
        record["iteration"] = iteration
        for field in ("CD", "CL", "CM_y"):
            record[field] = data[field]
        record["global"] = data["res_total_max"]
        record["near"] = data["near_res_total_max"]
        records.append(record)
        segment_metadata.append(
            {
                "offset": offset,
                "history": f"{path.parent.name}/{path.name}",
                "first_total_iteration": int(iteration[0]),
                "last_total_iteration": int(iteration[-1]),
            }
        )

    history = np.concatenate(records)
    if np.any(np.diff(history["iteration"]) <= 0):
        raise RuntimeError("assembled total-iteration coordinate is not strictly increasing")

    figure, axes = plt.subplots(2, 1, figsize=(11, 8), sharex=True, constrained_layout=True)
    ax = axes[0]
    colors = {"CL": "#1261a0", "CD": "#d95f02", "CM_y": "#6a3d9a"}
    references = {"CL": args.nasa_cl, "CD": args.nasa_cd, "CM_y": args.nasa_cm}
    labels = {"CL": "$C_L$", "CD": "$C_D$", "CM_y": "$C_M$"}
    for field in ("CL", "CD", "CM_y"):
        ax.plot(history["iteration"], history[field], color=colors[field], lw=1.7, label=labels[field])
        ax.axhline(references[field], color=colors[field], ls="--", lw=1.0, alpha=0.65)
    ax.set_ylabel("aerodynamic coefficient")
    ax.set_title("NASA CRM Cartesian-AMR steady pseudo-iteration history")
    ax.grid(alpha=0.25)
    ax.legend(ncol=3)

    ax = axes[1]
    ax.semilogy(history["iteration"], history["global"], color="#1b9e77", lw=1.7, label="global max residual")
    ax.semilogy(history["iteration"], history["near"], color="#e7298a", lw=1.7, label="near-body max residual")
    ax.set_xlabel("cumulative pseudo-iteration")
    ax.set_ylabel("normalized residual")
    ax.grid(alpha=0.25, which="both")
    ax.legend()
    figure.suptitle("Dashed load lines: NASA USM3D reference; trajectory is not fixed-point converged")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.out, dpi=220)
    plt.close(figure)

    metadata = {
        "description": "continuous restart-segment history; steady pseudo-iterations, not physical-time samples",
        "segments": segment_metadata,
        "points": int(history.size),
        "first_total_iteration": int(history["iteration"][0]),
        "last_total_iteration": int(history["iteration"][-1]),
        "endpoint": {field: float(history[field][-1]) for field in ("CD", "CL", "CM_y", "global", "near")},
        "NASA_USM3D": {"CD": args.nasa_cd, "CL": args.nasa_cl, "CM_y": args.nasa_cm},
    }
    args.out.with_suffix(".json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
