#!/usr/bin/env python3
"""Minimal airfoil CFD interface with locked, verified numerical presets."""

from __future__ import annotations

import argparse
from argparse import Namespace
import json
import os
from pathlib import Path
import tempfile

import quick_cfd


ROOT = Path(__file__).resolve().parent
MESH_PRESETS = {4700: "LF1", 9500: "LF2"}


def geometry_format(path: Path) -> str:
    with path.open(encoding="utf-8-sig") as stream:
        header = [value.strip().lower() for value in stream.readline().split(",")]
    if header == ["x", "y"]:
        return "gxxx"
    if header == ["x", "y", "patch"]:
        return "polygon"
    raise ValueError("geometry header must be x,y or x,y,patch")


def preset_patch(algorithm: str, mesh_cells: int) -> dict:
    common = {
        "spatial_order": 2,
        "preconditioning": "turkel",
        "mach_floor": 0.02,
        "adaptive_cfl": {"enabled": False},
        "mode": "steady",
        "stop_time_star": 150.0,
        "convergence": {
            "metric": "cl_cm_cp",
            "raw_tolerance": 0.01,
            "coefficient_change": 0.005,
            "cp_change": 0.005,
        },
    }
    if algorithm == "rk3":
        common.update(
            integrator="rk3",
            backend="cuda",
            local_pseudo_time=True,
            budget_seconds=55.0,
        )
    elif algorithm == "dplur":
        is_fine = mesh_cells == 9500
        common.update(
            integrator="dplur",
            backend="cpu",
            local_pseudo_time=False,
            budget_seconds=75.0 if is_fine else 15.0,
            dplur={
                "CFL": 1.0 if is_fine else 4.0,
                "sweeps": 5,
                "omega": 0.8,
                "block_coupling": True,
                "maximum_iterations": 10000 if is_fine else 5000,
                "check_interval": 100,
            },
        )
    else:
        raise ValueError("algorithm must be rk3 or dplur")
    return {
        "solver": common,
        "output": {"formats": ["npz"], "audit_budgets": False},
    }


def run(arguments: argparse.Namespace) -> int:
    geometry = Path(arguments.geometry).resolve()
    if not geometry.is_file():
        raise ValueError(f"geometry file not found: {geometry}")
    output = Path(arguments.out).resolve()
    profile = MESH_PRESETS[arguments.mesh_cells]
    patch = preset_patch(arguments.algorithm, arguments.mesh_cells)
    patch["freestream"] = {
        "Mach": arguments.mach,
        "AOA_degrees": arguments.aoa,
        "p_Pa": arguments.pressure_pa,
        "T_K": arguments.temperature_k,
    }

    with tempfile.TemporaryDirectory(prefix="quick_cfd_simple_") as temporary:
        config_path = Path(temporary) / "preset.json"
        config_path.write_text(json.dumps(patch, indent=2) + "\n")
        solver_arguments = Namespace(
            command="run",
            points=str(geometry),
            config=str(config_path),
            out=str(output),
            executable=os.environ.get("QUICK_CFD_EXECUTABLE"),
            point_format=geometry_format(geometry),
            levels=[profile],
            mach=None,
            aoa=None,
            pressure=None,
            temperature=None,
            chord=None,
            cfl=None,
            dt=None,
            stop_time=None,
            physical_time=None,
            mode=None,
            backend=None,
            turbulence=None,
            viscosity=None,
            fields=None,
            formats=None,
            region=None,
            surface=None,
            loads=None,
            volume=None,
            load_fields=None,
        )
        status = quick_cfd.run(solver_arguments)

    pipeline = json.loads((output / "pipeline.json").read_text())
    result = pipeline["cases"][0]
    iterations = result.get("pseudo_iteration", 0) or result["steps"]
    summary = {
        "interface": "simple_airfoil_v1",
        "user_inputs": {
            "geometry": geometry.name,
            "algorithm": arguments.algorithm,
            "mesh_cells_nominal": arguments.mesh_cells,
            "Mach": arguments.mach,
            "AOA_degrees": arguments.aoa,
            "pressure_Pa": arguments.pressure_pa,
            "temperature_K": arguments.temperature_k,
        },
        "locked_preset": {
            "profile": profile,
            "backend": patch["solver"]["backend"],
            "spatial_order": 2,
            "turbulence": "sa-neg",
            "preconditioning": "turkel",
            "convergence_gate": "raw residual plus stable CL/CM/Cp for five reports",
        },
        "result": {
            "actual_control_volumes": result["fluid_control_volumes"],
            "iterations": iterations,
            "process_seconds": result["process_seconds"],
            "pipeline_seconds": pipeline["pipeline_seconds"],
            "converged": result["converged"],
            "stop_reason": result["stop_reason"],
            "CL": result["Cl"],
            "CD": result["Cd"],
            "CM_CCW_positive": result["Cm"],
            "Re_chord": pipeline["scales"]["Re_chord"],
        },
    }
    (output / "simple_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    value = summary["result"]
    print(
        f"simple result: {value['actual_control_volumes']} CV, "
        f"{value['iterations']} iterations, {value['process_seconds']:.3f} s, "
        f"CL/CD/CM={value['CL']:.6f}/{value['CD']:.6f}/{value['CM_CCW_positive']:.6f}, "
        f"converged={value['converged']}"
    )
    print(output / "simple_summary.json")
    return status


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--geometry", required=True, help="x,y GXXX CSV or ordered x,y,patch polygon CSV")
    parser.add_argument("--algorithm", required=True, choices=["rk3", "dplur"], help="rk3 selects CUDA; dplur selects the CPU reference")
    parser.add_argument("--mesh-cells", required=True, type=int, choices=sorted(MESH_PRESETS), help="nominal cell preset; actual cut-cell count is reported")
    parser.add_argument("--mach", required=True, type=float)
    parser.add_argument("--aoa", required=True, type=float, help="angle of attack in degrees")
    parser.add_argument("--pressure-pa", required=True, type=float)
    parser.add_argument("--temperature-k", required=True, type=float)
    parser.add_argument("--out", required=True, help="new output directory")
    return run(parser.parse_args())


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, KeyError, TypeError, OSError, RuntimeError) as error:
        print(f"ERROR: {error}")
        raise SystemExit(2)
