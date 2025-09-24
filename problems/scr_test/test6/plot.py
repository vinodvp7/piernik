#!/usr/bin/env python3
"""
Stitch PIERNIK .h5 blocks found in the **current directory** and make a 3×1 plot:

  (0) density line-out at the earliest requested snapshot
  (1) E_c * (v_A)^{p} at a later snapshot, with p=4/3 by default
  (2) E_c line-outs at three snapshots (e.g. indices 1, 2, 10)

Files are discovered in CWD. By default the script looks for:
  scr_tst_0001.h5, scr_tst_0002.h5, scr_tst_0010.h5
and falls back to the first three *.h5 files (sorted) if they are missing.

Time labels are parsed from the filename index times --time-scale (default: 100),
so: 0001 -> t=100, 0002 -> t=200, 0010 -> t=1000.

Usage examples:
  python plot_1d_profiles.py
  python plot_1d_profiles.py --field escr_01 --bx mag_field_x --indices 1 5 20 --time-scale 0.1 --xlim 0 600 --out plot.png
"""

from __future__ import annotations

import argparse
import glob
import os
import sys
from typing import List, Tuple

import h5py
import numpy as np
import matplotlib.pyplot as plt


# ------------------------------- I/O utilities -------------------------------

def get_field_names(h5_file: h5py.File) -> List[str]:
    """Return dataset names from the first block under '/data'."""
    if "data" not in h5_file:
        raise ValueError("HDF5 file is missing group '/data'.")
    try:
        first_block = next(iter(h5_file["data"]))
    except StopIteration as exc:
        raise ValueError("HDF5 group '/data' is empty.") from exc
    return list(h5_file["data"][first_block].keys())


def load_and_stitch_data(fname: str):
    """
    Load and stitch 3D cell-centered fields from all blocks.

    Returns
    -------
    stitched_data : dict[str, np.ndarray]  # field -> (Nz, Ny, Nx)
    global_cell_dims : np.ndarray          # [Nx, Ny, Nz]
    origin : np.ndarray                    # [x0, y0, z0]
    spacing : np.ndarray                   # [dx, dy, dz]
    """
    stitched_data = {}
    with h5py.File(fname, "r") as f:
        field_names = get_field_names(f)
        block_names = list(f["data"].keys())

        all_off = np.array([f["data"][n].attrs["off"] for n in block_names], dtype=int)
        all_dims = np.array([f["data"][n].attrs["n_b"] for n in block_names], dtype=int)

        global_cell_dims = np.max(all_off + all_dims, axis=0).astype(int)

        origin = np.array(f["simulation_parameters"].attrs["domain_left_edge"], dtype=float)
        domain_right = np.array(f["simulation_parameters"].attrs["domain_right_edge"], dtype=float)
        spacing = (domain_right - origin) / np.maximum(1.0, global_cell_dims.astype(float))

        # On-disk arrays are (z, y, x); stitch in that order.
        for field in field_names:
            vol = np.empty((global_cell_dims[2], global_cell_dims[1], global_cell_dims[0]), dtype=np.float32)
            for i, name in enumerate(block_names):
                block = f["data"][name][field][:]
                off = all_off[i]
                dims = all_dims[i]
                slc = np.s_[off[2]: off[2] + dims[2], off[1]: off[1] + dims[1], off[0]: off[0] + dims[0]]
                vol[slc] = block
            stitched_data[field] = vol

    return stitched_data, global_cell_dims, origin, spacing


def discover_files(indices: List[int]) -> List[str]:
    """
    Try to find 'scr_tst_{idx:04d}.h5' for requested indices.
    If not all are present, fall back to first three *.h5 files (sorted).
    """
    requested = [f"scr_tst_{i:04d}.h5" for i in indices]
    present = [f for f in requested if os.path.exists(f)]
    if len(present) == len(indices):
        return present

    # Fallback: any *.h5 in CWD
    any_h5 = sorted(glob.glob("*.h5"))
    if len(any_h5) < 3:
        print("Need at least three .h5 files in the current directory.", file=sys.stderr)
        sys.exit(1)
    print("Requested files not all found; falling back to first three *.h5:", any_h5[:3], file=sys.stderr)
    return any_h5[:3]


def parse_time_from_name(fname: str, time_scale: float) -> float | None:
    """Parse 'scr_tst_0002.h5' -> 2 * time_scale; None if no match."""
    base = os.path.basename(fname)
    if base.startswith("scr_tst_") and base.endswith(".h5"):
        num = base.replace("scr_tst_", "").replace(".h5", "")
        try:
            return int(num) * time_scale
        except Exception:
            return None
    return None


# ---------------------------------- Main -------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description="Make 1D profiles from PIERNIK .h5 files in the current directory.")
    parser.add_argument("--field", type=str, default="escr_01", help="Scalar field for panels (default: escr_01).")
    parser.add_argument("--rho", type=str, default="density", help="Density field name (default: density).")
    parser.add_argument("--bx", type=str, default="mag_field_x", help="Magnetic field component to use for v_A (default: mag_field_x).")
    parser.add_argument("--alfven-exp", type=float, default=4.0 / 3.0, help="Exponent p in E_c * (v_A)^p (default: 4/3).")
    parser.add_argument("--indices", type=int, nargs=3, default=[1, 2, 10], help="Three snapshot indices, e.g. 1 2 10.")
    parser.add_argument("--time-scale", type=float, default=100.0, help="Multiply file index by this for time labels (default: 100).")
    parser.add_argument("--xlim", type=float, nargs=2, default=[0.0, 600.0], help="x-axis limits (default: 0 600).")
    parser.add_argument("--out", type=str, default="plot.png", help="Output PNG filename (default: plot.png).")
    args = parser.parse_args()

    plt.rcParams.update({
        "figure.dpi": 150, "savefig.dpi": 720,
        "axes.linewidth": 2.5, "font.size": 12,
        "xtick.major.size": 5, "ytick.major.size": 5,
        "xtick.major.width": 1.0, "ytick.major.width": 1.0,
        "xtick.direction": "out", "ytick.direction": "out",
        "legend.frameon": True, "legend.edgecolor": "black", "legend.framealpha": 1.0,
        "figure.figsize": (8, 6),
    })

    # Find three files in CWD
    f1, f2, f3 = discover_files(args.indices)

    # Build x from first file
    first, cell_dims, origin, spacing = load_and_stitch_data(f1)
    nx = int(cell_dims[0])
    dx = float(spacing[0])
    x0 = float(origin[0])
    xe = x0 + np.arange(nx + 1) * dx
    x = 0.5 * (xe[:-1] + xe[1:])  # centers

    # Safety checks
    needed0 = [args.rho, args.field]
    for nm in needed0:
        if nm not in first:
            print(f"{f1} missing '{nm}'. Available: {list(first.keys())}", file=sys.stderr)
            sys.exit(2)

    fig, ax = plt.subplots(3, 1, figsize=(8, 6), sharex=True)

    # Panel 0: density at f1
    rho1 = first[args.rho][0, 0, :]
    ax[0].plot(x, rho1, color="k", linewidth=1.5)
    t1 = parse_time_from_name(f1, args.time_scale)
    label_t1 = f"t={t1:g}" if t1 is not None else os.path.basename(f1)
    ax[0].set_ylabel(r"$\boldsymbol{\rho}$", fontweight="bold", labelpad=10)
    ax[0].set_title(label_t1, fontsize=11)

    # Panel 1: E_c * (v_A)^p at f3
    third, _, _, _ = load_and_stitch_data(f3)
    for nm in (args.field, args.bx, args.rho):
        if nm not in third:
            print(f"{f3} missing '{nm}'. Available: {list(third.keys())}", file=sys.stderr)
            sys.exit(2)
    ec3 = third[args.field][0, 0, :]
    bx3 = third[args.bx][0, 0, :]
    rho3 = third[args.rho][0, 0, :]

    # Alfvén speed proxy v_A ~ Bx / sqrt(rho); no 4π factor (consistent with code units)
    with np.errstate(divide="ignore", invalid="ignore"):
        vA3 = bx3 / np.sqrt(np.maximum(rho3, 1e-30))
    prof = ec3 * np.power(vA3, args.alfven_exp)

    t3 = parse_time_from_name(f3, args.time_scale)
    label_t3 = f"t={t3:g}" if t3 is not None else os.path.basename(f3)
    ax[1].plot(x, prof, color="k", linewidth=1.6, label=label_t3)
    ax[1].set_ylabel(r"$\mathbf{E_c\,v_A^{\,p}}$", fontweight="bold", labelpad=4)
    ax[1].legend(fontsize="small")
    # Simple dynamic ticks (4 evenly spaced)
    ymin, ymax = float(np.nanmin(prof)), float(np.nanmax(prof))
    if np.isfinite(ymin) and np.isfinite(ymax) and ymax > ymin:
        ax[1].set_yticks(np.linspace(ymin, ymax, 4))

    # Panel 2: E_c at f1, f2, f3
    second, _, _, _ = load_and_stitch_data(f2)
    for nm in (args.field,):
        if nm not in second:
            print(f"{f2} missing '{nm}'. Available: {list(second.keys())}", file=sys.stderr)
            sys.exit(2)

    ec1 = first[args.field][0, 0, :]
    ec2 = second[args.field][0, 0, :]
    ec3 = third[args.field][0, 0, :]

    t2 = parse_time_from_name(f2, args.time_scale)
    label_t2 = f"t={t2:g}" if t2 is not None else os.path.basename(f2)

    ax[2].plot(x, ec1, color="k", linewidth=1.5, label=label_t1)
    ax[2].plot(x, ec2, color="r", linewidth=1.5, linestyle="dashed", label=label_t2)
    ax[2].plot(x, ec3, color="b", linewidth=1.5, linestyle="dashed", label=label_t3)
    ax[2].legend(fontsize="small")
    ax[2].set_xlabel(r"$\mathbf{x}$", fontweight="bold")
    ax[2].set_ylabel(r"$\mathbf{E_c}$", fontweight="bold", labelpad=12)

    # Common formatting
    ax[2].set_xlim(args.xlim[0], args.xlim[1])
    for a in ax:
        for s in a.spines.values():
            s.set_linewidth(2.0)

    plt.tight_layout()
    plt.savefig(args.out, dpi=720)
    print(f"Saved plot to ./{args.out}")


if __name__ == "__main__":
    main()
