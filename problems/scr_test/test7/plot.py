#!/usr/bin/env python3
"""
Stitch PIERNIK .h5 blocks found in the **current directory** and make a 2×2 figure:

  (0,0) empty panel (reserved)
  (0,1) line-outs of E_c at three snapshots
  (1,0) line-outs of density at three snapshots
  (1,1) line-outs of v_x at three snapshots

By default the script looks for files:
  scr_tst_0000.h5, scr_tst_0001.h5, scr_tst_0003.h5
in the current directory. If any are missing, it falls back to the first three
'*.h5' files (sorted).

Time labels are parsed from the filename index times --time-scale (default: 0.02),
so 0001 -> t=0.02, 0003 -> t=0.06, etc.

Usage examples:
  python plot_1d_tripanel.py
  python plot_1d_tripanel.py --field escr_01 --rho density --vx velocity_x \
      --indices 0 1 3 --time-scale 0.02 --xlim -1 1 --out plot.png
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
    stitched_data: dict[str, np.ndarray] = {}
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
    Try to find 'scr_tst_{idx:04d}.h5' for requested indices in CWD.
    If not all are present, fall back to first three '*.h5' files (sorted).
    """
    requested = [f"scr_tst_{i:04d}.h5" for i in indices]
    present = [f for f in requested if os.path.exists(f)]
    if len(present) == len(indices):
        return present

    any_h5 = sorted(glob.glob("*.h5"))
    if len(any_h5) < 3:
        print("Need at least three .h5 files in the current directory.", file=sys.stderr)
        sys.exit(1)
    print("Requested files not all found; falling back to first three *.h5:", any_h5[:3], file=sys.stderr)
    return any_h5[:3]


def parse_time_from_name(fname: str, time_scale: float) -> float | None:
    """Parse 'scr_tst_0003.h5' -> 3 * time_scale; None if no match."""
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
    parser = argparse.ArgumentParser(description="Make 1D line-outs from PIERNIK .h5 files in the current directory.")
    parser.add_argument("--field", type=str, default="escr_01", help="Scalar field for the top-right panel (default: escr_01).")
    parser.add_argument("--rho", type=str, default="density", help="Density field name (default: density).")
    parser.add_argument("--vx", type=str, default="velocity_x", help="Velocity-x field name (default: velocity_x).")
    parser.add_argument("--indices", type=int, nargs=3, default=[0, 1, 3], help="Three snapshot indices, e.g. 0 1 3.")
    parser.add_argument("--time-scale", type=float, default=0.02, help="Multiply file index by this for time labels (default: 0.02).")
    parser.add_argument("--xlim", type=float, nargs=2, default=[-1.0, 1.0], help="x-axis limits (default: -1 1).")
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
    f0, f1, f2 = discover_files(args.indices)
    print(f"Using files: {f0}, {f1}, {f2}")

    # Build x from first file
    first, cell_dims, origin, spacing = load_and_stitch_data(f0)
    nx = int(cell_dims[0]); dx = float(spacing[0]); x0 = float(origin[0])
    xe = x0 + np.arange(nx + 1) * dx
    x = 0.5 * (xe[:-1] + xe[1:])  # centers

    # Sanity checks for fields
    for nm in (args.field, args.rho, args.vx):
        if nm not in first:
            print(f"{f0} missing '{nm}'. Available: {list(first.keys())}", file=sys.stderr)
            sys.exit(2)

    # Load the other two
    second, _, _, _ = load_and_stitch_data(f1)
    third,  _, _, _ = load_and_stitch_data(f2)
    for nm in (args.field, args.rho, args.vx):
        if nm not in second:
            print(f"{f1} missing '{nm}'. Available: {list(second.keys())}", file=sys.stderr)
            sys.exit(2)
        if nm not in third:
            print(f"{f2} missing '{nm}'. Available: {list(third.keys())}", file=sys.stderr)
            sys.exit(2)

    # Time labels
    t0 = parse_time_from_name(f0, args.time_scale)
    t1 = parse_time_from_name(f1, args.time_scale)
    t2 = parse_time_from_name(f2, args.time_scale)
    lbl0 = f"t={t0:g}" if t0 is not None else os.path.basename(f0)
    lbl1 = f"t={t1:g}" if t1 is not None else os.path.basename(f1)
    lbl2 = f"t={t2:g}" if t2 is not None else os.path.basename(f2)

    fig, ax = plt.subplots(2, 2, figsize=(8, 6))

    # (0,0) empty / reserved
    ax[0, 0].axis("off")

    # (0,1) E_c at three times
    ec0 = first[args.field][0, 0, :]
    ec1 = second[args.field][0, 0, :]
    ec2 = third[args.field][0, 0, :]
    ax[0, 1].plot(x, ec0, linewidth=1.6, color="k", label=lbl0)
    ax[0, 1].plot(x, ec1, linewidth=1.6, color="r", label=lbl1)
    ax[0, 1].plot(x, ec2, linewidth=1.6, color="g", label=lbl2)
    ax[0, 1].set_xlabel(r"$\mathbf{x}$", fontweight="bold")
    ax[0, 1].set_ylabel(r"$\mathbf{E_c}$", fontweight="bold")
    ax[0, 1].legend(fontsize="small")
    ax[0, 1].set_xlim(args.xlim[0], args.xlim[1])

    # (1,0) density at three times
    rho0 = first[args.rho][0, 0, :]
    rho1 = second[args.rho][0, 0, :]
    rho2 = third[args.rho][0, 0, :]
    ax[1, 0].plot(x, rho0, linewidth=1.6, color="k", label=lbl0)
    ax[1, 0].plot(x, rho1, linewidth=1.6, color="r", label=lbl1)
    ax[1, 0].plot(x, rho2, linewidth=1.6, color="g", label=lbl2)
    ax[1, 0].set_xlabel(r"$\mathbf{x}$", fontweight="bold")
    ax[1, 0].set_ylabel(r"$\boldsymbol{\rho}$", fontweight="bold")
    ax[1, 0].legend(fontsize="small")
    ax[1, 0].set_xlim(args.xlim[0], args.xlim[1])

    # (1,1) v_x at three times
    vx0 = first[args.vx][0, 0, :]
    vx1 = second[args.vx][0, 0, :]
    vx2 = third[args.vx][0, 0, :]
    ax[1, 1].plot(x, vx0, linewidth=1.6, color="k", label=lbl0)
    ax[1, 1].plot(x, vx1, linewidth=1.6, color="r", label=lbl1)
    ax[1, 1].plot(x, vx2, linewidth=1.6, color="g", label=lbl2)
    ax[1, 1].set_xlabel(r"$\mathbf{x}$", fontweight="bold")
    ax[1, 1].set_ylabel(r"$\mathbf{v_x}$", fontweight="bold")
    ax[1, 1].legend(fontsize="small")
    ax[1, 1].set_xlim(args.xlim[0], args.xlim[1])

    plt.tight_layout()
    for a in (ax[0, 1], ax[1, 0], ax[1, 1]):
        for s in a.spines.values():
            s.set_linewidth(3)

    plt.savefig(args.out, dpi=720)
    print(f"Saved plot to ./{args.out}")


if __name__ == "__main__":
    main()
