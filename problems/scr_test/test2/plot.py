#!/usr/bin/env python3
"""
Stitch PIERNIK .h5 blocks in CWD and plot a 1D line-out with an analytical curve.

By default:
  - searches for 'scr_tst_*.h5' (or any '*.h5' if none) in the current directory,
  - plots field 'escr_01' from z=0, y=0 at the first two times found,
  - overlays an analytical curve at the second time,
  - writes <field>.png unless --out is given.

Examples
--------
# Use defaults (field escr_01, time-scale 0.06)
python plot.py

# Explicit field and time-scale, custom output file
python plot.py --field escr_01 --time-scale 0.06 --out plot.png

# Force analytical time if filenames don't encode it
python plot.py --analytic-time 0.06
"""

from __future__ import annotations

import argparse
import glob
import os
import sys

import h5py
import numpy as np
import matplotlib.pyplot as plt


def get_field_names(h5_file: h5py.File) -> list[str]:
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
    stitched_data : dict[str, np.ndarray]
        Field -> ndarray (Nz, Ny, Nx).
    global_cell_dims : np.ndarray
        Integer array [Nx, Ny, Nz].
    origin : np.ndarray
        Physical left edge [x0, y0, z0].
    spacing : np.ndarray
        Cell spacing [dx, dy, dz].
    """
    stitched_data: dict[str, np.ndarray] = {}
    with h5py.File(fname, "r") as f:
        field_names = get_field_names(f)
        block_names = list(f["data"].keys())

        all_off = np.array([f["data"][n].attrs["off"] for n in block_names], dtype=int)
        all_dims = np.array([f["data"][n].attrs["n_b"] for n in block_names], dtype=int)

        global_cell_dims = np.max(all_off + all_dims, axis=0).astype(int)

        origin = np.array(
            f["simulation_parameters"].attrs["domain_left_edge"], dtype=float
        )
        domain_right = np.array(
            f["simulation_parameters"].attrs["domain_right_edge"], dtype=float
        )
        spacing = (domain_right - origin) / np.maximum(
            1.0, global_cell_dims.astype(float)
        )

        for field in field_names:
            vol = np.empty(
                (global_cell_dims[2], global_cell_dims[1], global_cell_dims[0]),
                dtype=np.float32,
            )
            for i, name in enumerate(block_names):
                block = f["data"][name][field][:]
                off = all_off[i]
                dims = all_dims[i]
                slc = np.s_[off[2]: off[2] + dims[2], off[1]: off[1] + dims[1], off[0]: off[0] + dims[0]]
                vol[slc] = block
            stitched_data[field] = vol

    return stitched_data, global_cell_dims, origin, spacing


def pick_input_files() -> list[str]:
    """Return at least two .h5 files: prefer 'scr_tst_*.h5', else any '*.h5'."""
    files = sorted(glob.glob("scr_tst_*.h5"))
    if not files:
        files = sorted(glob.glob("*.h5"))
    if len(files) < 2:
        print("Need at least two HDF5 files in CWD.", file=sys.stderr)
        sys.exit(1)
    return files[:2]


def parse_time_from_name(fname: str, time_scale: float) -> float | None:
    """
    Parse time from 'scr_tst_0001.h5' -> 1 * time_scale.
    Returns None if pattern doesn't match.
    """
    base = os.path.basename(fname)
    if base.startswith("scr_tst_") and base.endswith(".h5"):
        num = base.replace("scr_tst_", "").replace(".h5", "")
        try:
            return int(num) * time_scale
        except Exception:
            return None
    return None


def analytic_profile(x: np.ndarray, t: float) -> np.ndarray:
    """
    Piecewise analytic: for |x| <= x_m use constant plateau,
    else linear decay with |x|.

    x_m = sqrt((1 + 4/3 t)^2 + 8/3 t - 1)
    y   = 2 + 4/3 t - x_m          for |x| <= x_m
          2 + 4/3 t - |x|          for |x| >  x_m
    """
    xm = np.sqrt((1.0 + (4.0 / 3.0) * t) ** 2 + (8.0 / 3.0) * t - 1.0)
    y = np.full_like(x, 2.0 + (4.0 / 3.0) * t - xm)
    mask = (x < -xm) | (x > xm)
    y[mask] = 2.0 + (4.0 / 3.0) * t - np.abs(x[mask])
    return y


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Stitch PIERNIK .h5 outputs in CWD and plot 1D line-out vs. analytic."
    )
    parser.add_argument(
        "--field", type=str, default="escr_01", help="Field name to plot (default: escr_01)."
    )
    parser.add_argument(
        "--time-scale",
        type=float,
        default=0.06,
        help="Multiply file index by this to interpret time (default: 0.06).",
    )
    parser.add_argument(
        "--analytic-time",
        type=float,
        default=None,
        help="Override time for analytic curve (default: parsed from second file).",
    )
    parser.add_argument(
        "--out", type=str, default=None, help="Output PNG filename (default: <field>.png)."
    )
    args = parser.parse_args()

    files = pick_input_files()

    # Build x from the first file grid
    first, cell_dims, origin, spacing = load_and_stitch_data(files[0])
    nx = int(cell_dims[0])
    dx = float(spacing[0])
    x0 = float(origin[0])
    xe = x0 + np.arange(nx + 1) * dx
    x = 0.5 * (xe[:-1] + xe[1:])

    field = args.field
    if field not in first:
        print(
            f"Field '{field}' not found in {files[0]}. "
            f"Available: {list(first.keys())}",
            file=sys.stderr,
        )
        sys.exit(2)

    # Plot configuration
    plt.rcParams.update(
        {
            "figure.dpi": 150,
            "savefig.dpi": 720,
            "axes.linewidth": 2.5,
            "font.size": 12,
            "xtick.major.size": 5,
            "ytick.major.size": 5,
            "xtick.major.width": 1.0,
            "ytick.major.width": 1.0,
            "xtick.direction": "out",
            "ytick.direction": "out",
            "legend.frameon": True,
            "legend.edgecolor": "black",
            "legend.framealpha": 1.0,
            "figure.figsize": (8, 6),
        }
    )

    plt.clf()
    ax = plt.gca()

    # Numerical: t0 and t1 from first two files
    colors = ["k", "r"]
    labels = []
    y_curves = []
    for i, fname in enumerate(files[:2]):
        stitched, _, _, _ = load_and_stitch_data(fname)
        y = stitched[field][0, 0, :]
        t_val = parse_time_from_name(fname, args.time_scale)
        label = f"t={t_val:g}" if t_val is not None else os.path.basename(fname)
        ax.plot(x, y, label=label, linewidth=1.8, color=colors[i % 2])
        labels.append(label)
        y_curves.append(y)

    # Analytic at the second time
    t_analytic = args.analytic_time
    if t_analytic is None:
        t_analytic = parse_time_from_name(files[1], args.time_scale) or 0.0
    y_ana = analytic_profile(x, t_analytic)
    ax.plot(
        x,
        y_ana,
        label=f"t={t_analytic:g} (analytical)",
        linewidth=1.8,
        color="b",
        linestyle="dashed",
    )

    # Error norms vs. the second numerical curve, if present
    if len(y_curves) >= 2:
        y_num = y_curves[1]
        l1 = float(np.mean(np.abs(y_num - y_ana)))
        l2 = float(np.sqrt(np.mean((y_num - y_ana) ** 2)))
        msg = rf"$L_1 = {l1:.3e}$" + "\n" + rf"$L_2 = {l2:.3e}$"
        ax.text(
            0.02,
            0.98,
            msg,
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=11,
            bbox=dict(
                boxstyle="round,pad=0.30",
                facecolor="white",
                edgecolor="black",
                lw=0.8,
                alpha=0.9,
            ),
            zorder=5,
        )

    ax.set_xlabel(r"$\mathbf{x}$", fontweight="bold", labelpad=6)
    ax.set_ylabel(r"$\mathbf{E_c}$", fontweight="bold", labelpad=6)
    ax.set_xlim(x.min(), x.max())
    ax.legend()
    plt.tight_layout()

    outname = args.out if args.out is not None else f"{field}.png"
    plt.savefig(outname, dpi=720)
    print(f"Saved plot to {outname}")


if __name__ == "__main__":
    main()
