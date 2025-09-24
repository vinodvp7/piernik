#!/usr/bin/env python3
"""
Compare 'diffusion' vs 'streaming' runs by stitching PIERNIK .h5 blocks in CWD.

If both exist: draws a 2x2 figure (diffusion row, streaming row).
If only one exists: draws a 1x2 figure and prints a message prompting the other run.

Colorbars/ticks/ranges are fully automatic (no manual vmin/vmax or tick settings).
If the right-panel field is strictly positive, a LogNorm is used automatically;
otherwise a linear scale is used.

Axes limits are forced to x,y in [-0.5, 0.5].
"""

from __future__ import annotations

import argparse
import glob
import os
import sys
from typing import Optional, Tuple

import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.patches as mpatches


def get_field_names(h5_file: h5py.File) -> list[str]:
    if "data" not in h5_file:
        raise ValueError("HDF5 file is missing group '/data'.")
    try:
        first_block = next(iter(h5_file["data"]))
    except StopIteration as exc:
        raise ValueError("HDF5 group '/data' is empty.") from exc
    return list(h5_file["data"][first_block].keys())


def load_and_stitch_data(fname: str):
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
                slc = np.s_[
                    off[2]: off[2] + dims[2],
                    off[1]: off[1] + dims[1],
                    off[0]: off[0] + dims[0],
                ]
                vol[slc] = block
            stitched_data[field] = vol

    return stitched_data, global_cell_dims, origin, spacing


def parse_time_from_name(fname: str, time_scale: float) -> Optional[float]:
    base = os.path.basename(fname)
    if base.startswith("scr_") and base.endswith(".h5"):
        stem = base.replace(".h5", "")
        parts = stem.split("_")
        if len(parts) >= 3 and parts[-1].isdigit():
            return int(parts[-1]) * time_scale
    return None


def find_latest(pattern: str) -> Optional[str]:
    files = sorted(glob.glob(pattern))
    return files[-1] if files else None


def build_xy(fname: str) -> Tuple[np.ndarray, np.ndarray]:
    _, cell_dims, origin, spacing = load_and_stitch_data(fname)
    nx = int(cell_dims[0])
    ny = int(cell_dims[1])
    dx = float(spacing[0])
    dy = float(spacing[1])
    x0 = float(origin[0])
    y0 = float(origin[1])
    xe = x0 + np.arange(nx + 1) * dx
    ye = y0 + np.arange(ny + 1) * dy
    x = 0.5 * (xe[:-1] + xe[1:])
    y = 0.5 * (ye[:-1] + ye[1:])
    return x, y


def plot_case(
    ax_l: plt.Axes,
    ax_r: plt.Axes,
    fname: str,
    x: np.ndarray,
    y: np.ndarray,
    label: str,
    field_scalar: str,
    bx_name: str,
    by_name: str,
    time_scale: float,
    fig: plt.Figure,
) -> None:
    stitched, _, _, _ = load_and_stitch_data(fname)

    for needed in ("density", "velocity_x", "velocity_y", field_scalar, bx_name, by_name):
        if needed not in stitched:
            print(
                f"[{label}] Missing field '{needed}' in {fname}. "
                f"Available: {list(stitched.keys())}",
                file=sys.stderr,
            )
            sys.exit(2)

    z0 = 0
    rho = stitched["density"][z0, :, :]
    vx = stitched["velocity_x"][z0, :, :]
    vy = stitched["velocity_y"][z0, :, :]
    scalar = stitched[field_scalar][z0, :, :]
    bx = stitched[bx_name][z0, :, :]
    by = stitched[by_name][z0, :, :]

    t_val = parse_time_from_name(fname, time_scale)
    title = f"t={t_val:g},{label}" if t_val is not None else f"{label} ({os.path.basename(fname)})"

    # Left: density (automatic scaling)
    im_l = ax_l.imshow(
        rho,
        extent=[x[0], x[-1], y[0], y[-1]],
        origin="lower",
        cmap="RdGy_r",
        aspect="auto",
    )
    ax_l.streamplot(x, y, vx, vy, density=0.8, color="k")
    ax_l.set_xlabel(r"$\mathbf{x}$", fontweight="bold")
    ax_l.set_ylabel(r"$\mathbf{y}$", fontweight="bold")
    ax_l.set_title(title, fontweight="bold")
    # Enforce requested axis limits
    ax_l.set_xlim(-0.5, 0.5)
    ax_l.set_ylim(-0.5, 0.5)
    fig.colorbar(im_l, ax=ax_l, fraction=0.046, pad=0.04).set_label(
        r"$\boldsymbol{\rho}$", fontweight="bold", labelpad=6
    )

    # Right: scalar field (automatic scaling; LogNorm if strictly positive)
    norm = LogNorm() if np.all(np.isfinite(scalar) & (scalar > 0.0)) else None
    im_r = ax_r.imshow(
        scalar,
        extent=[x[0], x[-1], y[0], y[-1]],
        origin="lower",
        cmap="RdGy_r",
        norm=norm,
        aspect="auto",
    )
    ax_r.streamplot(x, y, bx, by, density=0.5, color="white")
    ax_r.set_xlabel(r"$\mathbf{x}$", fontweight="bold")
    ax_r.set_ylabel(r"$\mathbf{y}$", fontweight="bold")
    ax_r.set_title(title, fontweight="bold")
    # Enforce requested axis limits
    ax_r.set_xlim(-0.5, 0.5)
    ax_r.set_ylim(-0.5, 0.5)
    fig.colorbar(im_r, ax=ax_r, fraction=0.046, pad=0.04).set_label(
        r"$\mathbf{E_c}$", fontweight="bold", labelpad=6
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Stitch PIERNIK .h5 in CWD and compare diffusion vs streaming panels."
    )
    parser.add_argument("--field", type=str, default="escr_01", help="Scalar field for right panels (default: escr_01).")
    parser.add_argument("--bx", type=str, default="mag_field_x", help="B_x field name.")
    parser.add_argument("--by", type=str, default="mag_field_y", help="B_y field name.")
    parser.add_argument("--time-scale", type=float, default=0.1, help="Index multiplier to label time (default: 0.1).")
    parser.add_argument("--out", type=str, default=None, help="Output PNG (default: <field>.png).")
    args = parser.parse_args()

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
            "figure.figsize": (8, 6),
        }
    )

    dif_file = find_latest("scr_dif_*.h5")
    stm_file = find_latest("scr_stm_*.h5")

    if not dif_file and not stm_file:
        print("No 'scr_dif_*.h5' or 'scr_stm_*.h5' files found in CWD.", file=sys.stderr)
        sys.exit(1)

    ref_file = dif_file or stm_file
    x, y = build_xy(ref_file)

    if dif_file and stm_file:
        fig, axes = plt.subplots(2, 2, figsize=(8, 6))
        (ax00, ax01), (ax10, ax11) = axes
        plot_case(ax00, ax01, dif_file, x, y, "diffusion", args.field, args.bx, args.by, args.time_scale, fig)
        plot_case(ax10, ax11, stm_file, x, y, "streaming", args.field, args.bx, args.by, args.time_scale, fig)
    else:
        fig, (ax_l, ax_r) = plt.subplots(1, 2, figsize=(8, 6))
        if dif_file:
            print("Streaming file missing. Please run the streaming case (scr_stm_*.h5).", file=sys.stderr)
            plot_case(ax_l, ax_r, dif_file, x, y, "diffusion", args.field, args.bx, args.by, args.time_scale, fig)
        else:
            print("Diffusion file missing. Please run the diffusion case (scr_dif_*.h5).", file=sys.stderr)
            plot_case(ax_l, ax_r, stm_file, x, y, "streaming", args.field, args.bx, args.by, args.time_scale, fig)

    plt.tight_layout()

    # Thicken spines of all axes with data
    for ax in fig.axes:
        if isinstance(ax, plt.Axes) and ax.has_data():
            for spine in ax.spines.values():
                spine.set_linewidth(3)

    rect = mpatches.Rectangle(
        (0.005, 0.005), 0.99, 0.99, transform=fig.transFigure, fill=False, lw=4, ec="black"
    )
    fig.add_artist(rect)

    outname = args.out if args.out is not None else f"{args.field}.png"
    plt.savefig(outname, dpi=720)
    print(f"Saved plot to {outname}")


if __name__ == "__main__":
    main()
