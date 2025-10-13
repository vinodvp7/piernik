#!/usr/bin/env python3
"""
Stitch PIERNIK .h5 blocks in CWD and plot a 2×2 panel with streamlines.

By default:
  - searches for 'scr_tst_*.h5' (or any '*.h5' if none) in the current directory,
  - plots scalar field 'escr_01' at the first three times found,
  - overlays streamlines from 'mag_field_x'/'mag_field_y' (z=0),
  - writes <field>.png unless --out is given.

Examples
--------
python plot.py
python plot.py --field escr_01 --bx mag_field_x --by mag_field_y --time-scale 0.2
"""

from __future__ import annotations

import argparse
import glob
import os
import sys

import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.patches as mpatches


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

        # On-disk arrays are (z, y, x); stitch in that order.
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
    """Return up to three .h5 files: prefer 'scr_tst_*.h5', else any '*.h5'."""
    files = sorted(glob.glob("scr_tst_*.h5"))
    if not files:
        files = sorted(glob.glob("*.h5"))
    if not files:
        print("No HDF5 files found in CWD.", file=sys.stderr)
        sys.exit(1)
    return files[:3]


def parse_time_from_name(fname: str, time_scale: float) -> float | None:
    """
    Parse time from 'scr_tst_0002.h5' -> 2 * time_scale.
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


def build_xy(first_fname: str):
    """Construct cell-centered x,y from the first file's grid."""
    _, cell_dims, origin, spacing = load_and_stitch_data(first_fname)
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


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Stitch PIERNIK .h5 outputs in CWD and plot scalar + streamlines."
    )
    parser.add_argument(
        "--field", type=str, default="escr_01", help="Scalar field to imshow (default: escr_01)."
    )
    parser.add_argument(
        "--bx", type=str, default="mag_field_x", help="X-component field for streamplot."
    )
    parser.add_argument(
        "--by", type=str, default="mag_field_y", help="Y-component field for streamplot."
    )
    parser.add_argument(
        "--time-scale",
        type=float,
        default=0.2,
        help="Multiply file index by this to label time (default: 0.2).",
    )
    parser.add_argument(
        "--vmin", type=float, default=1e-6, help="LogNorm vmin for imshow (default: 1e-6)."
    )
    parser.add_argument(
        "--vmax", type=float, default=1.0, help="LogNorm vmax for imshow (default: 1)."
    )
    parser.add_argument(
        "--out", type=str, default=None, help="Output PNG (default: <field>.png)."
    )
    args = parser.parse_args()

    # Style
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

    files = pick_input_files()
    x, y = build_xy(files[0])

    fig, axes = plt.subplots(2, 2, figsize=(8, 6))
    axes = np.asarray(axes)
    panels = [(0, 0), (0, 1), (1, 0)]
    titles: list[str] = []

    for idx, (ii, jj) in enumerate(panels):
        if idx >= len(files):
            break
        fname = files[idx]
        stitched, _, _, _ = load_and_stitch_data(fname)

        for needed in (args.field, args.bx, args.by):
            if needed not in stitched:
                print(
                    f"Field '{needed}' not found in {fname}. "
                    f"Available: {list(stitched.keys())}",
                    file=sys.stderr,
                )
                sys.exit(2)

        z0 = 0  # z-index for the slice
        scalar = stitched[args.field][z0, :, :]
        bx = stitched[args.bx][z0, :, :]
        by = stitched[args.by][z0, :, :]

        ax = axes[ii, jj]
        im = ax.imshow(
            scalar,
            extent=[x[0], x[-1], y[0], y[-1]],
            origin="lower",
            cmap="RdGy_r",
            norm=LogNorm(vmin=args.vmin, vmax=args.vmax),
            aspect="auto",
        )
        ax.streamplot(x, y, bx, by, density=0.1, color="b")

        t_val = parse_time_from_name(fname, args.time_scale)
        if t_val is None:
            title = os.path.basename(fname)
        else:
            title = f"t={t_val:g}"
        titles.append(title)

        ax.set_xlabel(r"$\mathbf{x}$", fontweight="bold")
        ax.set_ylabel(r"$\mathbf{y}$", fontweight="bold")
        ax.set_title(title, fontweight="bold")

        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cb.set_label(r"$\mathbf{E}_c$", fontweight="bold", labelpad=6)

    # Turn off the empty fourth panel
    axes[1, 1].axis("off")

    # Layout and cosmetics
    plt.tight_layout()
    for a in (axes[0, 0], axes[0, 1], axes[1, 0]):
        for spine in a.spines.values():
            spine.set_linewidth(3)

    # Bold border around the whole figure
    rect = mpatches.Rectangle(
        (0.005, 0.005),
        0.99,
        0.99,
        transform=fig.transFigure,
        fill=False,
        lw=4,
        ec="black",
        zorder=10,
    )
    fig.add_artist(rect)

    outname = args.out if args.out is not None else f"{args.field}.png"
    plt.savefig(outname, dpi=720)
    print(f"Saved plot to {outname}")


if __name__ == "__main__":
    main()
