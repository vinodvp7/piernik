#!/usr/bin/env python3
"""
Stitch PIERNIK .h5 blocks from the current directory and make a 1D plot.

By default looks for files matching 'scr_tst_*.h5' (or any '*.h5' if none),
stitches blocks, extracts a (z=0, y=0) line-out of the requested field, and
saves a PNG.

Examples
--------
# default field 'escr_01', auto-picks up to 3 files, saves escr_01.png
python plot.py

# explicit field and output name; time scale 0.05 -> 0001 -> t=0.05
python plot.py --field escr_01 --time-scale 0.05 --out plot.png
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
    """
    Return available dataset names from the first block under '/data'.

    Parameters
    ----------
    h5_file : h5py.File
        Open HDF5 file handle.

    Returns
    -------
    list[str]
        Field names present in the first data block.
    """
    if "data" not in h5_file:
        raise ValueError("HDF5 file is missing group '/data'.")
    try:
        first_block_name = next(iter(h5_file["data"]))
    except StopIteration as exc:
        raise ValueError("HDF5 group '/data' is empty.") from exc
    return list(h5_file["data"][first_block_name].keys())


def load_and_stitch_data(fname: str):
    """
    Load and stitch 3D cell-centered fields from all blocks.

    Parameters
    ----------
    fname : str
        Path to HDF5 file.

    Returns
    -------
    stitched_data : dict[str, np.ndarray]
        Maps field name -> ndarray shaped (Nz, Ny, Nx).
    global_cell_dims : np.ndarray
        Integer array [Nx, Ny, Nz] of cell counts.
    origin : np.ndarray
        Physical left edge [x0, y0, z0].
    spacing : np.ndarray
        Cell spacing [dx, dy, dz].
    """
    stitched_data: dict[str, np.ndarray] = {}
    with h5py.File(fname, "r") as f:
        field_names = get_field_names(f)
        block_names = list(f["data"].keys())

        all_offsets = np.array(
            [f["data"][name].attrs["off"] for name in block_names], dtype=int
        )  # [nblocks, 3] (Nx, Ny, Nz offsets)
        all_dims = np.array(
            [f["data"][name].attrs["n_b"] for name in block_names], dtype=int
        )  # [nblocks, 3] (block sizes)

        # Global cell counts per axis (Nx, Ny, Nz)
        global_cell_dims = np.max(all_offsets + all_dims, axis=0).astype(int)

        # Domain geometry
        origin = np.array(
            f["simulation_parameters"].attrs["domain_left_edge"], dtype=float
        )
        domain_right = np.array(
            f["simulation_parameters"].attrs["domain_right_edge"], dtype=float
        )
        domain_size = domain_right - origin
        spacing = domain_size / np.maximum(1.0, global_cell_dims.astype(float))

        # On disk arrays are (z, y, x); build stitched volume in that order
        for field in field_names:
            stitched_field = np.empty(
                (global_cell_dims[2], global_cell_dims[1], global_cell_dims[0]),
                dtype=np.float32,
            )
            for i, name in enumerate(block_names):
                block = f["data"][name][field][:]
                off = all_offsets[i]
                dims = all_dims[i]
                slc = np.s_[off[2]: off[2] + dims[2], off[1]: off[1] + dims[1], off[0]: off[0] + dims[0]]
                stitched_field[slc] = block
            stitched_data[field] = stitched_field

    return stitched_data, global_cell_dims, origin, spacing


# Preferred rc parameters (readability; does not affect linting)
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


def pick_input_files() -> list[str]:
    """
    Choose which .h5 files to plot from current directory.

    Priority
    --------
    1) 'scr_tst_*.h5' (sorted)
    2) any '*.h5' (sorted)

    Returns up to 3 files, preferring canonical indices 0000, 0001, 0002 if present.
    """
    candidates = sorted(glob.glob("scr_tst_*.h5"))
    if not candidates:
        candidates = sorted(glob.glob("*.h5"))
    if not candidates:
        print("No .h5 files found in current directory.", file=sys.stderr)
        sys.exit(1)

    preferred = ["scr_tst_0000.h5", "scr_tst_0001.h5", "scr_tst_0002.h5"]
    chosen: list[str] = [c for c in preferred if c in candidates]
    if not chosen:
        chosen = candidates[:3]
    else:
        for c in candidates:
            if len(chosen) >= 3:
                break
            if c not in chosen:
                chosen.append(c)
    return chosen


def parse_time_from_name(fname: str, time_scale: float) -> float | None:
    """
    Parse time from filename like 'scr_tst_0001.h5' -> 1 * time_scale.

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


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Stitch PIERNIK .h5 outputs in CWD and plot a 1D line-out."
    )
    parser.add_argument(
        "--field",
        type=str,
        default="escr_01",
        help="Field name to plot (default: escr_01).",
    )
    parser.add_argument(
        "--time-scale",
        type=float,
        default=0.05,
        help="Multiply file index by this to label time (default: 0.05).",
    )
    parser.add_argument(
        "--out",
        type=str,
        default=None,
        help="Output PNG filename (default: <field>.png).",
    )
    args = parser.parse_args()

    files = pick_input_files()

    # Build x-coordinates from the first file grid
    first_data, cell_dims, origin, spacing = load_and_stitch_data(files[0])
    nx = int(cell_dims[2])
    dx = float(spacing[2])
    x0 = float(origin[2])
    xe = x0 + np.arange(nx + 1) * dx
    x = 0.5 * (xe[:-1] + xe[1:])

    field = args.field
    if field not in first_data:
        print(
            f"Field '{field}' not found in {files[0]}. "
            f"Available: {list(first_data.keys())}",
            file=sys.stderr,
        )
        sys.exit(2)

    plt.clf()
    ax = plt.gca()

    colors = ["k", "r", "g"]
    for idx, fname in enumerate(files):
        stitched, _, _, _ = load_and_stitch_data(fname)
        if field not in stitched:
            print(
                f"Field '{field}' not in {fname}. Skipping.",
                file=sys.stderr,
            )
            continue

        y = stitched[field][:, 0, 0]
        color = colors[idx % len(colors)]
        t_val = parse_time_from_name(fname, args.time_scale)
        label = f"t={t_val:g}" if t_val is not None else os.path.basename(fname)
        ax.plot(x, y, label=label, linewidth=1.8, color=color)

    ax.set_xlabel(r"$\mathbf{x}$", fontweight="bold", labelpad=6)
    ax.set_ylabel(r"$\mathbf{E_c}$", fontweight="bold", labelpad=6)
    ax.set_xlim(x.min(), x.max())
    ax.set_ylim(0, 1)
    ax.legend()
    plt.tight_layout()

    outname = args.out if args.out is not None else f"{field}.png"
    plt.savefig(outname, dpi=720)
    print(f"Saved plot to {outname}")


if __name__ == "__main__":
    main()
