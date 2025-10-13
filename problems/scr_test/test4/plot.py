#!/usr/bin/env python3
"""
Stitch PIERNIK .h5 blocks in CWD and plot a 1D line-out vs. analytical solution.

By default:
  - searches for 'scr_tst_*.h5' (or any '*.h5' if none) in the current directory,
  - plots field 'escr_01' from z=0, y=0 for up to three files,
  - overlays an analytical profile using parameters sigma and vx.

Examples
--------
python plot.py
python plot.py --field escr_01 --vx 0.15 --sigma 10.0 --time-scale 0.2 --out plot.png
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
    # Prefer canonical first three if present
    preferred = ["scr_tst_0000.h5", "scr_tst_0001.h5", "scr_tst_0002.h5"]
    chosen = [f for f in preferred if f in files]
    if not chosen:
        chosen = files[:3]
    else:
        for f in files:
            if len(chosen) >= 3:
                break
            if f not in chosen:
                chosen.append(f)
    return chosen


def parse_time_from_name(fname: str, time_scale: float) -> float | None:
    """Parse time from 'scr_tst_0002.h5' -> 2 * time_scale; None if no match."""
    base = os.path.basename(fname)
    if base.startswith("scr_tst_") and base.endswith(".h5"):
        num = base.replace("scr_tst_", "").replace(".h5", "")
        try:
            return int(num) * time_scale
        except Exception:
            return None
    return None


def analytical_profile(x: np.ndarray, t: float, sigma: float, vx: float) -> np.ndarray:
    """
    y(t, x) = 1 / sqrt(1 + 160 t / (3 sigma)) *
              exp(-40 * (x - vx t)^2 / (1 + 160 t / (3 sigma))).
    """
    denom = 1.0 + (160.0 * t) / (3.0 * sigma)
    shift = x - vx * t
    return 1.0 / np.sqrt(denom) * np.exp(-40.0 * (shift * shift) / denom)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Stitch PIERNIK .h5 in CWD and plot 1D line-out vs analytical."
    )
    parser.add_argument(
        "--field", type=str, default="escr_01", help="Field to plot (default: escr_01)."
    )
    parser.add_argument(
        "--time-scale",
        type=float,
        default=0.2,
        help="Multiply file index by this to label time (default: 0.2).",
    )
    parser.add_argument(
        "--sigma",
        type=float,
        default=10.0,
        help="Sigma parameter used in the analytical solution (default: 10).",
    )
    parser.add_argument(
        "--vx",
        type=float,
        default=0.0,
        help="Gas velocity vx used in the analytical solution (default: 0.0).",
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

    ts: list[float] = []
    l1s: list[float] = []
    l2s: list[float] = []

    plt.clf()
    ax = plt.gca()

    # Colors for up to three snapshots
    colors_num = ["k", "r", "g"]
    colors_ana = ["r", "b", "c"]

    for idx, fname in enumerate(files[:3]):
        stitched, _, _, _ = load_and_stitch_data(fname)
        y = stitched[field][0, 0, :]

        t_val = parse_time_from_name(fname, args.time_scale)
        if t_val is None:
            t_val = float(idx) * args.time_scale

        y_ana = analytical_profile(x, t_val, sigma=args.sigma, vx=args.vx)

        ax.plot(x, y, label=f"t={t_val:g}", linewidth=2.0, color=colors_num[idx % 3])
        ax.plot(
            x,
            y_ana,
            label=f"t={t_val:g} (analytical)",
            linewidth=2.0,
            color=colors_ana[idx % 3],
            linestyle="dashed",
        )

        l1 = float(np.mean(np.abs(y - y_ana)))
        l2 = float(np.sqrt(np.mean((y - y_ana) ** 2)))
        ts.append(t_val)
        l1s.append(l1)
        l2s.append(l2)

    ax.set_xlabel(r"$\mathbf{x}$", fontweight="bold", labelpad=6)
    ax.set_ylabel(r"$\mathbf{E_c}$", fontweight="bold", labelpad=6)
    ax.set_xlim(x.min(), x.max())
    ax.set_ylim(0.0, 1.0)
    ax.legend()
    plt.tight_layout()

    # Table overlay
    hdr = f"{'t':>6} {'L1':>12} {'L2':>12}"
    rows = [f"{t:>6.3f} {l1:>12.3e} {l2:>12.3e}" for t, l1, l2 in zip(ts, l1s, l2s)]
    table_text = "\n".join([hdr] + rows)
    ax.text(
        0.01,
        0.98,
        table_text,
        transform=ax.transAxes,
        va="top",
        ha="left",
        family="monospace",
        fontsize=8,
        bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="black", lw=2),
    )

    outname = args.out if args.out is not None else f"{field}.png"
    plt.savefig(outname, dpi=720)
    print(f"Saved plot to {outname}")


if __name__ == "__main__":
    main()
