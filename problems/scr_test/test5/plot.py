#!/usr/bin/env python3
"""
Stitch PIERNIK .h5 blocks found in the **current directory** and make a 2x2 figure:
  (0,0) numerical at t0
  (0,1) numerical at t1
  (1,0) analytical reference at t1
  (1,1) percent error |num-ana|/ana at t1

It prefers files 'scr_tst_0000.h5' and 'scr_tst_0001.h5'. If those are absent,
it falls back to the first two '*.h5' files (sorted). The figure is saved as
'escr_01.png' (or change `data_to_plot` below).
"""

from __future__ import annotations

import glob
import os
import sys
from typing import Tuple, List

import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import special as sp


# ------------------------------- I/O utilities -------------------------------

def get_field_names(h5_file: h5py.File) -> List[str]:
    if "data" not in h5_file:
        raise ValueError("HDF5 file is missing group '/data'.")
    try:
        first_block = next(iter(h5_file["data"]))
    except StopIteration as exc:
        raise ValueError("HDF5 group '/data' is empty.") from exc
    return list(h5_file["data"][first_block].keys())


def load_and_stitch_data(fname: str):
    """
    Load and stitch 3D fields from all blocks.

    Returns
    -------
    stitched_data : dict[name] -> ndarray (Nz, Ny, Nx)
    global_cell_dims : np.ndarray [Nx, Ny, Nz]
    origin : np.ndarray [x0, y0, z0]
    spacing : np.ndarray [dx, dy, dz]
    """
    stitched_data = {}
    with h5py.File(fname, "r") as f:
        field_names = get_field_names(f)

    with h5py.File(fname, "r") as f:
        block_names = list(f["data"].keys())
        all_offsets = np.array([f["data"][n].attrs["off"] for n in block_names], dtype=int)
        all_dims = np.array([f["data"][n].attrs["n_b"] for n in block_names], dtype=int)

        global_cell_dims = np.max(all_offsets + all_dims, axis=0)

        origin = f["simulation_parameters"].attrs["domain_left_edge"]
        domain_size = f["simulation_parameters"].attrs["domain_right_edge"] - origin
        spacing = domain_size / np.maximum(1, global_cell_dims)

        for field in field_names:
            vol = np.empty(
                (global_cell_dims[2], global_cell_dims[1], global_cell_dims[0]),
                dtype=np.float32,
            )
            for i, name in enumerate(block_names):
                arr = f["data"][name][field][:]
                off = all_offsets[i]
                dims = all_dims[i]
                vol[
                    off[2]: off[2] + dims[2],
                    off[1]: off[1] + dims[1],
                    off[0]: off[0] + dims[0],
                ] = arr
            stitched_data[field] = vol

    return stitched_data, global_cell_dims.astype(int), origin, spacing


def pick_two_files() -> Tuple[str, str]:
    """Pick two files from CWD: prefer scr_tst_0000.h5 and scr_tst_0001.h5, else first two *.h5."""
    preferred = ["scr_tst_0000.h5", "scr_tst_0001.h5"]
    cand = sorted(glob.glob("scr_tst_*.h5")) or sorted(glob.glob("*.h5"))
    cand = [c for c in cand if c.lower().endswith(".h5")]
    if len(cand) < 2:
        print("Need at least two .h5 files in the current directory.", file=sys.stderr)
        sys.exit(1)

    chosen = [p for p in preferred if p in cand]
    if len(chosen) >= 2:
        return chosen[0], chosen[1]

    # Otherwise pick the first two available
    return cand[0], cand[1]


# ------------------------------- Plot helpers --------------------------------

def _format_axes(ax: plt.Axes) -> None:
    ax.set_xlim(-1.0, 1.0)
    ax.set_ylim(-1.0, 1.0)
    ax.set_xticks([-0.5, 0.0, 0.5])
    ax.set_yticks([-0.5, 0.0, 0.5])


# ---------------------------------- Main -------------------------------------

def main() -> None:
    plt.rcParams.update({
        "figure.dpi": 150, "savefig.dpi": 720,
        "axes.linewidth": 2.5, "font.size": 12,
        "xtick.major.size": 5, "ytick.major.size": 5,
        "xtick.major.width": 1.0, "ytick.major.width": 1.0,
        "xtick.direction": "out", "ytick.direction": "out",
        "legend.frameon": True, "legend.edgecolor": "black", "legend.framealpha": 1.0,
        "figure.figsize": (8, 6),
    })

    data_to_plot = "escr_01"

    f0, f1 = pick_two_files()
    print(f"Using files: {f0} (t0), {f1} (t1)")

    # Build grid from the first file found in CWD
    data0, cell_dims, origin, spacing = load_and_stitch_data(f0)

    nx, dx, x0 = cell_dims[0], spacing[0], origin[0]
    ny, dy, y0 = cell_dims[1], spacing[1], origin[1]
    xe = x0 + np.arange(nx + 1) * dx
    ye = y0 + np.arange(ny + 1) * dy
    x = 0.5 * (xe[:-1] + xe[1:])
    y = 0.5 * (ye[:-1] + ye[1:])
    Xc, Yc = np.meshgrid(x, y, indexing="xy")

    fig, ax = plt.subplots(2, 2, figsize=(8, 6))

    # ---------------- Panel (0,0): t = 0.0
    if data_to_plot not in data0 or "mag_field_x" not in data0 or "mag_field_y" not in data0:
        print(f"{f0} missing fields for plotting; available: {list(data0.keys())}", file=sys.stderr)
        sys.exit(2)

    d0 = data0[data_to_plot][0, :, :]
    im = ax[0, 0].imshow(d0, extent=[x[0], x[-1], y[0], y[-1]], origin="lower", cmap="viridis")
    bx0 = np.nan_to_num(data0["mag_field_x"][0, :, :], nan=0.0, posinf=0.0, neginf=0.0)
    by0 = np.nan_to_num(data0["mag_field_y"][0, :, :], nan=0.0, posinf=0.0, neginf=0.0)
    ax[0, 0].streamplot(x, y, bx0, by0, density=0.4, color="k", maxlength=20.0, integration_direction="both")
    ax[0, 0].set_xlabel(r"$x$", fontweight="bold")
    ax[0, 0].set_ylabel(r"$y$", fontweight="bold")
    ax[0, 0].set_title("t=0.0", fontweight="bold")
    _format_axes(ax[0, 0])
    fig.colorbar(im, ax=ax[0, 0], fraction=0.046, pad=0.04).set_label(r"$E_c$", fontweight="bold", labelpad=6)

    # ---------------- Panel (0,1): t = 1st later snapshot
    data1, _, _, _ = load_and_stitch_data(f1)
    if data_to_plot not in data1 or "mag_field_x" not in data1 or "mag_field_y" not in data1:
        print(f"{f1} missing fields for plotting; available: {list(data1.keys())}", file=sys.stderr)
        sys.exit(2)

    d1 = data1[data_to_plot][0, :, :]
    im = ax[0, 1].imshow(d1, extent=[x[0], x[-1], y[0], y[-1]], origin="lower", cmap="viridis")
    bx1 = np.nan_to_num(data1["mag_field_x"][0, :, :], nan=0.0, posinf=0.0, neginf=0.0)
    by1 = np.nan_to_num(data1["mag_field_y"][0, :, :], nan=0.0, posinf=0.0, neginf=0.0)
    ax[0, 1].streamplot(x, y, bx1, by1, density=0.8, color="k", maxlength=10.0, integration_direction="both")
    ax[0, 1].set_xlabel(r"$x$", fontweight="bold")
    ax[0, 1].set_ylabel(r"$y$", fontweight="bold")
    ax[0, 1].set_title("t=Δt", fontweight="bold")
    _format_axes(ax[0, 1])
    fig.colorbar(im, ax=ax[0, 1], fraction=0.046, pad=0.04).set_label(r"$E_c$", fontweight="bold", labelpad=6)

    # ---------------- Panel (1,0): analytical reference at the 2nd time
    # Parameters for the ring-like analytical model (adjust as needed)
    phi0 = np.pi / 12.0
    sigma_par = 1.0
    t = 0.26   # if you want, infer from filename index×scale
    D = np.sqrt(4.0 * t / (3.0 * sigma_par))
    rin, rout = 0.5, 0.7
    r = np.hypot(Xc, Yc)
    phi = np.arctan2(Yc, Xc)
    Ec = np.full_like(r, 10.0, dtype=float)
    mask = (r > rin) & (r < rout)
    arg1 = ((phi - phi0) * r) / D
    arg2 = ((phi + phi0) * r) / D
    Ec[mask] += sp.erfc(arg1[mask]) - sp.erfc(arg2[mask])

    im = ax[1, 0].imshow(Ec, extent=[x[0], x[-1], y[0], y[-1]], origin="lower", cmap="viridis")
    ax[1, 0].set_xlabel(r"$x$", fontweight="bold")
    ax[1, 0].set_ylabel(r"$y$", fontweight="bold")
    ax[1, 0].set_title("analytical (t=Δt)", fontweight="bold")
    _format_axes(ax[1, 0])
    fig.colorbar(im, ax=ax[1, 0], fraction=0.046, pad=0.04).set_label(r"$E_c$", fontweight="bold", labelpad=6)

    # ---------------- Panel (1,1): percent error vs analytical
    l1 = float(np.mean(np.abs(Ec.ravel() - d1.ravel())))
    l2 = float(np.sqrt(np.mean((Ec.ravel() - d1.ravel()) ** 2)))
    err = 100.0 * np.abs(Ec - d1) / np.maximum(Ec, 1e-30)

    im = ax[1, 1].imshow(err, extent=[x[0], x[-1], y[0], y[-1]], origin="lower", cmap="viridis")
    ax[1, 1].set_xlabel(r"$x$", fontweight="bold")
    ax[1, 1].set_ylabel(r"$y$", fontweight="bold")
    ax[1, 1].set_title(fr"% error, $L_1={l1:.3e}$, $L_2={l2:.3e}$", fontweight="bold", fontsize=11)
    _format_axes(ax[1, 1])
    fig.colorbar(im, ax=ax[1, 1], fraction=0.046, pad=0.04).set_label(r"$E_c$", fontweight="bold", labelpad=6)

    plt.tight_layout()
    for a in (ax[0, 0], ax[0, 1], ax[1, 0], ax[1, 1]):
        for s in a.spines.values():
            s.set_linewidth(3)

    fig.add_artist(
        mpatches.Rectangle(
            (0.005, 0.005), 0.99, 0.99, transform=fig.transFigure,
            fill=False, lw=4, ec="black"
        )
    )

    outname = "escr_01.png"
    plt.savefig(outname, dpi=720)
    print(f"Saved plot to ./{outname}")


if __name__ == "__main__":
    main()
