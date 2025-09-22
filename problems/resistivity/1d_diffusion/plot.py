#!/usr/bin/env python3
"""
Converts PIERNIK .h5 output in the current directory to a stitched 1D plot.
Looks for res_tst_*.h5 files (or any .h5 if those don't exist) in CWD,
stitches data blocks, compares with analytical solution and writes plot.png.

Usage example:
  python plot.py --eta 0.05 --a 40 --b0 0.1 --time-scale 0.1 --field mag_field_y --out plot.png
"""

import sys
import glob
import h5py
import numpy as np
import matplotlib.pyplot as plt
import os
import argparse


def get_field_names(h5_file):
    if "data" not in h5_file:
        raise ValueError("HDF5 file does not contain a '/data' group.")
    try:
        first_block_name = next(iter(h5_file["data"]))
    except StopIteration:
        raise ValueError("HDF5 '/data' group is empty.")
    field_names = list(h5_file["data"][first_block_name].keys())
    return field_names


def load_and_stitch_data(fname):
    """
    Load and stitch 3D cell-centered data from all blocks in the HDF5 file.

    Returns:
        stitched_data: dict[field_name] -> ndarray shaped (Nz, Ny, Nx)
        global_cell_dims: ndarray (Nx, Ny, Nz) as integers
        origin: ndarray (x0, y0, z0) (left edge)
        spacing: ndarray (dx, dy, dz)
    """
    with h5py.File(fname, "r") as f:
        field_names = get_field_names(f)

        block_names = list(f["data"].keys())
        all_offsets = np.array([f["data"][name].attrs["off"] for name in block_names], dtype=int)
        all_dims = np.array([f["data"][name].attrs["n_b"] for name in block_names], dtype=int)

        # Determine global cell counts (cells per axis)
        global_cell_dims = np.max(all_offsets + all_dims, axis=0).astype(int)

        # Physical grid from simulation parameters
        origin = np.array(f["simulation_parameters"].attrs["domain_left_edge"], dtype=float)
        domain_right = np.array(f["simulation_parameters"].attrs["domain_right_edge"], dtype=float)
        domain_size = domain_right - origin

        # Spacing = domain_size / number_of_cells (avoid divide by zero)
        spacing = domain_size / np.maximum(1.0, global_cell_dims.astype(float))

        stitched_data = {}
        # On disk arrays are (z,y,x) according to the user's reference.
        for field in field_names:
            stitched_field = np.empty((global_cell_dims[2],
                                       global_cell_dims[1],
                                       global_cell_dims[0]), dtype=np.float32)

            for i, name in enumerate(block_names):
                block_data = f["data"][name][field][:]
                offset = all_offsets[i]
                dims = all_dims[i]

                slc = np.s_[offset[2]: offset[2] + dims[2],
                            offset[1]: offset[1] + dims[1],
                            offset[0]: offset[0] + dims[0]]
                stitched_field[slc] = block_data

            stitched_data[field] = stitched_field

    return stitched_data, global_cell_dims, origin, spacing


# Preferred rc parameters
plt.rcParams.update({
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
})


def y_analytical(t, x, b0, a, eta):
    """
    Analytical solution used in the original script:
      y(x,t) = b0 / sqrt(1 + 4 a eta t) * exp( - a x^2 / (1 + 4 a eta t) )
    """
    denom = 1.0 + 4.0 * a * eta * t
    pref = b0 / np.sqrt(denom)
    return pref * np.exp(-(a * x ** 2) / denom)


def pick_input_files():
    """
    Choose which .h5 files to plot from current directory.
    Priority:
      1) files matching res_tst_*.h5 (sorted)
      2) any .h5 files (sorted)
    Returns up to 3 files (ordered) to reproduce 0.0,0.3,0.5 convention if present;
    otherwise returns up to three earliest .h5 files.
    """
    candidates = sorted(glob.glob("res_tst_*.h5"))
    if not candidates:
        candidates = sorted(glob.glob("*.h5"))
    if not candidates:
        print("No .h5 files found in current directory.", file=sys.stderr)
        sys.exit(1)

    # Prefer canonical names if present
    chosen = []
    want = ["res_tst_0000.h5", "res_tst_0003.h5", "res_tst_0005.h5"]
    for w in want:
        if w in candidates:
            chosen.append(w)
    if not chosen:
        chosen = candidates[:3]
    else:
        for c in candidates:
            if len(chosen) >= 3:
                break
            if c not in chosen:
                chosen.append(c)
    return chosen


def parse_time_from_name(fname, time_scale):
    """
    Parse time from filename like res_tst_0003.h5 -> int(0003) * time_scale
    Default time_scale is 0.1 so 0003 -> 3 * 0.1 = 0.3
    """
    base = os.path.basename(fname)
    if base.startswith("res_tst_") and base.endswith(".h5"):
        num = base.replace("res_tst_", "").replace(".h5", "")
        try:
            return int(num) * time_scale
        except Exception:
            return None
    return None


def main():
    parser = argparse.ArgumentParser(description="Stitch PIERNIK .h5 outputs and plot 1D profile vs analytical.")
    parser.add_argument("--eta", type=float, default=0.05, help="Diffusivity eta (default: 0.05)")
    parser.add_argument("--a", type=float, default=40.0, help="Parameter a in analytical solution (default: 40)")
    parser.add_argument("--b0", type=float, default=0.1, help="Initial amplitude b0 (default: 0.1)")
    parser.add_argument("--out", type=str, default="plot.png", help="Output image filename (default: plot.png)")
    args = parser.parse_args()

    files = pick_input_files()
    data_to_plot = 'mag_field_y'
    b0 = float(args.b0)
    a = float(args.a)
    eta = float(args.eta)
    ts = []
    L1s = []
    L2s = []

    # Build x coordinates using the first file's grid
    first_data, first_cell_dims, first_origin, first_spacing = load_and_stitch_data(files[0])
    N = int(first_cell_dims[0])
    dx = float(first_spacing[0])
    x0 = float(first_origin[0])
    xe = x0 + np.arange(N + 1) * dx
    x = 0.5 * (xe[:-1] + xe[1:])

    plt.clf()
    ax = plt.gca()

    # requested color sequences:
    colors_solid = ["r", "g", "b"]   # rgb for solid (numerical)
    colors_dashed = ["cyan", "magenta", "red"]  # gbr for dashed (analytical)

    for idx, fname in enumerate(files):
        stitched, cell_dims, origin, spacing = load_and_stitch_data(fname)
        if data_to_plot not in stitched:
            print(f"Field '{data_to_plot}' not found in {fname}. Available fields: {list(stitched.keys())}",
                  file=sys.stderr)
            continue

        # select 1D slice: (z=0, y=0, :)
        y = stitched[data_to_plot][0, 0, :]

        # choose colors cycling through the provided sequences
        color_solid = colors_solid[idx % len(colors_solid)]
        color_dashed = colors_dashed[idx % len(colors_dashed)]

        # Derive time value and formatted label
        t_val = parse_time_from_name(fname, 0.1)
        if t_val is None:
            t_val = float(idx)
        time_label = f"t={t_val:g}"  # uses "g" to show 0, 0.3, 0.5 etc.

        # Plot numerical (solid) and analytical (dashed) with distinct colors
        numeric_label = f"{time_label} (numerical)"
        analytic_label = f"{time_label} (analytical)"

        ax.plot(x, y, label=numeric_label, linewidth=2.0, color=color_solid, linestyle="solid")
        ax.plot(x, y_analytical(t_val, x, b0=b0, a=a, eta=eta),
                label=analytic_label, linewidth=2.0, color=color_dashed, linestyle="dashed")

        L1 = np.mean(np.abs(y - y_analytical(t_val, x, b0=b0, a=a, eta=eta)))
        L2 = np.sqrt(np.mean((y - y_analytical(t_val, x, b0=b0, a=a, eta=eta)) ** 2))
        ts.append(t_val)
        L1s.append(L1)
        L2s.append(L2)

    ax.set_xlabel(r"$\mathbf{x}$", fontweight="bold", labelpad=6)
    ax.set_ylabel(r"$\mathbf{B_y}$", fontweight="bold", labelpad=6)
    ax.set_xlim(-2.0, 2.0)
    ax.set_ylim(0, 0.1)
    ax.set_xticks([-2, -1, 0, 1, 2])
    ax.legend()
    plt.tight_layout()

    # Compose table text
    hdr = f"{'t':>6} {'L1':>12} {'L2':>12}"
    rows = [f"{t:>6.3f} {l1:>12.3e} {l2:>12.3e}" for t, l1, l2 in zip(ts, L1s, L2s)]
    table_text = "\n".join([hdr] + rows)
    plt.text(
        0.01, 0.98, table_text,
        transform=ax.transAxes,
        va="top", ha="left",
        family="monospace",
        fontsize=8, bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="black", lw=0.5))

    outname = args.out
    plt.savefig(outname, dpi=720)
    print(f"Saved plot to {outname}")


if __name__ == "__main__":
    main()
