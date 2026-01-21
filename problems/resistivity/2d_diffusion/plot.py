#!/usr/bin/env python3
"""
PIERNIK 2D resistive sinusoidal B-field test:
Make a 2x2 figure in the xy-plane:
  (1) Numerical field (mag_field_x or mag_field_y)
  (2) Analytical field
  (3) Error = numerical - analytical
  (4) |Error| (useful dynamic range)

Analytical (v=0, uniform eta):
  Bx(x,y,t) = B0*ky*sin(kx x)*cos(ky y)*exp(-eta*(kx^2+ky^2)*t)
  By(x,y,t) = -B0*kx*cos(kx x)*sin(ky y)*exp(-eta*(kx^2+ky^2)*t)
Usage:
  python res2d_map.py --file res_tst_0005.h5 --eta 0.01 --field mag_field_x --out Bx_2d.png
  python res2d_map.py --eta 0.01 --field mag_field_y --time 0.5 --out By_2d.png
"""
import os
import glob
import argparse
import numpy as np
import h5py
import matplotlib.pyplot as plt


def get_field_names(h5_file):
    if "data" not in h5_file:
        raise ValueError("HDF5 file does not contain a '/data' group.")
    first_block_name = next(iter(h5_file["data"]))
    return list(h5_file["data"][first_block_name].keys())


def pick_input_file(pattern="res_tst_*.h5"):
    cands = sorted(glob.glob(pattern))
    if not cands:
        cands = sorted(glob.glob("*.h5"))
    if not cands:
        raise FileNotFoundError("No .h5 files found in current directory.")
    return cands[-1]


def load_and_stitch_data(fname):
    """
    Returns:
        stitched_data: dict[field] -> ndarray (Nz,Ny,Nx)
        global_cell_dims: (Nx,Ny,Nz) ints
        origin: (x0,y0,z0) left edge
        spacing: (dx,dy,dz)
        sim_par_attrs: dict of simulation_parameters attrs
    """
    with h5py.File(fname, "r") as f:
        field_names = get_field_names(f)
        block_names = list(f["data"].keys())

        all_offsets = np.array([f["data"][bn].attrs["off"] for bn in block_names], dtype=int)
        all_dims = np.array([f["data"][bn].attrs["n_b"] for bn in block_names], dtype=int)

        global_cell_dims = np.max(all_offsets + all_dims, axis=0).astype(int)  # (Nx,Ny,Nz)

        sp = f["simulation_parameters"].attrs
        origin = np.array(sp["domain_left_edge"], dtype=float)
        right = np.array(sp["domain_right_edge"], dtype=float)
        domain_size = right - origin
        spacing = domain_size / np.maximum(1.0, global_cell_dims.astype(float))

        sim_par_attrs = {k: sp[k] for k in sp.keys()}

        stitched_data = {}
        for field in field_names:
            arr = np.empty((global_cell_dims[2], global_cell_dims[1], global_cell_dims[0]), dtype=np.float32)  # (Nz,Ny,Nx)

            for off, dims, bn in zip(all_offsets, all_dims, block_names):
                block = f["data"][bn][field][:]
                slc = np.s_[off[2]:off[2] + dims[2],
                            off[1]:off[1] + dims[1],
                            off[0]:off[0] + dims[0]]
                arr[slc] = block

            stitched_data[field] = arr

    return stitched_data, global_cell_dims, origin, spacing, sim_par_attrs


def try_get_time(sim_par_attrs):
    for key in ["time", "t", "sim_time", "current_time"]:
        if key in sim_par_attrs:
            try:
                return float(sim_par_attrs[key])
            except Exception:
                pass
    return None


def parse_time_from_name(fname, time_scale=0.1):
    base = os.path.basename(fname)
    if base.startswith("res_tst_") and base.endswith(".h5"):
        num = base.replace("res_tst_", "").replace(".h5", "")
        try:
            return int(num) * time_scale
        except Exception:
            return None
    return None


def analytic_field_xy(xc, yc, t, B0, kx, ky, eta, which):
    """
    xc, yc are 2D arrays of cell centers (Ny,Nx)
    which: 'x' or 'y'
    """
    k2 = kx * kx + ky * ky
    decay = np.exp(-eta * k2 * t)

    if which == "x":
        return (B0 * ky * np.sin(kx * xc) * np.cos(ky * yc)) * decay
    elif which == "y":
        return (-B0 * kx * np.cos(kx * xc) * np.sin(ky * yc)) * decay
    else:
        raise ValueError("which must be 'x' or 'y'")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--file", type=str, default=None, help="Input .h5 file (default: pick latest res_tst_*.h5)")
    ap.add_argument("--eta", type=float, required=True, help="Resistivity eta used in the run")
    ap.add_argument("--B0", type=float, default=0.01, help="B0 used in IC (default 0.01)")
    ap.add_argument("--kx", type=float, default=6.283185307, help="kx (default 2pi)")
    ap.add_argument("--ky", type=float, default=6.283185307, help="ky (default 2pi)")
    ap.add_argument("--field", type=str, default="mag_field_x", help="mag_field_x or mag_field_y")
    ap.add_argument("--z-index", type=int, default=0, help="z-slab index (default 0)")
    ap.add_argument("--time", type=float, default=None, help="Override time t explicitly")
    ap.add_argument("--time-scale", type=float, default=0.1, help="If parsing from res_tst_####.h5, t=####*time_scale")
    ap.add_argument("--out", type=str, default="B_2d_compare.png", help="Output figure filename")
    ap.add_argument("--relative", action="store_true", help="Plot relative error panels instead of absolute")
    args = ap.parse_args()

    fname = args.file or pick_input_file()
    stitched, cell_dims, origin, spacing, sim_attrs = load_and_stitch_data(fname)

    if args.field not in stitched:
        raise KeyError(f"Field '{args.field}' not found. Available: {list(stitched.keys())}")

    Nx, Ny, Nz = int(cell_dims[0]), int(cell_dims[1]), int(cell_dims[2])
    dx, dy = float(spacing[0]), float(spacing[1])
    x0, y0 = float(origin[0]), float(origin[1])

    zi = int(args.z_index)
    if not (0 <= zi < Nz):
        raise ValueError(f"z-index out of range : {zi} not in [0,{Nz - 1}]")

    # Determine time
    t = args.time
    if t is None:
        t = try_get_time(sim_attrs)
    if t is None:
        t = parse_time_from_name(fname, time_scale=args.time_scale)
    if t is None:
        raise RuntimeError("Could not determine time from file; pass --time explicitly.")

    # Cell-center grids (Ny,Nx)
    x = x0 + (np.arange(Nx) + 0.5) * dx
    y = y0 + (np.arange(Ny) + 0.5) * dy
    xc, yc = np.meshgrid(x, y, indexing="xy")  # (Ny,Nx)

    # Numerical 2D slice
    num = stitched[args.field][zi, :, :].astype(np.float64)  # (Ny,Nx)

    # Analytical
    which = "x" if args.field.endswith("_x") else "y"
    ana = analytic_field_xy(xc, yc, t, args.B0, args.kx, args.ky, args.eta, which=which)

    # Errors
    err = num - ana
    tiny = 1e-30
    if args.relative:
        denom = np.maximum(np.abs(ana), tiny)
        err_plot = err / denom
        abs_err_plot = np.abs(err) / denom
        err_label = "relative error (num-ana)/|ana|"
        abs_err_label = "relative |error|"
    else:
        err_plot = err
        abs_err_plot = np.abs(err)
        err_label = "error (num - ana)"
        abs_err_label = "|error|"

    # Extent for imshow in physical coordinates
    # Use cell-center extents (approx). For better fidelity, use left/right edges:
    xL, xR = x0, x0 + Nx * dx
    yL, yR = y0, y0 + Ny * dy
    extent = [xL, xR, yL, yR]

    # Norms (over the plane)
    L1 = np.mean(np.abs(err))
    L2 = np.sqrt(np.mean(err * err))
    Linf = np.max(np.abs(err))
    print(f"file = {fname}")
    print(f"field = {args.field}, z-index = {zi}, t = {t:g}")
    print(f"Errors on full xy plane: L1={L1:.6e}  L2={L2:.6e}  Linf={Linf:.6e}")

    # Plot
    plt.rcParams.update({
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "axes.linewidth": 1.2,
        "font.size": 11,
        "figure.figsize": (10, 8),
    })

    fig, axs = plt.subplots(2, 2, constrained_layout=True)

    # Use shared vmin/vmax for num and ana for meaningful comparison
    vmin = min(np.min(num), np.min(ana))
    vmax = max(np.max(num), np.max(ana))

    im0 = axs[0, 0].imshow(num, origin="lower", extent=extent, aspect="auto", vmin=vmin, vmax=vmax)
    axs[0, 0].set_title(f"Numerical {args.field}")
    axs[0, 0].set_xlabel("x")
    axs[0, 0].set_ylabel("y")
    fig.colorbar(im0, ax=axs[0, 0], fraction=0.046, pad=0.04)

    im1 = axs[0, 1].imshow(ana, origin="lower", extent=extent, aspect="auto", vmin=vmin, vmax=vmax)
    axs[0, 1].set_title("Analytical")
    axs[0, 1].set_xlabel("x")
    axs[0, 1].set_ylabel("y")
    fig.colorbar(im1, ax=axs[0, 1], fraction=0.046, pad=0.04)

    # Center error colormap around 0 for signed error
    emax = np.max(np.abs(err_plot))
    im2 = axs[1, 0].imshow(err_plot, origin="lower", extent=extent, aspect="auto", vmin=-emax, vmax=emax)
    axs[1, 0].set_title(err_label)
    axs[1, 0].set_xlabel("x")
    axs[1, 0].set_ylabel("y")
    fig.colorbar(im2, ax=axs[1, 0], fraction=0.046, pad=0.04)

    im3 = axs[1, 1].imshow(abs_err_plot, origin="lower", extent=extent, aspect="auto")
    axs[1, 1].set_title(abs_err_label)
    axs[1, 1].set_xlabel("x")
    axs[1, 1].set_ylabel("y")
    fig.colorbar(im3, ax=axs[1, 1], fraction=0.046, pad=0.04)

    # Annotate norms
    txt = f"t={t:g}\nL1={L1:.3e}\nL2={L2:.3e}\nLinf={Linf:.3e}"
    axs[0, 0].text(0.02, 0.98, txt, transform=axs[0, 0].transAxes,
                   va="top", ha="left", family="monospace",
                   bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="black", lw=0.8))

    fig.suptitle(f"Resistive decay test: {args.field} (eta={args.eta:g}, kx={args.kx:g}, ky={args.ky:g})", y=1.02)
    plt.savefig(args.out)
    print(f"Saved figure to {args.out}")


if __name__ == "__main__":
    main()
