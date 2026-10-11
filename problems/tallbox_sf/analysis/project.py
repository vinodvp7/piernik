#!/usr/bin/env python3
"""
Projections of tallbox_sf snapshots along z (xy plane, face-on) or x / y (yz / xz plane, edge-on).

Reads the uniform-grid HDF5 outputs directly (PSM units: pc, Msun, Myr) and writes one multi-panel image per snapshot.

Usage:

    python project.py SNAP.h5 [SNAP.h5 ...] | RUNDIR  [--axis z|x|y] [--fields dens temp pres velz ...]
                      [--zrange ZMIN ZMAX] [--clim dens=1e-2,1e2 ...] [--every N] [--no-stars] [-o OUTDIR]

e.g. from piernik/runs:
    python ../problems/tallbox_sf/analysis/project.py tallbox_sf_A_20pc --axis x --fields dens temp pres vlos
    python ../problems/tallbox_sf/analysis/project.py tallbox_sf_A_20pc/tallbox_sf_A_h_0010.h5 --axis z --zrange -300 300

Fields (aliases, default projection in brackets):
    dens   gas surface density Sigma [sum, Msun/pc^2]
    temp   temperature [mass-weighted mean, K]
    pres   thermal pressure P/k_B [volume-weighted mean, K cm^-3]
    nH     hydrogen number density [mass-weighted mean, cm^-3]
    velx, vely, velz, vlos (velocity along the line of sight) [mass-weighted mean, km/s]
    sne    SN count in the column [sum]       sfr   SFR_n in the column [sum]
Any other dataset of the snapshot (e.g. magnetic or CR fields) can be given by its HDF5 name. The projection can be
overridden with NAME:MODE, MODE one of sum (column integral), mean (volume-weighted), mwmean (mass-weighted), max,
min, slice (the plane through the origin), e.g.  temp:max  dens:slice.
Star particles: formation_time is -10 while a sink accretes its first mass_SN*n_SN and is reset to t at every new
loadout (which starts the delayed-SN clock), so it is not a birth time. Particles are drawn as accreting sinks
(white), loadout younger than --young (cyan) and older (grey); Msun in stars is printed in the title.
--zrange restricts z to [ZMIN, ZMAX] pc: the plotted height for x/y projections, the integrated layer for z.
"""

import argparse
import glob
import os

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
import numpy as np

# PSM code units -> physical (same constants as compare.py)
NH_PER_RHO = 6.770e-23 / (1.4 * 1.6726e-24)      # 1 Msun/pc^3 -> n_H [cm^-3] (1.4 m_H per H)
POK_PER_P = 6.47e-13 / 1.3807e-16                 # 1 Msun/(pc Myr^2) -> P/k_B [K cm^-3]
KMS_PER_V = 0.9778                                # 1 pc/Myr -> km/s

AXES = {"x": 2, "y": 1, "z": 0}                   # array axis of each direction (datasets are stored as [z, y, x])

# alias -> (dataset(s), default mode, factor, label, colour map, log scale, default colour limits)
FIELDS = {
    "dens": ("density", "sum", 1.0, r"$\Sigma_\mathrm{gas}$ [M$_\odot$ pc$^{-2}$]", "viridis", True, (1e-3, 1e2)),
    "nH":   ("density", "mwmean", NH_PER_RHO, r"$\langle n_\mathrm{H} \rangle$ [cm$^{-3}$]", "viridis", True, (1e-4, 1e2)),
    "temp": (("temperature", "Temperature"), "mwmean", 1.0, r"$\langle T \rangle$ [K]", "inferno", True, (1e2, 1e7)),
    "pres": ("pressure", "mean", POK_PER_P, r"$\langle P/k_\mathrm{B} \rangle$ [K cm$^{-3}$]", "magma", True, (1e-1, 1e5)),
    "velx": ("velocity_x", "mwmean", KMS_PER_V, r"$\langle v_x \rangle$ [km/s]", "RdBu_r", False, None),
    "vely": ("velocity_y", "mwmean", KMS_PER_V, r"$\langle v_y \rangle$ [km/s]", "RdBu_r", False, None),
    "velz": ("velocity_z", "mwmean", KMS_PER_V, r"$\langle v_z \rangle$ [km/s]", "RdBu_r", False, None),
    "sne":  ("SNe_n", "sum", 1.0, "SNe in column", "cividis", False, None),
    "sfr":  ("SFR_n", "sum", 1.0, "SFR_n in column", "cividis", False, None),
}
VLOS = {"x": "velx", "y": "vely", "z": "velz"}
MODES = ("sum", "mean", "mwmean", "max", "min", "slice")


# ----------------------------------------------------------------------------- reading

def read_snapshot(fn, names):
    """Full-domain arrays [z, y, x] of the requested datasets, plus geometry, time and star particles."""
    with h5py.File(fn, "r") as f:
        sp = f["simulation_parameters"].attrs
        nx, ny, nz = (int(v) for v in sp["domain_dimensions"])
        left = np.array(sp["domain_left_edge"], dtype=float)
        right = np.array(sp["domain_right_edge"], dtype=float)
        if np.any(f["grid_level"][:] != 0):
            raise SystemExit(f"{fn}: refined levels found, only uniform grids are supported")
        out = {n: np.empty((nz, ny, nx), dtype=np.float32) for n in names}
        stars = {k: [] for k in ("x", "y", "z", "mass", "tform")}
        gli = f["grid_left_index"][:]
        gdim = f["grid_dimensions"][:]
        for g, (i0, j0, k0), (di, dj, dk) in zip(sorted(f["data"].keys()), gli, gdim):
            grp = f["data"][g]
            for n in names:
                out[n][k0:k0 + dk, j0:j0 + dj, i0:i0 + di] = grp[n][:]
            if "particles/stars" in grp and grp["particles/stars/mass"].shape[0] > 0:
                ps = grp["particles/stars"]
                for k, d in (("x", "position_x"), ("y", "position_y"), ("z", "position_z"), ("mass", "mass"),
                             ("tform", "formation_time")):
                    stars[k].append(ps[d][:])
        time = float(f.attrs["time"][0])
    stars = {k: (np.concatenate(v) if v else np.empty(0)) for k, v in stars.items()}
    dl = (right - left) / np.array([nx, ny, nz])
    return out, left, right, dl, time, stars


def available_datasets(fn):
    with h5py.File(fn, "r") as f:
        return set(f["data"][sorted(f["data"].keys())[0]].keys())


def resolve_fields(specs, axis, avail):
    """List of (title, dataset, mode, factor, label, cmap, log, clim) for the requested field specs."""
    out = []
    for spec in specs:
        name, _, mode = spec.partition(":")
        alias = VLOS[axis] if name == "vlos" else name
        if alias in FIELDS:
            ds, dmode, fac, label, cmap, log, clim = FIELDS[alias]
            if isinstance(ds, tuple):
                ds = next((d for d in ds if d in avail), ds[0])
        else:
            ds, dmode, fac, label, cmap, log, clim = name, "mwmean", 1.0, name, "viridis", None, None
        if ds not in avail:
            raise SystemExit(f"field '{name}': dataset '{ds}' not in the snapshot (has: {', '.join(sorted(avail))})")
        mode = mode or dmode
        if mode not in MODES:
            raise SystemExit(f"field '{spec}': unknown projection '{mode}' (use one of {', '.join(MODES)})")
        if mode != dmode and alias in FIELDS:
            label = f"{label}  ({mode})"
            if mode == "sum" or dmode == "sum":
                clim = None                        # the default limits belong to the default projection
        out.append((spec, ds, mode, fac, label, cmap, log, clim))
    return out


# ----------------------------------------------------------------------------- projection

def project(arr, rho, mode, ax, dl):
    if mode == "sum":
        return arr.sum(axis=ax, dtype=np.float64) * dl
    if mode == "mean":
        return arr.mean(axis=ax, dtype=np.float64)
    if mode == "mwmean":
        return (arr * rho).sum(axis=ax, dtype=np.float64) / rho.sum(axis=ax, dtype=np.float64)
    if mode == "max":
        return arr.max(axis=ax)
    if mode == "min":
        return arr.min(axis=ax)
    if mode == "slice":
        return np.take(arr, arr.shape[ax] // 2, axis=ax)
    raise ValueError(mode)


def plot_snapshot(fn, args, outdir):
    avail = available_datasets(fn)
    fields = resolve_fields(args.fields, args.axis, avail)
    names = sorted({f[1] for f in fields} | {"density"})
    data, left, right, dl, time, stars = read_snapshot(fn, names)

    # restrict z
    zc = left[2] + (np.arange(data["density"].shape[0]) + 0.5) * dl[2]
    zlo, zhi = (args.zrange if args.zrange else (left[2], right[2]))
    ksel = (zc >= zlo) & (zc <= zhi)
    if not ksel.any():
        raise SystemExit(f"--zrange {zlo:g} {zhi:g} contains no cells")
    data = {k: v[ksel] for k, v in data.items()}
    zlo_e, zhi_e = zc[ksel][0] - dl[2] / 2, zc[ksel][-1] + dl[2] / 2

    ax = AXES[args.axis]
    lmb = {"x": dl[0], "y": dl[1], "z": dl[2]}[args.axis]
    if args.axis == "z":                     # image [y, x]
        extent = (left[0], right[0], left[1], right[1])
        hv = ("x", "y")
    elif args.axis == "x":                   # image [z, y]
        extent = (left[1], right[1], zlo_e, zhi_e)
        hv = ("y", "z")
    else:                                    # image [z, x]
        extent = (left[0], right[0], zlo_e, zhi_e)
        hv = ("x", "z")

    width, height = extent[1] - extent[0], extent[3] - extent[2]
    n = len(fields)
    pw = 3.2 if args.axis == "z" else min(3.2, max(1.6, 8.0 * width / height))
    ph = pw * height / width
    fig, axs = plt.subplots(1, n, figsize=((pw + 0.3) * n + 0.6, ph + 1.6), squeeze=False, layout="constrained")
    axs = axs[0]

    if stars["mass"].size and args.stars:
        if args.zrange:
            inz = (stars["z"] >= zlo) & (stars["z"] <= zhi)
            stars = {k: v[inz] for k, v in stars.items()}
        accreting = stars["tform"] < 0.0                                # sink below mass_SN*n_SN: no SN loadout yet
        young = ~accreting & ((time - stars["tform"]) < args.young)      # SN loadout triggered within --young Myr
        old = ~accreting & ~young

    rho = data["density"]
    for a, (spec, ds, mode, fac, label, cmap, log, clim) in zip(axs, fields):
        img = project(data[ds], rho, mode, ax, lmb) * fac
        clim = args.clim.get(spec.partition(":")[0], clim)
        if log is None:                     # generic dataset: log scale if positive and spanning > 2 decades
            pos = img[img > 0]
            log = bool(img.min() > 0 and pos.size and pos.max() / pos.min() > 100)
        if clim is None:
            vals = img[np.isfinite(img)]
            vals = vals[vals > 0] if log else vals
            clim = (np.percentile(vals, 0.5), np.percentile(vals, 99.5)) if vals.size else (0, 1)
            if not log and cmap == "RdBu_r":       # velocities: symmetric limits
                m = max(abs(clim[0]), abs(clim[1]))
                clim = (-m, m)
            if clim[0] == clim[1]:
                clim = (clim[0] - 0.5, clim[1] + 0.5) if not log else (clim[0] / 2, clim[1] * 2)
        norm = LogNorm(*clim, clip=True) if log else Normalize(*clim)
        im = a.imshow(np.where(img > 0, img, np.nan) if log else img, origin="lower", extent=extent, cmap=cmap,
                      norm=norm, interpolation="nearest", aspect="equal")
        a.set_facecolor(plt.get_cmap(cmap)(0.0))
        fig.colorbar(im, ax=a, location="top", shrink=0.9, pad=0.01).set_label(label, fontsize="small")
        a.set_xlabel(f"{hv[0]} [pc]")
        if a is axs[0]:
            a.set_ylabel(f"{hv[1]} [pc]")
        else:
            a.tick_params(labelleft=False)
        if stars["mass"].size and args.stars:
            h, v = stars[hv[0]], stars[hv[1]]
            a.scatter(h[accreting], v[accreting], s=4, c="w", edgecolors="k", linewidths=0.2, alpha=0.8,
                      label=f"accreting sink ({accreting.sum()})")
            a.scatter(h[young], v[young], s=7, c="cyan", edgecolors="k", linewidths=0.2,
                      label=f"SN loadout < {args.young:g} Myr ago ({young.sum()})")
            a.scatter(h[old], v[old], s=1, c="0.7", edgecolors="none", alpha=0.4, label=f"older ({old.sum()})")

    if stars["mass"].size and args.stars:
        fig.legend(*axs[0].get_legend_handles_labels(), loc="outside lower center", ncol=3, fontsize="small",
                   markerscale=2, frameon=False)
    run_id = os.path.basename(fn).rsplit("_", 1)[0]
    nstar = stars["mass"].size if args.stars else 0
    fig.suptitle(f"{run_id}   t = {time:.1f} Myr   projection along {args.axis}"
                 + (f"   {nstar} star particles, {stars['mass'].sum():.3g} M$_\\odot$" if nstar else ""),
                 fontsize="medium")
    out = os.path.join(outdir, f"proj_{args.axis}_{os.path.splitext(os.path.basename(fn))[0]}.png")
    fig.savefig(out, dpi=args.dpi, bbox_inches="tight")
    plt.close(fig)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", nargs="+", help="snapshot files, or run directories (all their snapshots)")
    ap.add_argument("--axis", choices=("z", "x", "y"), default="x",
                    help="projection axis: z gives the xy plane, x the yz plane, y the xz plane (default: x)")
    ap.add_argument("--fields", nargs="+", default=["dens", "temp", "pres", "vlos"],
                    help="fields to project (default: dens temp pres vlos), see above")
    ap.add_argument("--zrange", type=float, nargs=2, metavar=("ZMIN", "ZMAX"), help="restrict z to this range [pc]")
    ap.add_argument("--clim", nargs="+", default=[], metavar="FIELD=LO,HI", help="colour limits, e.g. dens=1e-2,1e2")
    ap.add_argument("--every", type=int, default=1, help="use every N-th snapshot of a run directory")
    ap.add_argument("--young", type=float, default=10.0,
                    help="highlight particles whose last SN loadout (formation_time) is younger than this [Myr]")
    ap.add_argument("--no-stars", dest="stars", action="store_false", help="do not overplot star particles")
    ap.add_argument("--dpi", type=int, default=150)
    ap.add_argument("-o", "--outdir", default="proj", help="output directory (default: proj)")
    args = ap.parse_args()

    clim = {}
    for c in args.clim:
        try:
            k, v = c.split("=")
            lo, hi = (float(x) for x in v.split(","))
        except ValueError:
            raise SystemExit(f"--clim {c}: expected FIELD=LO,HI")
        clim[k] = (lo, hi)
    args.clim = clim

    files = []
    for p in args.inputs:
        if os.path.isdir(p):
            files += sorted(glob.glob(os.path.join(p, "*_[0-9][0-9][0-9][0-9].h5")))[::args.every]
        else:
            files.append(p)
    if not files:
        raise SystemExit("no snapshots found")

    os.makedirs(args.outdir, exist_ok=True)
    for fn in files:
        print("wrote", plot_snapshot(fn, args, args.outdir))


if __name__ == "__main__":
    main()
