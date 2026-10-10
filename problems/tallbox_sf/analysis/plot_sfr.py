#!/usr/bin/env python3
"""
Star formation rate of tallbox_sf runs, averaged over time bins of width dt.

The stellar mass formed comes from the dmass_stars_tot column of the .tsl files (cumulative mass converted into star
particles, Msun). It is not saved in restart files, so it restarts from 0 after every restart: increments are taken
within each .tsl file, and where a restarted run overlaps an earlier file the later file wins. Each increment is put in
the bin containing the end of its tsl interval, so dt should be well above dt_tsl.

Usage:

    python plot_sfr.py RUNDIR [RUNDIR ...] --dt 1 10 [--per-area] [--tmin T] [--tmax T] [--labels A B ...] [-o sfr.png]

e.g. from piernik/runs:
    python ../problems/tallbox_sf/analysis/plot_sfr.py tallbox_sf_A_20pc tallbox_sf_D_20pc --dt 1 10 --per-area
"""

import argparse
import glob
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

MASS_COL = "dmass_stars_tot"
LINESTYLES = ["-", "--", ":", "-."]


def read_tsl_files(rundir):
    """List of (time, cumulative stellar mass) arrays, one per .tsl file, in file order."""
    out = []
    for fn in sorted(glob.glob(os.path.join(rundir, "*.tsl"))):
        cols, rows = None, []
        with open(fn) as f:
            for line in f:
                if line.startswith("#"):
                    if cols is None and len(line.split()) > 2:
                        cols = line.lstrip("#").split()
                    continue
                if line.strip():
                    rows.append([float(x) for x in line.split()])
        if not rows:
            continue
        if MASS_COL not in cols:
            raise SystemExit(f"{fn}: no '{MASS_COL}' column (built without NBODY/star formation?)")
        data = np.array(rows)
        out.append((data[:, cols.index("time")], data[:, cols.index(MASS_COL)]))
    if not out:
        raise SystemExit(f"{rundir}: no .tsl files")
    return out


def formed_mass_increments(segments):
    """Times and stellar mass formed in each tsl interval, merged over restarts (a later file overrides overlapping times)."""
    t_all, dm_all = [], []
    for n, (t, m) in enumerate(segments):
        dm = np.diff(m, prepend=m[0])
        keep = np.ones_like(t, dtype=bool)
        if n > 0:
            keep[0] = False                    # first row of a restarted file: its increment is taken from the previous file
        if n + 1 < len(segments):
            keep &= t <= segments[n + 1][0][0]  # drop what the next (restarted) file recomputes
        t_all.append(t[keep])
        dm_all.append(dm[keep])
    t = np.concatenate(t_all)
    dm = np.concatenate(dm_all)
    order = np.argsort(t, kind="stable")
    return t[order], dm[order]


def domain_area_kpc2(rundir):
    """Lx * Ly in kpc^2 from the run's problem.par (PSM units: pc)."""
    fn = os.path.join(rundir, "problem.par")
    with open(fn) as f:
        txt = f.read()
    v = {}
    for key in ("xmin", "xmax", "ymin", "ymax"):
        m = re.search(rf"^\s*{key}\s*=\s*([-+0-9.eEdD]+)", txt, re.M)
        if m is None:
            raise SystemExit(f"{fn}: cannot find {key}")
        v[key] = float(m.group(1).replace("d", "e").replace("D", "e"))
    return (v["xmax"] - v["xmin"]) * (v["ymax"] - v["ymin"]) / 1.0e6


def binned_sfr(t, dm, dt, tmin, tmax, keep_partial):
    """Bin edges and SFR [Msun/yr] in bins of width dt [Myr] starting at tmin."""
    nbin = int(np.floor((tmax - tmin) / dt + 1e-9))
    if keep_partial and tmin + nbin * dt < tmax - 1e-9:
        nbin += 1
    if nbin < 1:
        return None, None
    edges = tmin + dt * np.arange(nbin + 1)
    sel = (t > tmin + 1e-9 * dt) & (t <= edges[-1] + 1e-9 * dt)
    idx = np.minimum(np.ceil((t[sel] - tmin) / dt - 1e-9).astype(int) - 1, nbin - 1)   # tolerance: a time on a bin edge closes that bin
    mass = np.bincount(np.maximum(idx, 0), weights=dm[sel], minlength=nbin)
    width = np.diff(edges)
    width[-1] = min(width[-1], tmax - edges[-2])    # partial last bin: average over the time actually covered
    return edges, mass / (width * 1.0e6)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("rundirs", nargs="+", help="run directories containing the .tsl files")
    ap.add_argument("--dt", type=float, nargs="+", default=[10.0], help="bin width(s) in Myr (default: 10)")
    ap.add_argument("--per-area", action="store_true", help="plot Sigma_SFR [Msun/yr/kpc^2] instead of SFR [Msun/yr]")
    ap.add_argument("--tmin", type=float, default=0.0, help="start of the first bin [Myr] (default: 0)")
    ap.add_argument("--tmax", type=float, default=None, help="end of the binned range [Myr] (default: last tsl time)")
    ap.add_argument("--keep-partial", action="store_true", help="also plot the incomplete last bin")
    ap.add_argument("--log", action="store_true", help="logarithmic y axis")
    ap.add_argument("--labels", nargs="+", help="legend labels, one per run directory")
    ap.add_argument("-o", "--output", default="sfr.png", help="output image (default: sfr.png)")
    args = ap.parse_args()

    if args.labels and len(args.labels) != len(args.rundirs):
        raise SystemExit("--labels needs one label per run directory")
    if any(dt <= 0 for dt in args.dt):
        raise SystemExit("--dt must be positive")

    fig, ax = plt.subplots(figsize=(8, 4.5))
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    for r, rundir in enumerate(args.rundirs):
        label = args.labels[r] if args.labels else os.path.basename(os.path.normpath(rundir))
        t, dm = formed_mass_increments(read_tsl_files(rundir))
        tmax = t[-1] if args.tmax is None else min(args.tmax, t[-1])
        norm = domain_area_kpc2(rundir) if args.per_area else 1.0
        print(f"{label}: t = {t[0]:.3g} .. {t[-1]:.4g} Myr, stellar mass formed {dm[(t > args.tmin) & (t <= tmax)].sum():.4g} Msun")
        for d, dt in enumerate(sorted(args.dt)):
            edges, sfr = binned_sfr(t, dm, dt, args.tmin, tmax, args.keep_partial)
            if edges is None:
                print(f"  dt = {dt:g} Myr: no complete bin in [{args.tmin:g}, {tmax:.4g}] Myr, skipped")
                continue
            ax.stairs(sfr / norm, edges, color=colors[r % len(colors)], linestyle=LINESTYLES[d % len(LINESTYLES)],
                      label=f"{label}, $\\Delta t$ = {dt:g} Myr")

    ax.set_xlabel("t [Myr]")
    ax.set_ylabel(r"$\Sigma_\mathrm{SFR}$ [M$_\odot$ yr$^{-1}$ kpc$^{-2}$]" if args.per_area else r"SFR [M$_\odot$ yr$^{-1}$]")
    if args.log:
        ax.set_yscale("log")
    ax.legend(fontsize="small")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(args.output, dpi=150)
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
