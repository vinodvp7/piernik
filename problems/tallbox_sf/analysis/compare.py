#!/usr/bin/env python3
"""
Compare tallbox_sf runs (HD / MHD / diffusive CR / two-moment CR / streaming CR).

Time series come from the .tsl files, snapshot diagnostics from the HDF5 outputs (uniform grid, PSM units:
pc, Msun, Myr). Usage:

    python compare.py RUNDIR [RUNDIR ...] [-o OUTDIR] [-t 100 250 500] [--labels A B ...]

e.g. from piernik/runs:
    python ../problems/tallbox_sf/analysis/compare.py tallbox_sf_A_hd tallbox_sf_B_mhd tallbox_sf_C1_crdiff \
           tallbox_sf_C2_scrdiff tallbox_sf_D_scr -t 100 250 500 -o compare_20pc
"""

import argparse
import glob
import os

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# PSM code units -> physical
NH_PER_RHO = 6.770e-23 / (1.4 * 1.6726e-24)      # 1 Msun/pc^3 -> n_H [cm^-3] (1.4 m_H per H)
POK_PER_P = 6.47e-13 / 1.3807e-16                 # 1 Msun/(pc Myr^2) -> P/k_B [K cm^-3]
GAMMA_CR = 4.0 / 3.0
COLORS = {}                                        # label -> colour, fixed for all panels (filled in main)


# ----------------------------------------------------------------------------- readers

def read_tsl(rundir):
    """Concatenate all *.tsl files of a run (restarts write new ones) into a dict of arrays, sorted by time."""
    files = sorted(glob.glob(os.path.join(rundir, "*.tsl")))
    if not files:
        return None
    cols, rows = None, []
    for fn in files:
        with open(fn) as f:
            for line in f:
                if line.startswith("#"):
                    if cols is None and len(line.split()) > 2:
                        cols = line.lstrip("#").split()
                    continue
                if line.strip():
                    rows.append([float(x) for x in line.split()])
    data = np.array(rows)
    order = np.argsort(data[:, cols.index("time")], kind="stable")
    data = data[order]
    _, keep = np.unique(data[:, cols.index("time")][::-1], return_index=True)   # last occurrence wins after a restart
    data = data[len(data) - 1 - keep]
    return {c: data[:, i] for i, c in enumerate(cols)}


def snapshot_times(rundir):
    out = []
    for fn in sorted(glob.glob(os.path.join(rundir, "*_[0-9][0-9][0-9][0-9].h5"))):
        with h5py.File(fn, "r") as f:
            out.append((float(f["simulation_parameters"].attrs["current_time"][0]), fn))
    return out


def read_snapshot(fn, fields):
    """Assemble uniform-grid fields (z, y, x order) plus geometry and star particles."""
    with h5py.File(fn, "r") as f:
        sp = f["simulation_parameters"].attrs
        nd = sp["domain_dimensions"][::-1]
        le, re = sp["domain_left_edge"], sp["domain_right_edge"]
        out = {"t": float(sp["current_time"][0])}
        avail = set()
        for name, g in f["data"].items():
            avail |= set(k for k in g.keys() if isinstance(g[k], h5py.Dataset))
        for fld in fields:
            if fld in avail:
                out[fld] = np.zeros(nd)
        stars = {"mass": [], "position_x": [], "position_y": [], "position_z": [], "formation_time": []}
        for name, g in f["data"].items():
            o = g.attrs["off"][::-1]
            n = g.attrs["n_b"].astype(int)[::-1]
            sl = tuple(slice(o[d], o[d] + n[d]) for d in range(3))
            for fld in fields:
                if fld in out:
                    out[fld][sl] = g[fld][...]
            if "particles" in g and "stars" in g["particles"]:
                for k in stars:
                    if k in g["particles"]["stars"]:
                        stars[k].append(g["particles"]["stars"][k][...])
        out["stars"] = {k: (np.concatenate(v) if v else np.zeros(0)) for k, v in stars.items()}
    dx = (re - le) / nd[::-1]
    out["dx"] = dx
    out["x"] = le[0] + (np.arange(nd[2]) + 0.5) * dx[0]
    out["y"] = le[1] + (np.arange(nd[1]) + 0.5) * dx[1]
    out["z"] = le[2] + (np.arange(nd[0]) + 0.5) * dx[2]
    return out


def cr_pressure(s):
    """CR pressure from whichever CR field the run has (diffusive 'cr_p+' or two-moment 'escr_01')."""
    for k in ("cr_p+", "escr_01"):
        if k in s:
            return (GAMMA_CR - 1.0) * s[k]
    return None


# ----------------------------------------------------------------------------- plots

def plot_timeseries(runs, outdir):
    panels = [("SFR10_kpc2", r"$\Sigma_{\rm SFR}$ (10 Myr) [M$_\odot$ yr$^{-1}$ kpc$^{-2}$]", True),
              ("Mstar", r"$M_*$ [M$_\odot$]", True),
              ("H_gas", r"$\langle|z|\rangle_m$ [pc]", False),
              ("sigma_v", r"$\sigma_v$ (mass-weighted) [pc/Myr]", False),
              ("Mdot_z2", r"$\dot M$ through $|z|$ = z_flux(2) [M$_\odot$/Myr]", False),
              ("eta", r"mass loading $\dot M(z_2)/\dot M_*$", True),
              ("E_cr/E_th", r"$E_{\rm CR}/E_{\rm th}$", True),
              ("emag", r"$E_{\rm mag}$ [code]", True)]
    fig, axs = plt.subplots(4, 2, figsize=(12, 14), sharex=True)
    for ax, (key, lab, logy) in zip(axs.flat, panels):
        for label, ts in runs.items():
            if ts is None:
                continue
            t = ts["time"]
            if key == "eta":
                if "Mdot_z2" not in ts or "SFR40_kpc2" not in ts:
                    continue
                sfr_myr = ts["SFR40_kpc2"] * 1e6 * ts.get("area_kpc2", 1.0)     # SFR columns are per kpc^2
                y = np.where(sfr_myr > 0, ts["Mdot_z2"] / np.maximum(sfr_myr, 1e-30), np.nan)
            elif key == "E_cr/E_th":
                if "E_cr" not in ts or "eint" not in ts or not np.any(ts["E_cr"] > 0):
                    continue
                y = ts["E_cr"] / ts["eint"]
            else:
                if key not in ts:
                    continue
                y = ts[key]
            if logy:
                y = np.where(y > 0, y, np.nan)
            ax.plot(t, y, label=label, color=COLORS.get(label), lw=1.2)
        ax.set_ylabel(lab)
        if logy:
            ax.set_yscale("log")
        ax.grid(alpha=0.3)
    for ax in axs[-1]:
        ax.set_xlabel("t [Myr]")
    axs[0, 0].legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "timeseries.png"), dpi=130)
    plt.close(fig)


def plot_profiles(snaps, outdir, tag):
    fig, axs = plt.subplots(2, 3, figsize=(15, 8), sharex=True)
    axs = axs.flat
    for label, s in snaps.items():
        if s is None:
            continue
        z = s["z"]
        rho = s["density"]
        m = rho.sum(axis=(1, 2))
        axs[0].semilogy(z, rho.mean(axis=(1, 2)) * NH_PER_RHO, label=label, color=COLORS.get(label))
        if "temperature" in s:
            axs[1].semilogy(z, (s["temperature"] * rho).sum(axis=(1, 2)) / m, label=label, color=COLORS.get(label))
        if "pressure" in s:
            axs[2].semilogy(z, s["pressure"].mean(axis=(1, 2)) * POK_PER_P, label=label, color=COLORS.get(label))
        if "mag_field_x" in s:
            pmag = 0.5 * (s["mag_field_x"] ** 2 + s["mag_field_y"] ** 2 + s["mag_field_z"] ** 2)
            axs[3].semilogy(z, pmag.mean(axis=(1, 2)) * POK_PER_P, label=label, color=COLORS.get(label))
        pcr = cr_pressure(s)
        if pcr is not None:
            axs[4].semilogy(z, pcr.mean(axis=(1, 2)) * POK_PER_P, label=label, color=COLORS.get(label))
        if "velocity_z" in s:
            axs[5].plot(z, (rho * s["velocity_z"]).sum(axis=(1, 2)) / m, label=label, color=COLORS.get(label))
    for ax, lab in zip(axs, [r"$\langle n_H \rangle$ [cm$^{-3}$]", r"$\langle T \rangle_m$ [K]", r"$\langle P_{\rm th}\rangle/k_B$ [K cm$^{-3}$]",
                             r"$\langle P_{\rm mag}\rangle/k_B$", r"$\langle P_{\rm CR}\rangle/k_B$", r"$\langle v_z \rangle_m$ [pc/Myr]"]):
        ax.set_ylabel(lab)
        ax.grid(alpha=0.3)
    for ax in list(axs)[3:]:
        ax.set_xlabel("z [pc]")
    axs[0].legend(fontsize=9)
    fig.suptitle(f"horizontally averaged profiles, {tag}")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, f"profiles_{tag}.png"), dpi=130)
    plt.close(fig)


def plot_slices(snaps, outdir, tag):
    labels = [l for l, s in snaps.items() if s is not None]
    if not labels:
        return
    fig, axs = plt.subplots(3, len(labels), figsize=(3.2 * len(labels), 13), squeeze=False)
    for j, label in enumerate(labels):
        s = snaps[label]
        ny = s["density"].shape[1] // 2
        nz = s["density"].shape[0] // 2
        ext_xz = [s["x"][0], s["x"][-1], s["z"][0], s["z"][-1]]
        ext_xy = [s["x"][0], s["x"][-1], s["y"][0], s["y"][-1]]
        im0 = axs[0, j].imshow(np.log10(s["density"][:, ny, :] * NH_PER_RHO), origin="lower", extent=ext_xz, aspect="auto", cmap="viridis", vmin=-4, vmax=2)
        axs[0, j].set_title(label)
        if "temperature" in s:
            im1 = axs[1, j].imshow(np.log10(s["temperature"][:, ny, :]), origin="lower", extent=ext_xz, aspect="auto", cmap="inferno", vmin=1.5, vmax=7.5)
        im2 = axs[2, j].imshow(np.log10(s["density"][nz, :, :] * NH_PER_RHO), origin="lower", extent=ext_xy, cmap="viridis", vmin=-3, vmax=2)
        st = s["stars"]
        if st["mass"].size:
            young = (s["t"] - st["formation_time"]) < 40.0
            axs[2, j].scatter(st["position_x"][young], st["position_y"][young], s=3, c="w", lw=0)
        axs[2, j].set_xlabel("x [pc]")
    axs[0, 0].set_ylabel("z [pc] (edge-on, y=0)")
    axs[1, 0].set_ylabel("z [pc] (edge-on, y=0)")
    axs[2, 0].set_ylabel("y [pc] (face-on, z=0)")
    fig.colorbar(im0, ax=axs[0, :].tolist(), label=r"log $n_H$ [cm$^{-3}$]")
    if "temperature" in snaps[labels[0]]:
        fig.colorbar(im1, ax=axs[1, :].tolist(), label="log T [K]")
    fig.colorbar(im2, ax=axs[2, :].tolist(), label=r"log $n_H$ [cm$^{-3}$] (white: stars < 40 Myr)")
    fig.suptitle(f"slices, {tag}")
    fig.savefig(os.path.join(outdir, f"slices_{tag}.png"), dpi=110)
    plt.close(fig)


def plot_phase(snaps, outdir, tag):
    labels = [l for l, s in snaps.items() if s is not None and "temperature" in s]
    if not labels:
        return
    fig, axs = plt.subplots(1, len(labels), figsize=(3.6 * len(labels), 3.6), squeeze=False, sharey=True)
    for j, label in enumerate(labels):
        s = snaps[label]
        n = np.log10(np.maximum(s["density"] * NH_PER_RHO, 1e-8)).ravel()
        T = np.log10(np.maximum(s["temperature"], 1.0)).ravel()
        h, xe, ye = np.histogram2d(n, T, bins=[np.linspace(-5, 3, 120), np.linspace(1, 8, 120)], weights=s["density"].ravel())
        axs[0, j].pcolormesh(xe, ye, np.log10(np.maximum(h.T / h.sum(), 1e-8)), cmap="magma", vmin=-6, vmax=-1)
        axs[0, j].set_title(label)
        axs[0, j].set_xlabel(r"log $n_H$ [cm$^{-3}$]")
    axs[0, 0].set_ylabel("log T [K]")
    fig.suptitle(f"mass-weighted phase diagram, {tag}")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, f"phase_{tag}.png"), dpi=130)
    plt.close(fig)


# ----------------------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("rundirs", nargs="+")
    ap.add_argument("--labels", nargs="+", help="legend labels (default: run directory names)")
    ap.add_argument("-t", "--times", nargs="+", type=float, default=[100.0, 250.0, 500.0], help="snapshot times [Myr] (nearest output is used)")
    ap.add_argument("-o", "--outdir", default="compare_out")
    args = ap.parse_args()

    labels = args.labels or [os.path.basename(os.path.normpath(r)) for r in args.rundirs]
    COLORS.update({l: f"C{i}" for i, l in enumerate(labels)})
    os.makedirs(args.outdir, exist_ok=True)

    runs = {l: read_tsl(r) for l, r in zip(labels, args.rundirs)}
    avail = {l: snapshot_times(r) for l, r in zip(labels, args.rundirs)}
    for l in labels:
        if runs[l] is not None and avail[l]:
            with h5py.File(avail[l][0][1], "r") as f:
                sp = f["simulation_parameters"].attrs
                L = (sp["domain_right_edge"] - sp["domain_left_edge"]) / 1000.0
                runs[l]["area_kpc2"] = L[0] * L[1]
    plot_timeseries(runs, args.outdir)

    fields = ["density", "temperature", "pressure", "velocity_z", "mag_field_x", "mag_field_y", "mag_field_z", "cr_p+", "escr_01"]
    for t in args.times:
        snaps = {}
        for l in labels:
            if not avail[l]:
                snaps[l] = None
                continue
            tt, fn = min(avail[l], key=lambda p: abs(p[0] - t))
            snaps[l] = read_snapshot(fn, fields)
            print(f"{l}: t_req = {t:g} Myr -> {os.path.basename(fn)} (t = {tt:.2f} Myr)")
        tag = f"t{t:g}Myr"
        plot_profiles(snaps, args.outdir, tag)
        plot_slices(snaps, args.outdir, tag)
        plot_phase(snaps, args.outdir, tag)

    print(f"plots written to {args.outdir}/")


if __name__ == "__main__":
    main()
