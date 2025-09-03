import h5py
import numpy as np

def _get_attr_any(g, *names, default=None):
    for n in names:
        if n in g.attrs: return g.attrs[n]
    return default

def load_amr_uniform(fname, field_names=None):
    """
    Rasterize a PIERNIK AMR file onto a uniform finest-resolution grid.

    Returns:
      data: dict[field] -> ndarray (Nz, Ny, Nx) at finest dx
      origin: 3-vector left edge
      spacing: 3-vector finest cell size
    """
    data = {}
    with h5py.File(fname, "r") as f:
        sim = f["simulation_parameters"]
        left_edge  = sim.attrs["domain_left_edge"].astype(float)
        right_edge = sim.attrs["domain_right_edge"].astype(float)
        domain_size = right_edge - left_edge

        blocks = [f["data"][k] for k in f["data"].keys()]
        if field_names is None:
            # Take fields from the first block
            field_names = list(blocks[0].keys())

        # ---- collect block metadata ----
        meta = []
        levels = []
        dxs = []
        for g in blocks:
            # level key can be 'level', 'lev', or 'l'
            lev = int(_get_attr_any(g, "level", "lev", "l", default=0))
            levels.append(lev)

            # Try direct geometry; else infer from dims and level
            left = _get_attr_any(g, "left_edge", default=None)
            right = _get_attr_any(g, "right_edge", default=None)
            nb = np.array(_get_attr_any(g, "n_b", default=g[next(iter(g.keys()))].shape[::-1]), dtype=int)  # (Nx,Ny,Nz)
            ng = int(_get_attr_any(g, "ng", "n_ghost", "ghost", default=0))
            # Per-level dx
            dxyz = _get_attr_any(g, "delta", "dl", default=None)
            if dxyz is not None:
                dxyz = np.array(dxyz, dtype=float)
            else:
                # Fallback: infer dx_lev from base dims and level
                # Assume base-level dims are max nb at the minimum level
                # (this matches typical block-structured AMR)
                # Use x-dim only (assumed isotropic)
                pass

            if left is None or right is None:
                # Fall back to offset+dx if present
                off = np.array(_get_attr_any(g, "off", default=[0,0,0]), dtype=int)
                # Need dx for this block:
                if dxyz is None:
                    # Infer dx_lev from level using base resolution
                    # Find minimum level across all blocks on first pass later
                    pass
                # compute physical edges from off & dx
                # IMPORTANT: off is in *active* cells; shift by ng to get interior
                left = left_edge + (off)* (dxyz if dxyz is not None else 1.0)
                right = left + nb * (dxyz if dxyz is not None else 1.0)

            meta.append(dict(group=g, level=lev, nb=nb, ng=ng,
                             left=np.array(left, float), right=np.array(right, float),
                             dxyz=np.array(dxyz if dxyz is not None else 0.0, float)))

            if dxyz is not None:
                dxs.append(dxyz)
        # Finest spacing
        if dxs:
            finest_dx = np.min(np.vstack(dxs), axis=0)
        else:
            # Last resort: infer finest spacing from highest level and domain size
            Lmax = max(levels) if levels else 0
            # estimate base dims from any level-0 block:
            base_nb = None
            for m in meta:
                if m["level"] == 0:
                    base_nb = m["nb"]; break
            if base_nb is None:
                base_nb = meta[0]["nb"]
            finest_dx = domain_size / (base_nb * (2**Lmax))

        # Global dims on finest grid
        global_dims = np.round(domain_size / finest_dx).astype(int)

        # Allocate outputs
        for fld in field_names:
            data[fld] = np.full((global_dims[2], global_dims[1], global_dims[0]), np.nan, dtype=np.float32)

        # Composite from coarse to fine (fine overwrites)
        order = np.argsort(levels)
        for idx in order:
            m = meta[idx]
            g = m["group"]
            lev = m["level"]
            dx = m["dxyz"] if np.any(m["dxyz"]) else finest_dx * (2**(max(levels)-lev))
            # integer upscale factor from this level to finest
            scale = np.round(dx / finest_dx).astype(int)

            # Compute finest-grid index window for this block (interior only)
            # Skip ghost if present
            nb_int = m["nb"]  # assume stored as interior count
            left_phys = m["left"]
            right_phys = m["right"]

            start = np.floor((left_phys - left_edge) / finest_dx).astype(int)
            stop  = np.floor((right_phys - left_edge) / finest_dx).astype(int)

            # slices on finest grid
            z0,z1 = start[2], stop[2]
            y0,y1 = start[1], stop[1]
            x0,x1 = start[0], stop[0]

            for fld in field_names:
                arr = g[fld][:]
                # Expect on-disk order (z,y,x). Remove ghosts if they’re present in the dataset:
                gz = max(0, arr.shape[0] - nb_int[2]) // 2 if arr.shape[0] != nb_int[2] else 0
                gy = max(0, arr.shape[1] - nb_int[1]) // 2 if arr.shape[1] != nb_int[1] else 0
                gx = max(0, arr.shape[2] - nb_int[0]) // 2 if arr.shape[2] != nb_int[0] else 0
                arr = arr[gz:arr.shape[0]-gz, gy:arr.shape[1]-gy, gx:arr.shape[2]-gx]

                # Upsample to finest grid by integer repetition
                az = np.repeat(arr, scale[2], axis=0)
                ayz = np.repeat(az, scale[1], axis=1)
                axyz = np.repeat(ayz, scale[0], axis=2)

                # Trim in case boundaries are not exact multiples
                tz = z1 - z0; ty = y1 - y0; tx = x1 - x0
                axyz = axyz[:tz, :ty, :tx]

                data[fld][z0:z1, y0:y1, x0:x1] = axyz  # fine overwrites coarse

        return data, left_edge, finest_dx, global_dims

#%%
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# --- import the AMR loader you defined earlier ---
# from your_module import load_amr_uniform

plt.rcParams.update({
    'figure.dpi': 720, 'savefig.dpi': 720,
    'axes.linewidth': 2.5, 'font.size': 12,
    'xtick.major.size': 5, 'ytick.major.size': 5,
    'xtick.major.width': 1.0, 'ytick.major.width': 1.0,
    'xtick.direction': 'out', 'ytick.direction': 'out',
    'legend.frameon': True, 'legend.edgecolor': 'black', 'legend.framealpha': 1.0,
    'figure.figsize': (8, 6),
})

# ---- helpers ---------------------------------------------------------------
def pick_key(keys_available, *candidates):
    """Return the first candidate present in keys_available (raise if none)."""
    for c in candidates:
        if c in keys_available:
            return c
    raise KeyError(f"None of {candidates} found. Available: {sorted(keys_available)}")

def build_coords(origin, dx, dims):
    """Return x,y center coords for a (Nz,Ny,Nx) array rasterized at finest dx."""
    Nx, Ny, Nz = dims[0], dims[1], dims[2]
    x = origin[0] + (np.arange(Nx) + 0.5) * dx[0]
    y = origin[1] + (np.arange(Ny) + 0.5) * dx[1]
    return x, y

# ---- files (use different files for diffusion/streaming if you have them) --
file_diff = '/home/vinodvp/simdir/piernik/runs/test10/scr_tst_0000.h5'
file_stream = '/home/vinodvp/simdir/piernik/runs/test10/scr_tst_0005.h5'

# ---- load (AMR-aware) ------------------------------------------------------
dataD, originD, dxD, dimsD = load_amr_uniform(file_diff)
dataS, originS, dxS, dimsS = load_amr_uniform(file_stream)

# coords (finest grid)
xD, yD = build_coords(originD, dxD, dimsD)
xS, yS = build_coords(originS, dxS, dimsS)

# mid-plane index in z
zD = 0 if dataD[next(iter(dataD))].shape[0] == 1 else dataD[next(iter(dataD))].shape[0] // 2
zS = 0 if dataS[next(iter(dataS))].shape[0] == 1 else dataS[next(iter(dataS))].shape[0] // 2

# resolve field names
keysD = set(dataD.keys())
keysS = set(dataS.keys())

k_rhoD = pick_key(keysD, 'density', 'dens')
k_vxD  = pick_key(keysD, 'velocity_x', 'velx', 'vx')
k_vyD  = pick_key(keysD, 'velocity_y', 'vely', 'vy')
k_bxD  = pick_key(keysD, 'mag_field_x', 'magx', 'bx', 'b_x') if keysD & {'mag_field_x','magx','bx','b_x'} else None
k_byD  = pick_key(keysD, 'mag_field_y', 'magy', 'by', 'b_y') if keysD & {'mag_field_y','magy','by','b_y'} else None
k_ecD  = pick_key(keysD, 'escr_01', 'escr')

k_rhoS = pick_key(keysS, 'density', 'dens')
k_vxS  = pick_key(keysS, 'velocity_x', 'velx', 'vx')
k_vyS  = pick_key(keysS, 'velocity_y', 'vely', 'vy')
k_bxS  = pick_key(keysS, 'mag_field_x', 'magx', 'bx', 'b_x') if keysS & {'mag_field_x','magx','bx','b_x'} else None
k_byS  = pick_key(keysS, 'mag_field_y', 'magy', 'by', 'b_y') if keysS & {'mag_field_y','magy','by','b_y'} else None
k_ecS  = pick_key(keysS, 'escr_01', 'escr')

# ---- plot ------------------------------------------------------------------
fig, ax = plt.subplots(2, 2, figsize=(8, 6))

# (0,0) density + velocity (diffusion)
im = ax[0,0].imshow(dataD[k_rhoD][zD,:,:], extent=[xD[0], xD[-1], yD[0], yD[-1]],
                    origin='lower', cmap='RdGy_r', aspect='equal')
# ax[0,0].streamplot(xD, yD, dataD[k_vxD][zD,:,:], dataD[k_vyD][zD,:,:], density=0.8, color='k')
ax[0,0].set(title='t=0.1, diffusion', xlabel=r'$\mathbf{x}$', ylabel=r'$\mathbf{y}$')
cb = fig.colorbar(im, ax=ax[0,0], fraction=0.046, pad=0.04); cb.set_label(r'$\mathbf{\rho}$', fontweight='bold')

# (0,1) Ec + (optional) B streamlines (diffusion)
im = ax[0,1].imshow(dataD[k_ecD][zD,:,:], extent=[xD[0], xD[-1], yD[0], yD[-1]],
                    origin='lower', cmap='RdGy_r', aspect='equal')
if k_bxD and k_byD:
    ax[0,1].streamplot(xD, yD, dataD[k_bxD][zD,:,:], dataD[k_byD][zD,:,:], density=0.5, color='white')
ax[0,1].set(title='t=0.1, diffusion', xlabel=r'$\mathbf{x}$', ylabel=r'$\mathbf{y}$')
cb = fig.colorbar(im, ax=ax[0,1], fraction=0.046, pad=0.04); cb.set_label(r'$\mathbf{E_c}$', fontweight='bold')

# (1,0) density + velocity (streaming)
im = ax[1,0].imshow(dataS[k_rhoS][zS,:,:], extent=[xS[0], xS[-1], yS[0], yS[-1]],
                    origin='lower', cmap='RdGy_r', aspect='equal')
# ax[1,0].streamplot(xS, yS, dataS[k_vxS][zS,:,:], dataS[k_vyS][zS,:,:], density=0.8, color='k')
ax[1,0].set(title='t=0.1, streaming', xlabel=r'$\mathbf{x}$', ylabel=r'$\mathbf{y}$')
cb = fig.colorbar(im, ax=ax[1,0], fraction=0.046, pad=0.04); cb.set_label(r'$\mathbf{\rho}$', fontweight='bold')
cb.set_ticks(np.arange(0.15,1.20,0.15))
# (1,1) Ec + (optional) B streamlines (streaming)
im = ax[1,1].imshow(dataS[k_ecS][zS,:,:], extent=[xS[0], xS[-1], yS[0], yS[-1]],
                    origin='lower', cmap='RdGy_r', aspect='equal')
# if k_bxS and k_byS:
    # ax[1,1].streamplot(xS, yS, dataS[k_bxS][zS,:,:], dataS[k_byS][zS,:,:], density=0.5, color='white')
ax[1,1].set(title='t=0.1, streaming', xlabel=r'$\mathbf{x}$', ylabel=r'$\mathbf{y}$')
cb = fig.colorbar(im, ax=ax[1,1], fraction=0.046, pad=0.04); cb.set_label(r'$\mathbf{E_c}$', fontweight='bold')

plt.tight_layout()
for a in ax.ravel():
    for s in a.spines.values(): s.set_linewidth(3)

fig.add_artist(mpatches.Rectangle((0.005, 0.005), 0.99, 0.99,
                                  transform=fig.transFigure, fill=False, lw=4, ec='black'))
plt.savefig('escr_01.png', dpi=720)
