import sys
import h5py
import numpy as np

# --- VTK Imports ---
# We wrap this in a try/except block to provide a helpful error message
# if the user doesn't have VTK installed.


def get_field_names(h5_file):
    """
    Scans the first data block in the HDF5 file to find all available
    dataset names (e.g., 'density', 'pressure').

    Args:
        h5_file (h5py.File): An open H5py file object.

    Returns:
        list: A list of strings with the names of the data fields.
    """
    print("Scanning for available data fields...")
    field_names = []
    # Get the name of the first grid block (e.g., 'grid_0000000000')
    first_block_name = list(h5_file["data"].keys())[0]
    for name in h5_file["data"][first_block_name]:
        field_names.append(name)
    print(f"Found fields: {field_names}")
    return field_names


def load_and_stitch_data(fname):
    """
    Loads and stitches a 3D field from multiple blocks in an HDF5 file.
    This function is based on the user-provided reference code and the
    HDF5 inspection report.

    Args:
        fname (str): Path to the HDF5 file.
        field_names (list): A list of all field names to extract.

    Returns:
        tuple: A tuple containing:
            - dict: A dictionary of the stitched 3D data arrays for each field.
            - np.ndarray: The global dimensions of the stitched cell grid (Nx, Ny, Nz).
            - np.ndarray: The physical coordinates of the grid origin (x0, y0, z0).
            - np.ndarray: The cell spacing in each direction (dx, dy, dz).
    """
    stitched_data = {}
    with h5py.File(fname, "r") as f:
        field_names = get_field_names(f)
    with h5py.File(fname, "r") as f:
        # --- 1. Read metadata from the HDF5 file ---
        block_names = list(f["data"].keys())
        all_offsets = np.array([f["data"][name].attrs["off"] for name in block_names], dtype=int)
        all_dims = np.array([f["data"][name].attrs["n_b"] for name in block_names], dtype=int)

        # Determine the full size of the stitched grid of cells
        global_cell_dims = np.max(all_offsets + all_dims, axis=0)
        print(f"Global cell grid dimensions determined to be: {global_cell_dims}")

        # Get global grid physical properties from the root attributes
        origin = f["simulation_parameters"].attrs["domain_left_edge"]
        domain_size = f["simulation_parameters"].attrs["domain_right_edge"] - origin

        # For cell-centered data, spacing is domain size / number of cells.
        spacing = domain_size / np.maximum(1, global_cell_dims)

        print(f"Global grid origin: {origin}")
        print(f"Cell spacing: {spacing}")

        # --- 2. Load and stitch each data field ---
        for field in field_names:
            print(f"  Processing field: {field}")
            # Create an empty array to hold the full stitched data
            # The data on disk is (z, y, x), so we create the numpy array in that order.
            stitched_field = np.empty((global_cell_dims[2], global_cell_dims[1], global_cell_dims[0]), dtype=np.float32)

            # Iterate over each block and place its data into the large array
            for i, name in enumerate(block_names):
                block_data = f["data"][name][field][:]  # This is already (z, y, x)
                offset = all_offsets[i]
                dims = all_dims[i]

                # Define the slice where this block goes in the big array
                slc = np.s_[
                    offset[2]: offset[2] + dims[2],
                    offset[1]: offset[1] + dims[1],
                    offset[0]: offset[0] + dims[0],
                ]
                stitched_field[slc] = block_data

            stitched_data[field] = stitched_field

    return stitched_data, global_cell_dims.astype(int), origin, spacing
#%%
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize  # optional
import matplotlib.patches as mpatches
from scipy import special as sp

plt.rcParams.update({
    'figure.dpi': 150, 'savefig.dpi': 1080,
    'axes.linewidth': 2.5,
    'font.size': 12,
    'xtick.major.size': 5, 'ytick.major.size': 5,
    'xtick.major.width': 1.0, 'ytick.major.width': 1.0,
    'xtick.direction': 'out', 'ytick.direction': 'out',
    'legend.frameon': True, 'legend.edgecolor': 'black', 'legend.framealpha': 1.0,
    'figure.figsize': (8, 6)
})


data_to_plot = 'escr_01'

file = '/home/vinodvp/simdir/piernik/runs/test11/scr_dif_0005.h5'
data, cell_dims, origin, spacing = load_and_stitch_data(file)
N, dx, x0 = cell_dims[0], spacing[0], origin[0]
xe = x0 + np.arange(N+1)*dx                   # assume x0 is left edge
x  = 0.5*(xe[:-1] + xe[1:])                   # centers

N, dy, y0 = cell_dims[1], spacing[1], origin[1]
ye = y0 + np.arange(N+1)*dy                   # assume x0 is left edge
y  = 0.5*(ye[:-1] + ye[1:])                   # centers

Xc, Yc = np.meshgrid(x, y, indexing="xy")

fig,ax=plt.subplots(2,2,figsize=(8,6))

d0=data["density"][0,:,:]
ax[0,0].imshow(d0, extent=[x[0], x[-1], y[0], y[-1]],
               origin="lower",
               cmap="viridis")

ax[0,0].streamplot(x,y,
    data["velocity_x"][0,:,:],data["velocity_y"][0,:,:],
    density=0.8,color='k')

ax[0,0].set_xlabel(r'$\mathbf{x}$', fontweight='bold')
ax[0,0].set_ylabel(r'$\mathbf{y}$', fontweight='bold')
ax[0,0].set_title('t=0.1,diffusion', fontweight='bold')
cb = fig.colorbar(ax[0,0].images[-1], ax=ax[0,0], fraction=0.046, pad=0.04)
cb.set_label(r'$\mathbf{\rho}$', fontweight='bold', labelpad=6)  # bold E with subscript c
cb.set_ticks(np.arange(0.15,1.20,0.15))


d0=data["escr_01"][0,:,:]
ax[0,1].imshow(d0, extent=[x[0], x[-1], y[0], y[-1]],
               origin="lower",
               cmap="viridis", norm=LogNorm(0.1,10))

ax[0,1].streamplot(x,y,
    data["mag_field_x"][0,:,:],data["mag_field_y"][0,:,:],
    density=0.5,color='white')

ax[0,1].set_xlabel(r'$\mathbf{x}$', fontweight='bold')
ax[0,1].set_ylabel(r'$\mathbf{y}$', fontweight='bold')
ax[0,1].set_title('t=0.1,diffusion', fontweight='bold')
cb = fig.colorbar(ax[0,1].images[-1] ,ax=ax[0,1], fraction=0.046, pad=0.04)
cb.set_label(r'$\mathbf{E_c}$', fontweight='bold', labelpad=6)  # bold E with subscript c
cb.set_ticks((0.1,1,10))
ax[0,1].set_ylim(np.min(y),np.max(y))

file = '/home/vinodvp/simdir/piernik/runs/test11/scr_stm_0005.h5'
data, cell_dims, origin, spacing = load_and_stitch_data(file)


d0=data["density"][0,:,:]
ax[1,0].imshow(d0, extent=[x[0], x[-1], y[0], y[-1]],
               origin="lower",
               cmap="viridis",norm=Normalize(0.15,1.35))

ax[1,0].streamplot(x,y,
    data["velocity_x"][0,:,:],data["velocity_y"][0,:,:],
    density=0.8,color='k')

ax[1,0].set_xlabel(r'$\mathbf{x}$', fontweight='bold')
ax[1,0].set_ylabel(r'$\mathbf{y}$', fontweight='bold')
ax[1,0].set_title('t=0.1,streaming', fontweight='bold')
cb = fig.colorbar(ax[1,0].images[-1], ax=ax[1,0], fraction=0.046, pad=0.04)
cb.set_label(r'$\mathbf{\rho}$', fontweight='bold', labelpad=6)  # bold E with subscript c
cb.set_ticks(np.arange(0.15,1.20,0.15))

d0=data["escr_01"][0,:,:]
ax[1,1].imshow(d0, extent=[x[0], x[-1], y[0], y[-1]],
               origin="lower",
               cmap="viridis", norm=LogNorm(np.min(d0),10))

ax[1,1].streamplot(x,y,
    data["mag_field_x"][0,:,:],data["mag_field_y"][0,:,:],
    density=0.5,color='white')

ax[1,1].set_xlabel(r'$\mathbf{x}$', fontweight='bold')
ax[1,1].set_ylabel(r'$\mathbf{y}$', fontweight='bold')
ax[1,1].set_title('t=0.1,streaming', fontweight='bold')
cb = fig.colorbar(ax[1,1].images[-1], ax=ax[1,1], fraction=0.046, pad=0.04)
cb.set_label(r'$\mathbf{E_c}$', fontweight='bold', labelpad=6)  # bold E with subscript c

ax[1,1].set_ylim(np.min(y),np.max(y))

plt.tight_layout()
for a in (ax[0,0], ax[0,1], ax[1,0], ax[1,1]):
    for s in a.spines.values(): s.set_linewidth(3)

# Bold border around the whole figure
fig.add_artist(mpatches.Rectangle((0.005, 0.005), 0.99, 0.99,
                                  transform=fig.transFigure, fill=False, lw=4, ec='black'))
plt.savefig(r'{}.png'.format(data_to_plot),dpi=720)