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

plt.rcParams.update({
    'figure.dpi': 720, 'savefig.dpi': 720,
    'axes.linewidth': 2.5,
    'font.size': 12,
    'xtick.major.size': 5, 'ytick.major.size': 5,
    'xtick.major.width': 1.0, 'ytick.major.width': 1.0,
    'xtick.direction': 'out', 'ytick.direction': 'out',
    'legend.frameon': True, 'legend.edgecolor': 'black', 'legend.framealpha': 1.0,
    'figure.figsize': (8, 6)
})


data_to_plot = 'escr_01'

file = '/home/vinodvp/simdir/piernik/runs/test2/scr1_tst_0000.h5'
data, cell_dims, origin, spacing = load_and_stitch_data(file)
N, dx, x0 = cell_dims[0], spacing[0], origin[0]
xe = x0 + np.arange(N+1)*dx                   # assume x0 is left edge
x  = 0.5*(xe[:-1] + xe[1:])                   # centers
y=data[data_to_plot][0,0,:]
plt.plot(x,y,label='t=0.0',linewidth=1.5,color='k')

file = '/home/vinodvp/simdir/piernik/runs/test2/scr1_tst_0001.h5'
data, cell_dims, origin, spacing = load_and_stitch_data(file)
y=data[data_to_plot][0,0,:]
plt.plot(x,y,label='t=0.06',linewidth=1.5,color='r')

t = 0.06

# xm calculation
xm = np.sqrt((1+4/3*t)**2 + 8/3*t - 1)


y_analytical = np.full_like(x, 2 + 4/3*t - xm)

mask = (x < -xm) | (x > xm)
y_analytical[mask] = 2 + 4/3*t - np.abs(x[mask])
plt.plot(x,y_analytical,label='t=0.06 (Analytical)',linewidth=1.5,color='b', linestyle='dashed')

L1 = np.mean(abs(y-y_analytical))
L2 =np.sqrt( np.mean((y-y_analytical)**2))

plt.xlabel(r'$\mathbf{x}$', fontweight='bold', labelpad=6)
plt.ylabel(r'$\mathbf{E_c}$', fontweight='bold', labelpad=6)
plt.xlim(-1.0,1.0)
plt.ylim(1.0,2.0)
plt.xticks([-1,-0.5,0,0.5,1])
plt.tight_layout()
plt.legend()
ax = plt.gca()

msg = (
    rf"$L_1 = {L1:.3e}$"     "\n"
    rf"$L_2 = {L2:.3e}$"
)

ax.text(
    0.02, 0.98, msg,
    transform=ax.transAxes,
    ha="left", va="top",
    fontsize=11,
    bbox=dict(boxstyle="round,pad=0.30",
              facecolor="white", edgecolor="black",
              lw=0.8, alpha=0.9),
    zorder=5,
)
plt.savefig(r'{}.png'.format(data_to_plot),dpi=720)












