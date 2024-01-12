import h5py

# Open the HDF5 file
with h5py.File(h5file, 'r') as file:
    # List all groups
    print("Keys: %s" % file.keys())
    a_group_key = list(file.keys())[0]

    # Get the data
    data = list(file[a_group_key])

    # Now, data contains the data stored in the first group.
    # You can process this data as needed.

h5file = 'dyncorr_results.h5'
