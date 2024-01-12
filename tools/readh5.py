import h5py

def hdf5_to_text(hdf5_file_path, output_text_file_path):
    with h5py.File(hdf5_file_path, 'r') as hdf, open(output_text_file_path, 'w') as output_file:
        def recursively_read_hdf5_group(hdf_group, indent_level=0):
            for key in hdf_group.keys():
                output_file.write('\t' * indent_level + f'Key: {key}\n')
                item = hdf_group[key]

                if isinstance(item, h5py.Dataset):
                    # Convert dataset to numpy array and write to file
                    data = item[()]
                    output_file.write('\t' * indent_level + f'Dataset shape: {data.shape}\n')
                    output_file.write('\t' * indent_level + f'Data: {data}\n\n')
                elif isinstance(item, h5py.Group):
                    # Recursively read group
                    recursively_read_hdf5_group(item, indent_level + 1)

        recursively_read_hdf5_group(hdf)

# Example usage
hdf5_file_path = 'corr/mdigest_results/dyncorr_results_gcc_allreplicas.h5'  # Replace with your HDF5 file path
output_text_file_path = 'output.txt'  # Replace with your desired output text file path
hdf5_to_text(hdf5_file_path, output_text_file_path)
