import h5py

def extract_distance_matrix(hdf5_file_path, protein1_interval, protein2_interval, output_text_file_path):
    with h5py.File(hdf5_file_path, 'r') as file:
        # Assuming there's only one dataset, get its name
        dataset_name = list(file.keys())[0]

        # Access the symmetric distance matrix dataset
        distance_matrix = file[dataset_name][()]

        # Extract the relevant submatrix
        # Assuming protein1_interval and protein2_interval are tuples (start, end)
        submatrix = distance_matrix[protein1_interval[0]:protein1_interval[1], protein2_interval[0]:protein2_interval[1]]

        # Save the extracted submatrix to a text file
        with open(output_text_file_path, 'w') as output_file:
            for row in submatrix:
                row_str = ' '.join(map(str, row))
                output_file.write(row_str + '\n')

# Example usage
hdf5_file_path = 'corr/mdigest_results/dyncorr_results_pcc_allreplicas.h5'  # Replace with your HDF5 file path
protein1_interval = (0, 100)  # Replace with the residue interval for protein 1
protein2_interval = (100, 200)  # Replace with the residue interval for protein 2
output_text_file_path = 'distance_matrix_output.txt'  # Replace with your desired output text file path

extract_distance_matrix(hdf5_file_path, protein1_interval, protein2_interval, output_text_file_path)
