import h5py
import matplotlib.pyplot as plt
import seaborn as sns

def extract_matrix(hdf5_file_path, protein1_interval, protein2_interval, output_text_file_path):
    def plot_heatmap(matrix):
        plt.figure(figsize=(2,8))
        sns.heatmap(matrix, annot=False, cmap='viridis')
        plt.xlabel("Residues")
        plt.ylabel("Residues")
        # plt.colorbar(label="Distance")
        # plt.show()
        plt.tight_layout()
        plt.savefig('test.pdf')

    with h5py.File(hdf5_file_path, 'r') as file:
        def recursively_read_hdf5_group(hdf_group, indent_level=0):
            for key in hdf_group.keys():
                print(key)
                item = hdf_group[key]

                if isinstance(item, h5py.Dataset):
                    # Convert dataset to numpy array and write to file
                    matrix = item[()]
                    print(matrix.shape)
                    submatrix = matrix[protein1_interval[0]:protein1_interval[1], protein2_interval[0]:protein2_interval[1]]
                    print(submatrix.shape)
                    plot_heatmap(submatrix)

                    # Save the extracted submatrix to a text file
                    with open(output_text_file_path, 'w') as output_file:
                        all_elements = ' '.join(map(str, submatrix.ravel()))
                        output_file.write(all_elements + '\n')

                elif isinstance(item, h5py.Group):
                    # Recursively read group
                    recursively_read_hdf5_group(item, indent_level + 1)

        recursively_read_hdf5_group(file)

# Example usage
hdf5_file_path = 'corr/eq_set_1/dyncorr_results_gcc_allreplicas.h5'  # Replace with your HDF5 file path
protein1_interval = (0, 179)  # Replace with the residue interval for protein 1
protein2_interval = (375, 385)  # Replace with the residue interval for protein 2
output_text_file_path = 'output_matrix.txt'  # Replace with your desired output text file path

extract_matrix(hdf5_file_path, protein1_interval, protein2_interval, output_text_file_path)
