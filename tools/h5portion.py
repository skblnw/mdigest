import h5py

def extract_matrix(hdf5_file_path, protein1_interval, protein2_interval, output_text_file_path):
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

                    # Save the extracted submatrix to a text file
                    with open(output_text_file_path, 'w') as output_file:
                        for row in submatrix:
                            row_str = ' '.join(map(str, row))
                            output_file.write(row_str + '\n')
                elif isinstance(item, h5py.Group):
                    # Recursively read group
                    recursively_read_hdf5_group(item, indent_level + 1)
        return submatrix

        submatrix = recursively_read_hdf5_group(file)

    return submatrix

# Example usage
hdf5_file_path = 'corr/mdigest_results/dyncorr_results_gcc_allreplicas.h5'  # Replace with your HDF5 file path
protein1_interval = (0, 374)  # Replace with the residue interval for protein 1
protein2_interval = (375, 384)  # Replace with the residue interval for protein 2
output_text_file_path = 'output_matrix.txt'  # Replace with your desired output text file path

submatrix = extract_matrix(hdf5_file_path, protein1_interval, protein2_interval, output_text_file_path)

import matplotlib.pyplot as plt
import seaborn as sns

def plot_heatmap(matrix, output_image_file_path):
    plt.figure(figsize=(10, 8))
    sns.heatmap(matrix, annot=False, cmap='viridis')
    plt.title("Heatmap of Extracted Submatrix")
    plt.xlabel("Residue Index Protein 2")
    plt.ylabel("Residue Index Protein 1")
    plt.colorbar(label="Distance")
    plt.show()

plot_heatmap(submatrix)
