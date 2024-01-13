from preparedata import DynCorrExtractor

def read_params(file_path):
    params = {}
    with open(file_path, 'r') as file:
        for line in file:
            key, value = [x.strip() for x in line.strip().split('=')]
            # Check if the value is a tuple
            if value.startswith('(') and value.endswith(')'):
                value = tuple(map(int, value[1:-1].split(',')))  # Convert to tuple of integers
            params[key] = value
    return params

# Example usage
if __name__ == "__main__":
    params = read_params('parameters.txt')  # Reading parameters from the file

    topology_file = params.get('topology_file')
    trajectory_file = params.get('trajectory_file')
    dividing_frame_index = int(params.get('dividing_frame_index'))
    num_selected_frames = int(params.get('num_selected_frames'))
    protein1_interval = params.get('protein1_interval')
    protein2_interval = params.get('protein2_interval')

    tool = DynCorrExtractor(topology_file, trajectory_file, dividing_frame_index, num_selected_frames, protein1_interval, protein2_interval)
    tool.create_new_trajectories()
