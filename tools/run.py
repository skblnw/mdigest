from preparedata import DynCorrExtractor

# Example usage
if __name__ == "__main__":
    # Replace with your actual file paths and parameters
    topology_file = 'new_ionized.prot.pdb'
    trajectory_file = 't2.prot.pf1000ps.xtc'
    dividing_frame_index = 500
    num_selected_frames = 10

    tool = DynCorrExtractor(topology_file, trajectory_file, dividing_frame_index, num_selected_frames)
    tool.create_new_trajectories()
