import MDAnalysis as mda
from MDAnalysis.coordinates.XTC import XTCWriter
import random

# Load the trajectory
topology_file = 'new_ionized.prot.pdb'  # Replace with your topology file
trajectory_file = 't2.prot.pf1000ps.xtc'  # Replace with your trajectory file
u = mda.Universe(topology_file, trajectory_file)

# Given frame index to divide the trajectory
dividing_frame_index = 500  # Replace with your specific frame index

# Determine the number of frames in each segment
num_frames_eq = dividing_frame_index + 1
num_frames_neq = len(u.trajectory) - num_frames_eq

# Number of frames to select from each segment for each set
num_selected_frames = 10  # Adjust as needed

# Create random sets of indices for each segment
random_indices_eq = [random.sample(range(num_frames_eq), num_selected_frames) for _ in range(10)]
random_indices_neq = [random.sample(range(num_frames_eq, len(u.trajectory)), num_selected_frames) for _ in range(10)]

# Function to write selected frames to a new trajectory and record indices
def write_trajectory_and_record(indices, segment, set_number):
    file_name = f'{segment}_set_{set_number}.xtc'
    with XTCWriter(file_name, n_atoms=u.atoms.n_atoms) as W:
        for index in sorted(indices):
            u.trajectory[index]
            W.write(u)
    # Record the indices used
    with open(f'{segment}_set_{set_number}_frames.txt', 'w') as f:
        f.write(f'Frames used for {file_name}:\n')
        f.write(', '.join(map(str, sorted(indices))))

# Create new trajectories and record frame indices for each set of indices
for i, indices in enumerate(random_indices_eq):
    write_trajectory_and_record(indices, 'eq', i + 1)

for i, indices in enumerate(random_indices_neq):
    write_trajectory_and_record(indices, 'neq', i + 1)

print("New trajectories for eq and neq segments have been created, along with records of the frames used.")
