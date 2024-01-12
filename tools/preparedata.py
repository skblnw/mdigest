"""#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# @author: kevin"""

import MDAnalysis as mda
from MDAnalysis.coordinates.XTC import XTCWriter
import random
import mdigest
from mdigest.core.parsetrajectory import MDS
from mdigest.core.correlation import DynCorr
import mdigest.core.savedata as sd
import os

class DynCorrExtractor:
    """General purpose class handling computation and extraction of different correlation metrics from atomic displacements sampled over MD trajectories."""

    def __init__(self, topology_file, trajectory_file, dividing_frame_index, num_selected_frames, num_samples=10):
        # Initialize parameters
        self.topology_file = topology_file
        self.trajectory_file = trajectory_file
        self.dividing_frame_index = dividing_frame_index
        self.num_selected_frames = num_selected_frames

        # Load the trajectory using MDAnalysis
        self.u = self.mda.Universe(topology_file, trajectory_file)

        # Initialize variables for new trajectories
        self.num_frames_eq = dividing_frame_index + 1
        self.num_frames_neq = len(self.u.trajectory) - self.num_frames_eq
        self.random_indices_eq = [self.random.sample(range(self.num_frames_eq), num_selected_frames) for _ in range(num_samples)]
        self.random_indices_neq = [self.random.sample(range(self.num_frames_eq, len(self.u.trajectory)), num_selected_frames) for _ in range(num_samples)]
        
        # Log file setup
        self.log_file_name = 'frame_indices_log.txt'
        self.traj_directory = 'traj'
        self.output_directory = 'corr'

    def write_trajectory_and_log(self, trajdir, indices, segment, set_number):
        file_name = f'{trajdir}/{segment}_set_{set_number}.xtc'
        with self.XTCWriter(file_name, n_atoms=self.u.atoms.n_atoms) as W:
            for index in sorted(indices):
                self.u.trajectory[index]
                W.write(self.u)
        with open(self.log_file_name, 'a') as log_file:
            log_file.write(f'Frames used for {file_name}: {sorted(indices)}\n')

    def create_new_trajectories(self):
        with open(self.log_file_name, 'w') as log_file:
            log_file.write("Log of frames used for each trajectory\n")
            log_file.write("=====================================\n\n")
        
        trajdir = self.traj_directory
        if not self.os.path.exists(trajdir):
            self.os.makedirs(trajdir)

        for i, indices in enumerate(self.random_indices_eq):
            self.write_trajectory_and_log(trajdir, indices, 'eq', i + 1)
            self.perform_computation('eq', i + 1)

        for i, indices in enumerate(self.random_indices_neq):
            self.write_trajectory_and_log(trajdir, indices, 'neq', i + 1)
            self.perform_computation('neq', i + 1)

        print("New trajectories for eq and neq segments have been created, and the frames used are logged.")

    def perform_computation(self, segment, set_number):
        savedir = self.os.path.join(self.output_directory, f'{segment}_set_{set_number}')
        if not self.os.path.exists(savedir):
            self.os.makedirs(savedir)

        mds = self.MDS()
        mds.set_num_replicas(1)
        mds.load_system(self.topology_file, f'{segment}_set_{set_number}.xtc')
        mds.align_traj(selection='name CA')
        mds.set_selection('protein and name CA', 'protein')
        mds.stride_trajectory(initial=0, final=-1, step=1)

        dyncorr = self.DynCorr(mds)
        dyncorr.parse_dynamics(scale=True, normalize=True, LMI='gaussian', MI='None', DCC=False, PCC=True)
        dyncorr.save_class(file_name_root=savedir)

        # print("Dynamic correlation computations are complete and results are saved.")
