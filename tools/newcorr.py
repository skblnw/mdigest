import mdigest
from mdigest.core.parsetrajectory import MDS
from mdigest.core.correlation import DynCorr
import mdigest.core.savedata as sd
import os

topology_path = 'new_ionized.prot.pdb'
trajectory_paths = ['eq_set_9.xtc', 'neq_set_9.xtc']  # List of trajectories
output_directory = 'corr'

mds = MDS()
mds.set_num_replicas(1)
mds.load_system(topology_path, trajectory_paths)
mds.align_traj(selection='name CA')
mds.set_selection('protein and name CA', 'protein')
mds.stride_trajectory(initial=0, final=-1, step=1)

dyncorr = DynCorr(mds)
dyncorr.parse_dynamics(scale=True, normalize=True, LMI='gaussian', MI='None', DCC=False, PCC=True)

savedir = os.path.join(output_directory, 'mdigest_results')
if not os.path.exists(savedir):
    os.makedirs(savedir)

# Replace 'dyncorr_results' with an appropriate file name
dyncorr.save_class(file_name_root=os.path.join(savedir, 'dyncorr_results'))
