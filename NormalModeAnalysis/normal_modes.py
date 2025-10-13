"""Compute autocorrelation function for Rouse modes.

The GSD file is assumed to be organized so that each linear chain of length M
is ordered consecutively. For example, the beads in chain 0 are 0, ..., M-1, the
beads in chain 1 are M, ..., 2M-1, etc.

"""
import argparse
import gsd.hoomd
import numpy as np
import scipy.fft

# parse arguments
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('-i', dest='infile', type=str, required=True, help='GSD filename to analyze')
parser.add_argument('-o', dest='outfile', type=str, required=True, help='Output filename')
parser.add_argument('-a', dest='ampfile', type=str, required=True, help='Amplitude filename')
parser.add_argument('-t', dest='sample_steps', type=int, required=True, help='Number of timesteps to compute ACF')
parser.add_argument('-M', dest='M', type=int, required=True, help='Number of beads in linear chain')
args = parser.parse_args()
sample_steps = args.sample_steps
M = args.M

# load trajectory
traj = None
timesteps = None
with gsd.hoomd.open(args.infile) as t:
    for i,s in enumerate(t):
        if traj is None:
            traj = np.empty((len(t),s.particles.N,3), dtype=np.float64)
            timesteps = np.empty(len(t), dtype=np.int64)
        traj[i] = s.particles.position + s.configuration.box[:3]*s.particles.image
        timesteps[i] = s.configuration.step

# assume linear chains of same length and compute Rouse modes by FFT
# Orthogonal DCT-II is the same as the definitions used in:
#   - https://doi.org/10.1063/1.474934
#   - https://doi.org/10.1039/c5sm00754b
traj = np.reshape(traj, (traj.shape[0],traj.shape[1]//M,M,3))
x = scipy.fft.dct(traj, type=2, norm='ortho', axis=-2)

# determine sampling window based on frame spacing in trajectory
dt_steps = timesteps[1:]-timesteps[:-1]
if not np.all(np.isclose(dt_steps,dt_steps[0])):
    raise ValueError('Trajectory must be evenly spaced')
dt_steps = dt_steps[0]
if sample_steps % dt_steps != 0:
    raise ValueError('Sample steps must be a multiple of trajectory spacing')
window = sample_steps//dt_steps

# evaluate autocorrelation function of Rouse modes
acf = np.zeros((window+1,M),dtype=np.float64)
counts = np.zeros((window+1,M),dtype=np.int64)
for i,x_i in enumerate(x):
    j = min(i+window+1,x.shape[0])
    delta = j-i
    x_j = x[i:j]

    acf[:delta] += np.mean(np.sum(x_i*x_j,axis=-1),axis=1)
    counts[:delta] += 1

# average the result and discard the first mode that corresponds to center of mass motion
acf = np.divide(acf, counts, where=(counts > 0))
np.savetxt(args.outfile, np.column_stack((dt_steps*np.arange(window+1),acf/acf[0])), header='timestep <X0(0).X0(timestep)> <X1(0).X1(timestep)> ...')

mode_ids = np.arange(0,M)
final_amplitudes = np.transpose(np.vstack((mode_ids, acf[0])))
np.savetxt(args.ampfile, final_amplitudes, header='<X0(0)> <X1(0)> ...')