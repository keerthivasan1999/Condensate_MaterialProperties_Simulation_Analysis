### This Python code calculates the following:

1) The autocorrelation functions (ACFs) of the **Rouse normal modes** for polymer chains from HOOMD-blue trajectories.  
2) The amplitudes of the Rouse modes at time zero.

The input GSD file is assumed to be organized so that each linear chain of length M is ordered consecutively. For example, the beads in chain 0 are 0, ..., M-1, the beads in chain 1 are M, ..., 2M-1, etc.

The code outputs:
* `<outfile>` - A text file containing the autocorrelation functions of all Rouse modes over the sampled timesteps.  
  The first column corresponds to the time (in simulation timesteps), and the subsequent columns correspond to the normalized autocorrelation functions of modes X0, X1, X2, ...  
* `<ampfile>` - A text file containing the initial amplitudes of each Rouse mode.


### Usage

Run the code as:

python normal_modes.py -i <trajectory.gsd> -o <outfile.txt> -a <ampfile.txt> -t <sample_steps> -M <chain_length>

---
This code was was previously used in the study:

> **Sundaravadivelu Devarajan, D., Rekhi, S., Nikoubashman, A., Kim, Y.C., Howard, M.P., and Mittal, J.**  
> _Effect of Charge Distribution on the Dynamics of Polyampholytic Disordered Proteins_,  
> **Macromolecules** (2022).  
> [https://pubs.acs.org/doi/full/10.1021/acs.macromol.2c01390](https://pubs.acs.org/doi/full/10.1021/acs.macromol.2c01390)

The current version is provided to support the analyses presented in our following manuscript on E–K sequence condensate dynamics:

> **Muthukumar, K., Sundaravadivelu Devarajan, D., Kim, Y.C., and Mittal, J.**  
> _Sticky Interactions Govern Sequence-Dependent Dynamics in Biomolecular Condensates_,  
> **Proceedings of the National Academy of Sciences** (2026).  
> [https://www.pnas.org/doi/10.1073/pnas.2518384122](https://www.pnas.org/doi/10.1073/pnas.2518384122)

---
