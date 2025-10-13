# Dense-Phase Martini 3 Simulations for Viscosity and Diffusivity

These set of files can be used to run the dense phase simulations for computing **viscosity** based on the **Green-Kubo relation**.  
The example set of files given here are for the **diblock E-K sequence (80 chains)**, simulated using the **Martini 3 model**.  
See `EKV_15.dat` for the sequence of diblock E-K.

Note that an extension of this run to longer times, as discussed in the manuscript, can be done for the purpose of computing **diffusion coefficients** at long times.

---

## File Details

- Initial coordinate file > `PRO_SOL_IONS.gro`
- Final coordinate file > `prod.gro`
- Topology file > `PRO_SOL_IONS.top`
- GROMACS input script file > `prod_Martini_NVT.mdp`

---

## Diffusion Coefficient Calculation

Command used to compute mean-squared displacement within GROMACS from which diffusion coefficient was extracted as discussed in the Methods
section of the manuscript: "gmx msd" within the GROMACS software package. This command is part of the GROMACS software package.

---

## Viscosity Calculation (Green-Kubo Relation)

For computing viscosity based on the Green-Kubo relation, we used the codes shared on a Github repository as part
of a recent Journal of Chemical Information and Modeling (JCIM) article by Prass et al., 2023 (link: https://doi.org/10.1021/acs.jcim.3c00947).
---

## Provenance

This simulation setup was originally developed by **Dinesh S. Devarajan** and was previously used in the study:

> **Sundaravadivelu Devarajan, D., Wang, J., Szała-Mendyk, B., Rekhi, S., Nikoubashman, A., Kim, Y.C., and Mittal, J.**  
> _Sequence-dependent material properties of biomolecular condensates and their relation to dilute phase conformations_,  
> **Nature Communications** (2024).  
> [https://www.nature.com/articles/s41467-024-46223-w](https://www.nature.com/articles/s41467-024-46223-w)

The current version is provided to support the **dense-phase Martini 3 simulations** presented in our following manuscript on E–K sequence condensate dynamics:

> **Muthukumar, K., Sundaravadivelu Devarajan, D., Kim, Y.C., and Mittal, J.**   
> _Sticky Interactions Govern Sequence-Dependent Dynamics in Biomolecular Condensates_,  
> **bioRxiv** (2025).  
> [https://doi.org/10.1101/2025.07.09.664001](https://doi.org/10.1101/2025.07.09.664001)

Please cite these works if you use the codes provided here, either as-is or in modified form.

