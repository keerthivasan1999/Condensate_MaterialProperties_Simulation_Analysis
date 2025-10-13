### Description

These files can be used to run the **initial equilibrium (P = 0 atm)** simulations to obtain the protein's natural dense-phase concentration, followed by **Langevin dynamics (LD)** simulations to study the dynamics of the condensed phase. From these simulations, one can compute **intermittent contact lifetimes**, **E2E vector relaxation times**, **radius of gyration**, and **diffusion coefficients**, for which the analysis codes are provided separately in the GitHub repository.

The example system provided here corresponds to **nSCD = 1** sequence of length 250 (120 chains) simulated using the **Kapcha–Rossky scale**.

See the file `nSCD1CL250.dat` for the nSCD = 1 (N = 250) sequence.

---

### File Details

| File | Description |
|------|--------------|
| `400.0x400.0x400.0_box.gsd` | Initial coordinate file for NPT run |
| `prod.gsd` | Final coordinate file after NPT run |
| `eq_pzero.py` | HOOMD-blue input script for NPT run |
| `hoomd_in.py` | HOOMD-blue input script for NVT run |
| `nSCD1CL250.dat` | Sequence file for nSCD = 1 (N = 250) |

---

### Simulation Workflow

1. **NPT (P = 0 atm) Equilibration Run**  
   The system is first equilibrated at *P = 0 atm* to reach the natural dense-phase concentration of the protein chains.
   Input script: eq_pzero.py
   
   Input coordinate file: 400.0x400.0x400.0_box.gsd
   
   Output coordinate file: prod.gsd

3. **Langevin Dynamics (NVT) Run**  
   The simulation is then switched to Langevin dynamics to model the condensed phase at the desired temperature (here, 300 K).
   Input script: hoomd_in.py
   
   Input coordinate file: prod.gsd
   
   Output coordinate file: long_dumps.gsd, short_dumps.gsd, veryshort_dumps.gsd

---

### Running the Simulation

Execute the simulations using HOOMD-blue (version 2.x) as:
python eq_pzero.py

python hoomd_in.py

