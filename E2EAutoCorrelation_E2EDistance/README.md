### This C++ code calculates the following:

1) End-to-end distance as a function of time in a file named "E2EDistance.txt"  
2) End-to-end autocorrelation for inter-chains as a function of time (only when there are multiple chains in the system) in a file named "Inter_E2EAutoCorr.txt"  
3) End-to-end autocorrelation for intra-chains as a function of time in a file named "Intra_E2EAutoCorr.txt"

Note that the "Intra_E2EAutoCorr.txt" will have 3 columns: The first column is the timestep, the second column is just cos(theta), and the third column is the dot product of the end-to-end vectors divided by the initial dot product value at time t = 0.

---

### Requirements

* This code requires dump files named "DumpFile.TimeStep" sorted based on atom-ids and has to be in the following format: "id mol type q mass xu yu zu ix iy iz"  
* xu, yu, zu refers to unwrapped coordinates

---

### Required files

* This code comes with 3 files: (1) one .cpp file, (2) one .h header file, (3) one parameter text file  
* In the parameter file, one needs to input the information about the system and the simulation details like trajectory file name, number of chains, number of beads per chain, start time, end time, interval time, autocorrelation end time.  Autocorrelation end time is nothing but let's say one has trajectory information for 10^8 steps but knows that the correlation goes to zero in 10^6 steps, one can then specify autocorrelation end time as 10^6 and the code will terminate once it has calculated the autocorrelation value until 10^6 steps.

---

### Compile and run the C++ code as below

* If using a GCC compiler:  
    gcc -lstdc++ -lm  EndtoEndOrientation.cpp  
    ./a.out
    
* If using a ICC compiler:  
    icc EndtoEndOrientation.cpp  
    ./a.out

---

### Author and Provenance

This code was originally written by **Dinesh S. Devarajan** and was previously used in the study:

> **Sundaravadivelu Devarajan, D., Wang, J., Szała-Mendyk, B., Rekhi, S., Nikoubashman, A., Kim, Y.C., and Mittal, J.**  
> _Sequence-dependent material properties of biomolecular condensates and their relation to dilute phase conformations_,  
> **Nature Communications** (2024).  
> [https://www.nature.com/articles/s41467-024-46223-w](https://www.nature.com/articles/s41467-024-46223-w)

The current version is provided to support the analyses presented in our following manuscript on E–K sequence condensate dynamics:

> **Muthukumar, K., Sundaravadivelu Devarajan, D., Kim, Y.C., and Mittal, J.**  
> _Sticky Interactions Govern Sequence-Dependent Dynamics in Biomolecular Condensates_,  
> **bioRxiv** (2025).  
> [https://doi.org/10.1101/2025.07.09.664001](https://doi.org/10.1101/2025.07.09.664001)

Please cite these works if you use the codes provided here, either as-is or in modified form.

---

### Other details

This code does moving average over all the frames and also averages over chains if there are more than 1 chain in the system.

