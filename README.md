# BSE
The Bethe-Salpeter equation (BSE) approach is the state-of-the-art in the calculation of ab initio optical spectra of materials in the solid state [1]. One important step in the calculation is that of obtaining the matrix element of the resolvent of the effective two-particle Hamiltonian. A direct method would involve the inversion of a very large operator/matrix that is in general non-Hermitian, which becomes computationally intensive and impractical for realistic systems.

Iterative methods are a preferred alternative as they make it possible to complete calculations within a reasonable amount of time. The Lanczos based methods are quite popular becasue of their ease of implementation and computational efficiency. The Tamm-Dancoff approximation (TDA) to the BSE reduces the computational complexity even further by ignoring certain blocks of the non-Hermitian full BSE Hamiltonian which correspond to negative energy electron-hole pairs, rendering the problem Hermitian, which can be readily solved by the Hermitian Lanczos algorithm.

# ELF and IMFP
In certain material systems with strong electron-hole interactions like silicon, the TDA has been known to be inadequate in accurately describing some quantities like the energy loss function (ELF). To get ELF spectra that closely matches the experimental one it is necessary to solve the full BSE Hamiltonian. One way to solve the full BSE Hamiltonian is by recasting it as a generalized eigenvalue problem and redefining the inner product of vectors, after which a Hermitian-like Lanczos algorithm with this new inner product can be used. This method however only works if the block structure of the full BSE is maintained. If other excitations like phonons are included, then the nice block structure no longer exists and such a method will fail. The generalized minimal residual (GMRES) algorithm, a well-known, highly accurate and stable algorithm works for Hermitian and non-Hermitian operators irrespective of any underlying structure. 

GMRES obtains the dielectric function and subsequently the ELF by solving a new linear system for each required energy. This becomes a challenge when we want to obtain ab initio ELF's that cover a large energy range, like say 0 to 100 eV, an important regime in the inelastic mean free path (IMFP) of materials. Current theoretical models used to calculate the IMFP of materials fail in this energy regime, and various experimental approaches give conflicting data. We believe that ab initio methods, by making the fewest assumptions and with their demonstrated utility in other calculations of material properties can shed light on this issue.

In order to be able to calculate the IMFP, one requires the ELF to be sampled reasonably well at different momentum transfer over the desired energy range (0 to 100 eV). We seek an iterative approach that is as versatile as the Hermitian Lanczos algorithm for the TDA and as accurate as GMRES for the non-Hermitian full BSE Hamiltonian. The non-Hermtian Lanczos algorithm fits this requirement and has been used before with a partial reorthogonalization strategy [2] to calculate the optical spectra of a small molecule Benzene [3,4].

# minPRBO
The non-Hermitian Lanczos algorithm requires two sets of Lanczos vectors and one speaks of rebiorthogonalization (full or partial) of these vectors to maintain stability of the algorithm. In the original implementation of the partial reorthogonalization strategy of Horst D. Simon applied to the non-Hermitian case, once it was deemed necessary to carry out a reorthogonalization of the Lanczos vectors, a full rebiorthogonalization step was carried out with respect to all the previous Lanczos vectors. A second full rebiorthogonalization step of the Lanczos vectors was carried out in the subsequent iteration to maintain semi-biorthogonality. This is considerably cheaper than carrying out a full rebiorthogonalization of the Lanczos vectors at each iteration which will bring the computational costs on par with a GMRES calculation.

We contend that this approach to partial reorthogonalization for the non-Hermitian case can still be improved upon especially given the application we have in mind, that of the calculatioin of the IMFP.
We propose to do the minimum number of rebiorthogonalizations of the Lanczos vectors whilst still maintaining their semi-biorthogonality. This is achieved by only rebiorthogonalizing the Lanczos vectors against those with which they have lost semi-biorthogonality. This is shown to be sufficient to maintain semi-biorthogonality and the calculated spectra are as accurate as the GMRES spectra, with errors on the order of less than one percent.

The sample programs provided carry out this minimal partial rebiorthogonalization (minPRBO) strategy on model Hamiltonians and on an ab initio Hamiltonian for silicon with complex energies that mimimcs the effects of phonons and calcualtes the dielectric function and compares them to the exact solution when possible and to the GMRES obtained solution.

# List of provided programs

## Model Hamiltonians
* lanczos.f90 - this contains implementations of the Hermitian Lanczos algorithm with full and partial reorthogonalization (PRO), the nonsymmetric Lanczos algorithm for real matrices, and the non Hermitian Lanczos algorithm with full rebiorthogonalization, minPRBO and partial rebiorthogonalization as originally proposed and the exact calculation of the dielectric function.
* planczos_full.f90 - this contains a parallel version of the full rebiorthogonalization strategy using ScaLAPACK and PBLAS libraries.
* planczos_minprbo.f90 -  this contains a parallel version of the minPRBO strategy using SCALAPACK and PBLAS libraries.

# Compilation
A Makefile file is provided showing the required libraries for compiling the different programs. All programs can be compiled by typing 

```
Make
```

# Execution
## serial programs
All serial programs can be run in the usual manner by running the executable file, for example
```
./lanczos
```
## parallel programs
For the parallel programs, mpirun or mpiexec can be used to run the compiled programs with the number of processors specified. For example
```
mpirun -n 4 planczos_full
```

# Example data
Example data obtained by the different strategies and the parameters used are in the EXAMPLES folder. A simple GNUPLOT program to plot the data and example plots are also provided.

# References
- [1] Onida, Giovanni, Lucia Reining, and Angel Rubio. "Electronic excitations: density-functional versus many-body Green’s-function approaches." Reviews of modern physics 74.2 (2002): 601.
- [2] Simon, Horst D. "The Lanczos algorithm with partial reorthogonalization." Mathematics of computation 42.165 (1984): 115-142.
- [3] Bai, Zhaojun, David Day, and Qiang Ye. "ABLE: an adaptive block Lanczos method for non-Hermitian eigenvalue problems." SIAM Journal on Matrix Analysis and Applications 20.4 (1999): 1060-1082.
- [4] Walker, Brent, et al. "Efficient Approach to Time-Dependent Density-Functional Perturbation Theory<? format?> for Optical Spectroscopy." Physical review letters 96.11 (2006): 113001.
