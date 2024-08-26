# minPRBO
minPRBO (minimal partial rebiorthogonalization) is a variant of the partial reorthogonalization strategy of Horst D. Simon that is particulary suited for the kinds of non-Hermitian Hamiltonians that occur witihin the context of ab initio calculations of  optical spectra of solids.

The Bethe-Salpeter equation (BSE) approach is the state-of-the-art in the calculation of ab initio optical spectra of materials at all scales (small molecules, extended systems, solids). One important step in the calculation is that of obtaining the matrix element of the resolvent of the effective Hamiltonian. A direct method would involve the inversion of a very large operator/matrix, that is in general non-Hermitian, which becomes computationally intensive and impractical for realistic systems.

Iterative methods are a preferred alternative as they make it possible to complete calculations within a reasonable amount of time. The Lanczos based methods are quite popular becasue of their ease of implementation and computational efficiency. The Tamm-Dancoff approximation (TDA) to the BSE reduces the computational complexity even further by ignoring certain blocks of the non-Hermitian full BSE Hamiltonian which correspond to negative energy electron-hole pairs, rendering the problem Hermitian, which can be readily solved by the Hermitian Lanczos algorithm.

# IMFP
In certain material systems like silicon with strong electron-hole interactions, the TDA has been known to be inadequate in accurately describing some quantities like the energy loss function (ELF). To get ELF spectra that closely matches the experimental one, it is crucial to solve the full BSE Hamiltonian. This non-Hermitian Hamiltonian can be solved with the generalized minimal residual (GMRES) algorithm, a well-known highly accurate and stable algorithm that works for Hermitian and non-Hermitian operators. 

GMRES obtains the dielectric function and subsequently the ELF by solving a new linear system for each required energy. This becomes a challenge when we want to obtain ab initio ELF's that cover a large energy range, like say 0 to 100 eV, and important regime in the inelastic mean free path (IMFP) of materials. Current models of obtaining the IMFP of materials fail in this energy regime, and various experimental approaches give conflicting data. We believe that ab initio methods, by making the fewest assumptions and their demonstrated utility in other calculations of material properties can shed light on this issue.

In order to be able to calculate the IMFP requires the ELF to be sampled reasonably well at different momentum transfer over the desired energy range (0 to 100 eV). We seek an iterative approach that is as versatile as the Hermitian Lanczos algorithm for the TDA and as accurate as GMRES for the non-Hermitian full BSE Hamiltonian. The non-Hermtian Lanczos algorithm fits this bill and has been used before to calculate the optical spectra of a small molecule Benzene.

In that implementation of partial reorthogonalization applied to the non-Hermitian case, when it was deemed necessary to carry out a reorthogonalization of the Lanczos vectors, a full reorthogonalization step was carried out with respect to all the previous Lanczos vectors, then a second full reorthogonalization step of the Lanczos vectors was carried out in the subsequent iteration to maintain semi-orthogonality. This is considerably cheaper than carrying out a full reorthogonalization of the Lanczos vectors at each iteration, which will bring the computational costs on par with a GMRES calculation.

We contend that this approach to partial reorthogonalization for the non-Hermitian case can still be improved upon, especially given the application to the calculatioin of the IMFP that we have in mind.
We propose to do the minimum number of rebiorthogonalizations of the Lanczos vectors whilst still maintaining their semi-orthogonality. This is achieved by only rebiorthogonalizing the Lanczos vectors to those with which they have lost semi-orthogonality. This is shown to be sufficient to maintain semi-orthogonality and also achieve calculated spectra that are as accurate as the GMRES obtained spectra, with errors on the order of less than one percent.

The sample programs provided carry out this minimal partial rebiorthogonalization (minPRBO) on model Hamiltonians and calcualtes the dielectric function and compares them to the exact solution.
