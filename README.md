# minPRBO
minPRBO (minimal partial rebiorthogonalization) is a variant of the partial reorthogonalization strategy of Horst D. Simon that is particulary suited for the kinds of non-Hermitian Hamiltonians that occur witihin the context of ab initio calculations of  optical spectra of solids.

The Bethe-Salpeter equation (BSE) approach is the state-of-the-art in the calculation of ab initio optical spectra of materials. One important step in the calculation is that of the matrix element of the resolvent of the effective Hamiltonian. A direct method would involve the inversion of a very large operator/matrix, that is in general non-Hermitian, which becomes computationally intensive and impractical for realistic systems.

Iterative methods are a preferred alternative as they make it possible to complete calculations within a reasonable amount of time. The Lanczos based methods are quite popular in the community becasue of their ease of implementation and computational efficiency. The Tamm-Dancoff approximation (TDA) to the BSE reduces the computational complexity even further by ignoring certain blocks of the full BSE Hamiltonian which correspond to negative energy electron-hole pairs, rendering the problem Hermitian, which can be readily solved by the Hermitian Lanczos algorithm.

# IMFP
In certain material systems like silicon with strong electron-hole interactions, the TDA has been known to be inadequate in accurately describing some quantities like the energy loss function (ELF). To get ELF spectra that closely matches the experimental one, it is crucial to solve the full BSE Hamiltonian. This non-Hermitian Hamiltonian can be solved with the generalized minimal residual (GMRES) algorithm, a well-known highly accurate and stable algorithm that works for Hermitian and non-Hermitian operators. 
