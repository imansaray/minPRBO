# minPRBO
minPRBO (minimal partial rebiorthogonalization) is a variant of the partial reorthogonalization strategy of Horst D. Simon that is particulary suited for the kinds of non-Hermitian Hamiltonians that occur witihin the context of ab initio calculations of  optical spectra of solids.

The Bethe-Salpeter equation (BSE) approach is the state-of-the-art in the calculation of ab initio optical spectra of materials. One important step in the calculation is that of the matrix element of the resolvent of the effective Hamiltonian. A direct method would involve the inversion of a very large operator/matrix, that is in general non-Hermitian, which becomes computationally intensive and impractical for realistic systems.

Iterative methods are a preferred alternative as they make it possible to complete calculations within a reasonable amount of time. The Lanczos based methods are quite popular in the community becasue of their ease of implementation and computational efficiency. The Tamm-Dancoff approximation (TDA) to the BSE reduces the computational complexity even further by ignoring certain blocks of the full BSE Hamiltonian which correspond to negative energy electron-hole pairs, rendering the problem Hermitian, which can be readily solved by the Hermitian Lanczos algorithm.

# minPRBO

Certain material systems like silicon 
