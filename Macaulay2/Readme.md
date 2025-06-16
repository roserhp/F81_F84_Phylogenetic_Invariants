* Tripod_propositions.m2: this file contains all computations regarding tripods:
  * Section 3: Propositions 3.2, 3.3, 3.4
  * Section 7: Propositions 7.1, 7.2
* F81Quartet_propositions.m2: this file contains all computations regarding quartets under the F81 model:
  * Section 4: Proposition 4.7
  * Section 5: Propositions 5.4, 5.5
  * Section 6: Propositions 6.4, 6.5, Theorem 6.7
  * Section 7: Theorem 7.6, Corollaries 7.8, 7.10
* F84Quartet_propositions.m2: this file contains all computations regarding quartets under the F84 model:
  * Section 4: Proposition 4.6
  * Section 5: Proposition 5.2
  * Section 7: Theorem 7.3, Corollaries 7.5, 7.10

For the sake of readability, some required computations needed for the previous files can be found in file Tripod.m2 and Quartet.m2. All imported txt files have been generated with Tripod.m2 and Quartet.m2 are available in this repository.

* Tripod.m2: this file contains the following computations regarding tripods:
  * No-evolution point $\rho=\psi(Id,Id,Id)$ in probability coordinates $p_{i_1i_2i_3}$
  * No-evolution point $\rho=\varphi(Id,Id,Id)$ in coordinates $\tilde{p}_{i_1i_2i_3}$
  * General point $p=\varphi(D_1,D_2,D_3)$ in coordinates $\tilde{p}_{i_1i_2i_3}$ for TN93
  * Vanishing ideal in coordinates $\tilde{p}_{i_1i_2i_3}$ for TN93
  * Complete intersection in coordinates $\tilde{p}_{i_1i_2i_3}$ for TN93

* Quartet.m2: this file contains the following computations regarding quartets:
  * No-evolution point $\rho=\psi(Id,Id,Id,Id,Id)$ and point $q=\psi(Id,Id,Id,Id,M)$ for TN93 in probability coordinates $p_{i_1i_2i_3i_4}$
  * Point $q=\varphi(Id,Id,Id,Id,D)$ in coordinates $\tilde{p}_{i_1i_2i_3i_4}$ for TN93
  * General point $p=\varphi(D_1,D_2,D_3,D_4,D)$ in coordinates $\tilde{p}_{i_1i_2i_3i_4}$ for TN93
  * Flattening matrix 12|34 for a general point $p=\varphi(D_1,D_2,D_3,D_4,D)$ in coordinates $\tilde{p}_{i_1i_2i_3i_4}$ for TN93
  * Complete intersection in coordinates $\tilde{p}_{i_1i_2i_3i_4}$ for TN93

AUXILIARY FILES
* Files generated with Tripod.m2:
  * Tripod_NonZeroEntries.txt
  * Tripod_NonZeroPtilde.txt
  * TripodCI_ptilde.txt
  * TripodVanishingIdeal_ptilde.txt
* Files generated with Quartet.m2:
  * qtilde1234.txt
  * Flat1234_ptilde.txt
  * QuartetCI_ptilde.txt
  * QuartetIntersectionVanishingIdeals_ptilde.txt
* Macaulay2 package Permutations (only needed for older M2 versions):
  * Permutations.m2
  * Folder Permutations, with auxiliary files 

