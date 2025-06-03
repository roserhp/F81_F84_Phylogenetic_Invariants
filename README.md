# Computing phylogenetic invariants for time-reversible models: from TN93 to its submodels
This repository contains files for the manuscript [Computing phylogenetic invariants for time-reversible models: from TN93 to its submodels](https://arxiv.org/abs/2505.20526) (arXiv:2505.20526) by Marta Casanellas, Jennifer Garbett, Roser Homs, Annachiara Korchmaros, Niharika Chakrabarty Paul.

## File descriptions

The repository contains two folders with code in the following to software systems:

* Macaulay2:
  * Tripod_propositions.m2: This file contains computational proofs for Propositions 3.2, 3.3, 3.4, 7.1 and 7.2.
  * Tripod.m2: This file contains necessary computations and generates the data files for Tripod_propositions.m2:
     * NonZeroPtilde.txt
     * NonZeroEntries.txt
     * TripodCI_ptilde.txt
     * TripodVanishingIdeal.txt
* Sage:
  * F84_nep.sage: This file contains the code for computations that are used in the proof of Theorem 7.3.
  * F84QuartetIdeals.ipynb: This file contains code to verify Corollary 7.5 computationally by computing relevant ideals.  
  * F81_nep.sage: This file contains the code for computations that are used in the proof of Theorem 7.6.
  * F81QuartetIdeals.ipynb: This file contains code to verify Theorem 7.6 and Corollary 7.8 computationally by computing relevant ideals. 
  
Under construction, completion expected by June 15th.
