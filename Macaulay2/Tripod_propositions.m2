-----------------------------------------------------------------------------------------------------------
-----------------------------------------------------------------------------------------------------------
--
-- Tripod_propositions.m2: this file contains all computational proofs regarding tripods:
--
-- Section 3: Propositions 3.2, 3.3, 3.4
-- Section 7: Propositions 7.1, 7.2
--
-- For the sake of readability, some required computations can be found in file Tripod.m2:
-- 1. No-evolution point rho=varphi(Id,Id,Id) in probability coordinates
-- 2. No-evolution point rho=varphi(Id,Id,Id) in ptilde coordinates
-- 3. General point p=varphi(D_1,D_2,D_3) in ptilde coordinates for TN93
-- 4. Vanishing ideal in ptilde coordinates for TN93
-- 5. Complete intersection for TN93
-- All imported txt files have been generated with Tripod.m2 and are available in this repository
--
----------------------------------------------------------------------------------------------------------- 
----------------------------------------------------------------------------------------------------------- 

-------------------------------------------------------------------
-- Setup
-------------------------------------------------------------------
restart

-- Ring declarations: 

-- variable l_(i,j) corresponds to eigenvalue j of transition matrix i
-- Since pi does not appear in the parametrization with ptilde coordinates, we omit them and work over QQ
R=QQ[l_(1,1)..l_(3,4)];

-- Retrieve non-zero coordinates of general point varphi(D_1,D_2,D_3) for TN93 in ptilde coordinates
nonZeroPtilde=value get "Auxiliary_files/Tripod_NonZeroPtilde.txt";
-- Retrieve indices of non-zero coordinates
nonZeroEntries=value get "Auxiliary_files/Tripod_NonZeroEntries.txt";
-- variables p_ijk, non-zero ptilde coordinates
varp=toList apply(nonZeroEntries,i->(symbol p)_i);
Rp=QQ[varp]; 

-- Vanishing ideal and complete intersection ideal for TN93 in ptilde coordinates
I=value get "Auxiliary_files/TripodVanishingIdeal_ptilde.txt";
CI=value get "Auxiliary_files/TripodCI_ptilde.txt";

-------------------------------------------------------------------
-- Propositions 3.2 and  3.3: vanishing ideal for F84
-------------------------------------------------------------------

-- In model F84, eigenvalues 3 and 4 are equal:
F84=toList apply(1..3,i->l_(i,3)=>l_(i,4))
-- Non-zero ptilde coordinates for F84
ptildeF84=matrix{toList apply(nonZeroPtilde,i->sub(i,F84))};
-- pullback of map varphi (in ptilde coordinates)
fF84=map(R,Rp,ptildeF84);
IF84=time trim kernel fF84;
betti IF84 -- 7 linear, 3 quadrics, 12 cubics
netList IF84_*

-------------------------------------------------------------------
-- Proposition 3.4: vanishing ideal for F81
-------------------------------------------------------------------

-- In model F81, eigenvalues 2, 3 and 4 are equal:
F81=flatten toList apply(1..3,i->{l_(i,3)=>l_(i,4),l_(i,2)=>l_(i,4)})
-- Non-zero ptilde coordinates for F81
ptildeF81=matrix{toList apply(nonZeroPtilde,i->sub(i,F81))};
-- pullback of map varphi (in ptilde coordinates)
fF81=map(R,Rp,ptildeF81);
IF81=time trim kernel fF81;
betti IF81 -- 14 linear, 1 cubics
netList IF81_*

-------------------------------------------------------------------
-- Proposition 7.1: complete intersection for F84
-------------------------------------------------------------------

--Symmetry equations for F84 (linear span of IF84)
S84=ideal{p_(1,4,4)-p_(1,3,3),p_(4,1,4)-p_(3,1,3),p_(4,4,1)-p_(3,3,1),p_(2,4,4)-p_(2,3,3),p_(4,2,4)-p_(3,2,3),p_(4,4,2)-p_(3,3,2),p_(4,4,4)-p_(3,3,3)};

--Complete intersection for F84
CIF84=S84+ideal{IF84_7,IF84_8,IF84_9,IF84_10,IF84_14}

betti CIF84
netList CIF84_*

codim CIF84==numgens CIF84 --true
codim CIF84==codim IF84 --true
-- This defines a complete intersection (not just locally, but a global CI)
-- Note that in Proposition 7.1 we omit linear generators because we work in the space of mixtures

-- V(IF84) is an irreducible component of V(CIF84),
-- saturation of V(CIF84) w.r.t. p_111*p_222*p_444 coincides with IF84
PDF84=primaryDecomposition CIF84;
length PDF84 --11
PDF84_10==IF84 --true
saturate(CIF84,sub(ideal{p_(1,1,1)*p_(2,2,2)*p_(4,4,4)},Rp))==IF84 --true


------------------------------------------------------------------------------------------
-- Proposition 7.2: obtaining ideals for F81 and F84 from ideals of TN93 plus symmetry equations
------------------------------------------------------------------------------------------

I84=trim(I+S84);
IF84==I84 --true

CI84=trim(CI+S84);
CI84==CIF84 --true

-- Symmetry equations for F81 (linear span of IF81)
S81=S84+ideal{p_(1,4,4)-p_(1,2,2),p_(4,1,4)-p_(2,1,2),p_(4,4,1)-p_(2,2,1),p_(2,4,4)-p_(4,4,4),p_(4,2,4)-p_(4,4,4),p_(4,4,2)-p_(4,4,4),p_(4,4,4)-p_(2,2,2)};

I81=trim(I+S81);
IF81==I81 --true
