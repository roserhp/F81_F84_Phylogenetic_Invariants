-----------------------------------------------------------------------------------------------------------
-----------------------------------------------------------------------------------------------------------
--
-- F81Quartet_propositions.m2: this file contains all computational proofs regarding tripods:
--
-- Section 4: 
-- Section 7: 7.6
--

-------------------------------------------------------------------
-- Setup
-------------------------------------------------------------------
restart

-- Ring declarations: 

-- p_1,p_2,p_3,p_4 are parameters representing the root distribution p=(p_1,p_2,p_3,p_4)
-- variable l_(i,j) corresponds to eigenvalue j of transition matrix i
K=frac(QQ[p_1,p_2,p_3,p_4]);
R=K[l_(1,1)..l_(5,4)];

-- Retrieve non-zero coordinates of general point varphi(D_1,D_2,D_3) for TN93 in ptilde coordinates
nonZeroPtilde=value get "1234NonZeroPtilde.txt";
-- Retrieve indices of non-zero coordinates
nonZeroEntries=value get "1234NonZeroEntries.txt";
-- variables p_ijk, non-zero ptilde coordinates
varp=toList apply(nonZeroEntries,i->(symbol p)_i);
Rp=K[varp]; 

-- Complete intersection ideal for TN93 in ptilde coordinates
CI=value get "Quartet1234CI_ptilde.txt";
betti CI --28 quadrics, 4 cubics, 14 quartics


-------------------------------------------------------------------
-- Theorem 7.6: vanishing ideal for F84
-------------------------------------------------------------------

-- In model F81, eigenvalues 2, 3 and 4 are equal:
F84=flatten toList apply(1..5,i->{l_(i,3)=>l_(i,4)})
-- Non-zero ptilde coordinates for F81
ptildeF84=matrix{toList apply(nonZeroPtilde,i->sub(i,F84))};


F=ideal{x_(1,1,4,4)-x_(1,1,2,2), -- (a)
        x_(1,4,1,4)-x_(1,2,1,2),
	x_(1,4,4,1)-x_(1,2,2,1),
	x_(4,1,1,4)-x_(2,1,1,2),
	x_(4,1,4,1)-x_(2,1,2,1),
	x_(4,4,1,1)-x_(2,2,1,1),
	x_(1,2,4,4)-x_(1,4,4,4), 
        x_(1,4,2,4)-x_(1,4,4,4),
	x_(1,4,4,2)-x_(1,4,4,4),
	x_(4,1,2,4)-x_(4,1,4,4),
	x_(4,1,4,2)-x_(4,1,4,4),
	x_(4,4,1,2)-x_(4,4,1,4),
	x_(2,1,4,4)-x_(4,1,4,4),
        x_(2,4,1,4)-x_(4,4,1,4),
	x_(2,4,4,1)-x_(4,4,4,1),
	x_(4,2,1,4)-x_(4,4,1,4),
	x_(4,2,4,1)-x_(4,4,4,1),
	x_(4,4,2,1)-x_(4,4,4,1),
        x_(1,2,2,2)-x_(1,4,4,4),
        x_(2,1,2,2)-x_(4,1,4,4),
	x_(2,2,1,2)-x_(4,4,1,4),
	x_(2,2,2,1)-x_(4,4,4,1),
	x_(4,2,4,4)-x_(2,4,4,4), -- (b)
	x_(4,4,2,4)-x_(2,4,4,4),
	x_(4,4,4,2)-x_(2,4,4,4),
	x_(3,3,4,4)-x_(4,4,3,3)} -- (c)
