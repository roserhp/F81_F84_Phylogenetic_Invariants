-----------------------------------------------------------------------------------------------------------
-----------------------------------------------------------------------------------------------------------
--
-- F84Quartet_propositions.m2: this file contains computations regarding quartets evolving under F84:
--
-- 0. Setup
-- 1. Parametrization 
--
-- Computational proofs (unless stated otherwise) of the following results:
--
-- 2. Proposition 4.6 (sanity check): the listed equations are indeed linear model equations
-- 3. Proposition 5.2: the listed equations are linear topology equations 
-- 4. Theorem 7.3: complete intersection
-- 5. Corollary 7.5: rank constraints in flattening matrices. See computations in Sage.
-- 6. Corollary 7.10 (partial proof): complete intersection of TN93 intersected with symmetry equations 
--
-------------------------------------------------------------------

-------------------------------------------------------------------
-- 0. Setup
-------------------------------------------------------------------
restart

-- Ring declaration: eigenvalues
-- p_1,p_2,p_3,p_4 are parameters representing the root distribution p=(p_1,p_2,p_3,p_4)
-- variable l_(i,j) corresponds to eigenvalue j of transition matrix i
K=frac(QQ[p_1,p_2,p_3,p_4]);
R=K[l_(1,1)..l_(5,4)];

-- ptilde coordinates for point q=varphi(Id,Id,Id,Id,D_5) for tree 12|34 for TN93
qtilde12=value get "Auxiliary_files/qtilde1234.txt";

-- permuting coordinates of q to obtain coordinates for the other tree topologies
st={1,2,3,4};
S=sort elements (set st)^**4/splice/splice;

load "Auxiliary_files/Permutations.m2"
--Package Permutations.m2 is available in new M2 versions

sigma13=permutation{1,3,2,4}
sigma14=permutation{1,4,3,2}
S13=toList apply(S,i->toSequence(sigma13*i));
S14=toList apply(S,i->toSequence(sigma14*i));
-- ptilde coordinates for point q=varphi(Id,Id,Id,Id,D_5) for trees 13|24 and 14|23 for TN93
qtilde13=transpose matrix{apply(S13,i->qtilde12_(position(S,j->j==i),0))};
qtilde14=transpose matrix{apply(S14,i->qtilde12_(position(S,j->j==i),0))};

--ptilde coordinates for general point p=varphi(D_1,D_2,D_3,D_4,D_5)
D1=diagonalMatrix toList(l_(1,1)..l_(1,4))
D2=diagonalMatrix toList(l_(2,1)..l_(2,4))
D3=diagonalMatrix toList(l_(3,1)..l_(3,4))
D4=diagonalMatrix toList(l_(4,1)..l_(4,4))
D=D1**D2**D3**D4;

ptilde12=D*qtilde12;
ptilde13=D*qtilde13;
ptilde14=D*qtilde14;

--non-zero coordinates for general point (at least for some tree topology)
nonZeroEntries=sort unique join(S_(positions(flatten entries ptilde12,i->i!=0)),
S_(positions(flatten entries ptilde13,i->i!=0)),
S_(positions(flatten entries ptilde14,i->i!=0)))

-- Ring declaration: ptildes
-- variables p_ijk, non-zero ptilde coordinates
varp=toList apply(nonZeroEntries,i->(symbol p)_i);
Rp=K[varp]; 

-----------------------------------------------------------------------
-- 1. Parametrization 
-----------------------------------------------------------------------

-- In model F84, eigenvalues 3 and 4 are equal:
F84=flatten toList apply(1..5,i->{l_(i,3)=>l_(i,4)})

-- Non-zero ptilde coordinates for F81 for tree topology 12|34
ptilde12F84=flatten entries sub(ptilde12^(apply(nonZeroEntries,i->position(S,j->j==i))),F84)
ptilde13F84=flatten entries sub(ptilde13^(apply(nonZeroEntries,i->position(S,j->j==i))),F84)
ptilde14F84=flatten entries sub(ptilde14^(apply(nonZeroEntries,i->position(S,j->j==i))),F84)

--Parametrization map
f12=map(R,Rp,ptilde12F84);
f13=map(R,Rp,ptilde13F84);
f14=map(R,Rp,ptilde14F84);

------------------------------------------------------------------
-------------------------------------------------------------------
-- 2. Proposition 4.7 symmetry equations for F81
-------------------------------------------------------------------
-------------------------------------------------------------------

-- Symmetry equations for F84

F=trim ideal{p_(1,1,4,4)-p_(1,1,3,3), -- (a.1)
        p_(1,4,1,4)-p_(1,3,1,3),
	p_(1,4,4,1)-p_(1,3,3,1),
	p_(4,1,1,4)-p_(3,1,1,3),
	p_(4,1,4,1)-p_(3,1,3,1),
	p_(4,4,1,1)-p_(3,3,1,1),
	p_(1,2,4,4)-p_(1,2,3,3), -- (a.2) 
        p_(1,4,2,4)-p_(1,3,2,3),
	p_(1,4,4,2)-p_(1,3,3,2),
	p_(2,1,4,4)-p_(2,1,3,3),
	p_(4,1,2,4)-p_(3,1,2,3),
	p_(4,1,4,2)-p_(3,1,3,2),	
	p_(2,4,1,4)-p_(2,3,1,3),
        p_(4,2,1,4)-p_(3,2,1,3),
	p_(4,4,1,2)-p_(3,3,1,2),
	p_(2,4,4,1)-p_(2,3,3,1),
	p_(4,2,4,1)-p_(3,2,3,1),
	p_(4,4,2,1)-p_(3,3,2,1),
	p_(1,3,3,3)-p_(1,4,4,4), -- (a.3)
        p_(3,1,3,3)-p_(4,1,4,4),
	p_(3,3,1,3)-p_(4,4,1,4),
	p_(3,3,3,1)-p_(4,4,4,1),
	p_(2,4,4,4)-p_(2,3,3,3), -- (b)
	p_(4,2,4,4)-p_(3,2,3,3), 
	p_(4,4,2,4)-p_(3,3,2,3),
	p_(4,4,4,2)-p_(3,3,3,2),
	p_(3,3,4,4)-p_(4,4,3,3), --(c)
	p_(3,4,3,4)-p_(4,3,4,3),
	p_(3,4,4,3)-p_(4,3,3,4)} 


--Sanity check: 29 symmetry equations
numcols mingens F==29

--Sanity check: 
-- the equations in F are symmetry equations of the model
-- note that we are not proving that they are the only symmetry equations

apply(flatten entries gens F,i->f12(i))
apply(flatten entries gens F,i->f13(i))
apply(flatten entries gens F,i->f14(i))

------------------------------------------------------------------------------
------------------------------------------------------------------------------
-- 3. Proposition 5.2: additional binomial equations for a specific tree topology
------------------------------------------------------------------------------
------------------------------------------------------------------------------

G=ideal{p_(3,4,3,4),
        p_(3,4,4,3),
	p_(2,4,2,4)-p_(2,3,2,3),
	p_(2,4,4,2)-p_(2,3,3,2),
	p_(4,2,4,2)-p_(3,2,3,2),
	p_(4,2,2,4)-p_(3,2,2,3)}

apply(flatten entries gens G,i->f12(i)) --0: Equations in G hold for 12|34
apply(flatten entries gens G,i->f13(i)) --not 0: Equations in G do not hold for 13|24
apply(flatten entries gens G,i->f14(i)) --not 0: Equations in G do not hold for 14|23


------------------------------------------------------------------------------
------------------------------------------------------------------------------
-- 4. Theorem 7.3: complete intersection
------------------------------------------------------------------------------
------------------------------------------------------------------------------

-- Step 1: reproduce the steps we did to obtain the equations

          -- (a) equations arising from tripods (adding 1's in tripod equations in 7.1)

x1=p_(1,2,2,2)*p_(1,4,4,1)-p_(1,2,2,1)*p_(1,4,4,2)
x2=p_(1,2,2,2)*p_(1,4,1,4)-p_(1,2,1,2)*p_(1,4,2,4)
x3=p_(1,1,4,4)*p_(1,2,2,2)-p_(1,1,2,2)*p_(1,2,4,4)
x4=p_(1,1,4,4)*p_(1,4,1,4)*p_(1,4,4,1)-p_(1,1,1,1)*p_(1,4,4,4)^2
x5=p_(1,2,4,4)*p_(1,4,2,4)*p_(1,4,4,2)-p_(1,2,2,2)*p_(1,4,4,4)^2

y1=p_(2,2,2,1)*p_(4,4,1,1)-p_(2,2,1,1)*p_(4,4,2,1)
y2=p_(2,2,2,1)*p_(4,1,4,1)-p_(2,1,2,1)*p_(4,2,4,1)
y3=p_(1,4,4,1)*p_(2,2,2,1)-p_(1,2,2,1)*p_(2,4,4,1)
y4=p_(1,4,4,1)*p_(4,1,4,1)*p_(4,4,1,1)-p_(1,1,1,1)*p_(4,4,4,1)^2
y5=p_(2,4,4,1)*p_(4,2,4,1)*p_(4,4,2,1)-p_(2,2,2,1)*p_(4,4,4,1)^2


          -- (b), (c), (d) equations arising from rank conditions on the flattening matrix	   

-- flattening matrix for tree topology 12|34 for TN93
flat=value get "Auxiliary_files/Flat1234_ptilde.txt"; 

-- adding (model and topology) symmetry constraints to the flattening matrix (i.e. work in space of mixtures for tree topology 12|34)
flatF84=sub(flat,apply(flatten entries gens (F+ideal{p_(2,3,2,3)-p_(2,4,2,4),p_(2,3,3,2)-p_(2,4,4,2),p_(3,2,2,3)-p_(4,2,2,4),p_(3,2,3,2)-p_(4,2,4,2)}),i->(terms i)_0=>-(terms i)_1));

--build rank blocks (same as for TN93, see Quartet.m2 or "A novel algebraic approach to time-reversible evolutionary models")
s={(1,1),(1,4),(4,1),(2,4),(4,2),(4,4),(2,2),(1,2),(2,1),(3,3),(1,3),(3,1),(2,3),(3,2),(3,4),(4,3)}
--Blocks of rk 1
B11=flatF84_{position(s,i->i==(1,4)),position(s,i->i==(4,1)),position(s,i->i==(2,4)),position(s,i->i==(4,2))}^{position(s,i->i==(1,4)),position(s,i->i==(4,1)),position(s,i->i==(2,4)),position(s,i->i==(4,2)),position(s,i->i==(4,4))}
B12=flatF84_{position(s,i->i==(1,3)),position(s,i->i==(3,1)),position(s,i->i==(3,2)),position(s,i->i==(2,3))}^{position(s,i->i==(1,3)),position(s,i->i==(3,1)),position(s,i->i==(2,3)),position(s,i->i==(3,2)),position(s,i->i==(3,3))}
B13=flatF84_{position(s,i->i==(1,2)),position(s,i->i==(2,1))}^{position(s,i->i==(1,2)),position(s,i->i==(2,1)),position(s,i->i==(2,2)),position(s,i->i==(3,3)),position(s,i->i==(4,4))}
--Block of rk 2
B2=flatF84_{position(s,i->i==(1,1)),position(s,i->i==(1,2)),position(s,i->i==(2,2))}^{position(s,i->i==(1,1)),position(s,i->i==(1,2)),position(s,i->i==(2,1)),position(s,i->i==(2,2)),position(s,i->i==(3,3)),position(s,i->i==(4,4))}
--Blocks of rk 3
B31=flatF84_{position(s,i->i==(1,1)),position(s,i->i==(1,2)),position(s,i->i==(1,4)),position(s,i->i==(4,4))}^{position(s,i->i==(1,1)),position(s,i->i==(1,2)),position(s,i->i==(1,4)),position(s,i->i==(4,1)),position(s,i->i==(2,4)),position(s,i->i==(4,2)),position(s,i->i==(2,1)),position(s,i->i==(2,2)),position(s,i->i==(3,3)),position(s,i->i==(4,4))}
B32=flatF84_{position(s,i->i==(1,1)),position(s,i->i==(1,2)),position(s,i->i==(1,3)),position(s,i->i==(3,3))}^{position(s,i->i==(1,1)),position(s,i->i==(1,2)),position(s,i->i==(1,3)),position(s,i->i==(3,1)),position(s,i->i==(2,3)),position(s,i->i==(3,2)),position(s,i->i==(2,1)),position(s,i->i==(2,2)),position(s,i->i==(3,3)),position(s,i->i==(4,4))}

--rk 0 conditions
CI0=ideal{p_(3,4,3,4),p_(4,3,4,3),p_(3,4,4,3),p_(4,3,3,4)}
--rk 1 conditions
L1=join(flatten for i to (numcols B11)-1 list (for j to (numrows B11)-1 list det B11_{0,i}^{0,j}),
flatten for i to (numcols B12)-1 list (for j to (numrows B12)-1 list det B12_{0,i}^{0,j}),
flatten for i to (numcols B13)-1 list (for j to (numrows B13)-1 list det B13_{0,i}^{0,j}));
CI1F84=sum apply(L1,i->trim ideal i);
netList CI1F84_*
betti CI1F84 --28 quadrics
betti trim CI1F84 --15
--rk 2 conditions
L2=flatten for j to (numrows B2)-1 list det B2_{0,1,2}^{0,1,j};
CI2F84=sum apply(L2,i->trim ideal i);
betti CI2F84
--rk 3 conditions
L3=join(flatten for j to (numrows B31)-1 list det B31^{0,1,2,j},flatten for j to (numrows B32)-1 list det B32^{0,1,2,j});
CI3F84=sum apply(L3,i->trim ideal i);
betti CI3F84 --14 quartics
betti trim CI3F84 --9 quartics
--all ranks
CIF84=CI0+CI1F84+CI2F84+CI3F84;
betti CIF84 -- 4 linear, 28 quadric, 4 cubic, 14 quartic
betti trim CIF84 --4 linear, 15 quadric, 4 cubic, 9 quartic

--factor out variables p1111,p1212,p1414
ciF84=sum apply(flatten entries gens CIF84,i->saturate(ideal{i},ideal{p_(1,1,1,1)*p_(1,2,1,2)*p_(1,4,1,4)}));
betti ciF84 --4 linear, 37 quadrics, 7 cubics, 2 quartics
betti trim ciF84 --4 linear, 20 quadrics, 6 cubics, 2 quartics

--NOTE: The command trim computes a minimal system of generators of the ideal. 
--Doing this step BEFORE factoring out the variables using saturation affects the result.
--To obtain the polynomials in Theorem 7.3 we need to do it in the order described in the code.

 -- Step 2: check that the equations cut a local complete intersection

X=ideal{x1,x2,x3,x4,x5,y1,y2,y3,y4,y5}+ideal select(flatten entries gens (trim ciF84),i->degree i>{1});
betti X --26 quadrics, 10 cubics, 2 quartics
netList X_*

--Rank of the Jacobian at the no-evolution point
Jac=jacobian X;
noEvol=matrix{toList apply(gens Rp,i->sub(f12(i),apply(gens R,i->i=>1)))}
JacNoEvol=sub(Jac,noEvol);
rank JacNoEvol --38


------------------------------------------------------------------------------
------------------------------------------------------------------------------
-- 5. Corollary 7.5: rank constraints in flattening matrices
------------------------------------------------------------------------------
------------------------------------------------------------------------------

--F defined in Proposition 4.6
--G defined in Proposition 5.2
--Corollary 6.3 proves that F+G cuts the space of phylogenetic mixtures for T=12|34
E12=trim(F+G);

JA=trim(ideal select(flatten entries gens (trim ciF84),i->degree i>{1})+E12);
betti JA --35 linear, 20 quadrics, 6 cubics, 2 quartics 
IA=trim(CIF84+F);
betti IA --31 linear, 15 quadrics, 4 cubics, 9 quartics

--Saturations are equal
--WARNING: saturations do not finish (not in a couple of minutes at least)
JAsat=time saturate(JA,ideal{p_(1,1,1,1)*p_(1,2,1,2)*p_(1,4,1,4)});
IAsat=time saturate(IA,ideal{p_(1,1,1,1)*p_(1,2,1,2)*p_(1,4,1,4)});
JAsat==IAsat

-- Check Sage computations for Corollary 7.5 (Sage folder in this repository). 

------------------------------------------------------------------------------
------------------------------------------------------------------------------
-- 6. Corollary 7.10: variety of TN93 intersected with symmetry equations
------------------------------------------------------------------------------
------------------------------------------------------------------------------

-- WARNING: there is a mistake in this proof that will be solved during the review process
-- The local complete intersection provided in Theorem 7.6 is not the intersection between
-- the local complete intersection for TN93 provided in "A novel algebraic approach for time-reversible evolutionary models"
-- and symmetry equations for F81
-- However, there exists A DIFFERENT local complete intersection that satisfies the statement

-- Complete intersection ideal for TN93 in ptilde coordinates
CI=value get "Auxiliary_files/QuartetCI_ptilde.txt";
betti CI --4 linear, 40 quadrics, 10 cubics, 14 quartics

-- Complete intersection ideal for TN93 intersected with symmetry equations for F81
CIF=trim(CI+F);
betti CIF -- 31 linear, 25 quadrics, 8 cubics, 9 quartics

-- Check that this ideal is equal to the local complete intersection from Theorem 7.6 up to saturation
CIFsat=trim sum apply(flatten entries gens CIF,i->saturate(ideal{i},ideal{p_(1,1,1,1)*p_(1,2,1,2)*p_(1,4,1,4)}));
betti CIFsat --35 linear, 30 quadrics, 8 cubics
--WARNING: saturations do not finish (not in a couple of minutes at least)
Xsat=time saturate(X+E12,ideal{p_(1,1,1,1)*p_(1,2,1,2)*p_(1,4,1,4)});
CIFsat==Xsat

-- Prove that CIF is also a local complete intersection
JacCIF=jacobian CIF;
rank JacCIF --73
noEvol=matrix{toList apply(gens Rp,i->sub(f12(i),apply(gens R,i->i=>1)))}
JacCIFNoEvol=sub(JacCIF,noEvol);
rank JacCIFNoEvol --73
