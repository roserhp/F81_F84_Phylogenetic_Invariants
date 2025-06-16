-----------------------------------------------------------------------------------------------------------
-----------------------------------------------------------------------------------------------------------
--
-- F81Quartet_propositions.m2: this file contains computations regarding quartets evolving under F81:
--
-- 0. Setup
-- 1. Parametrization and vanishing ideal for F81 with each tree topology
--
-- Computational proofs (unless stated otherwise) of the following results:
--
-- 2. Proposition 4.7 (sanity check): the listed equations are indeed linear model equations
-- 3. Proposition 5.4: the listed linear binomial equations hold for a specific tree topology
-- 4. Proposition 5.5: the listed linear non-binomial equations hold for a specific tree topology
-- 5. Proposition 6.4: space of mixtures of a specific tree topology
-- 6. Proposition 6.5: the listed equations are model equations
-- 7. Theorem 6.7: space of mixtures
-- 8. Theorem 7.6: complete intersection
-- 9. Corollary 7.8: rank constraints in flattening matrices
-- 10. Corollary 7.10 (partial proof): complete intersection of TN93 intersected with symmetry equations 
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
-- 1. Parametrization and vanishing ideal for F81 with each tree topology 
-----------------------------------------------------------------------

-- In model F81, eigenvalues for 2, 3 and 4 are equal
F81=flatten toList apply(1..5,i->{l_(i,3)=>l_(i,4),l_(i,2)=>l_(i,4)})

-- Non-zero ptilde coordinates for F81 for tree topology 12|34
ptilde12F81=flatten entries sub(ptilde12^(apply(nonZeroEntries,i->position(S,j->j==i))),F81)
ptilde13F81=flatten entries sub(ptilde13^(apply(nonZeroEntries,i->position(S,j->j==i))),F81)
ptilde14F81=flatten entries sub(ptilde14^(apply(nonZeroEntries,i->position(S,j->j==i))),F81)

--Parametrization map
f12=map(R,Rp,ptilde12F81);
f13=map(R,Rp,ptilde13F81);
f14=map(R,Rp,ptilde14F81);

--Vanishing ideal
I12=time trim kernel f12;
-- used 2.53137 seconds
I13=time trim kernel f13;
 -- used 2.83943 seconds
I14=time trim kernel f14;
-- used 2.73344 seconds

betti I12 -- 71 linear, 10 quadrics, 23 cubics
betti I13 -- 71 linear, 10 quadrics, 23 cubics
betti I14 --- 71 linear, 10 quadrics, 23 cubics



------------------------------------------------------------------
-------------------------------------------------------------------
-- 2. Proposition 4.7 symmetry equations for F81
-------------------------------------------------------------------
-------------------------------------------------------------------

-- Symmetry equations for F81
F=trim ideal{p_(1,1,4,4)-p_(1,1,2,2), -- (a)
        p_(1,4,1,4)-p_(1,2,1,2),
	p_(1,4,4,1)-p_(1,2,2,1),
	p_(4,1,1,4)-p_(2,1,1,2),
	p_(4,1,4,1)-p_(2,1,2,1),
	p_(4,4,1,1)-p_(2,2,1,1),
	p_(1,2,4,4)-p_(1,4,4,4), 
        p_(1,4,2,4)-p_(1,4,4,4),
	p_(1,4,4,2)-p_(1,4,4,4),
	p_(4,1,2,4)-p_(4,1,4,4),
	p_(4,1,4,2)-p_(4,1,4,4),
	p_(4,4,1,2)-p_(4,4,1,4),
	p_(2,1,4,4)-p_(4,1,4,4),
        p_(2,4,1,4)-p_(4,4,1,4),
	p_(2,4,4,1)-p_(4,4,4,1),
	p_(4,2,1,4)-p_(4,4,1,4),
	p_(4,2,4,1)-p_(4,4,4,1),
	p_(4,4,2,1)-p_(4,4,4,1),
        p_(1,2,2,2)-p_(1,4,4,4),
        p_(2,1,2,2)-p_(4,1,4,4),
	p_(2,2,1,2)-p_(4,4,1,4),
	p_(2,2,2,1)-p_(4,4,4,1),
	p_(1,1,3,3)-p_(1,1,2,2), 
        p_(1,3,1,3)-p_(1,2,1,2),
	p_(1,3,3,1)-p_(1,2,2,1),
	p_(3,1,1,3)-p_(2,1,1,2),
	p_(3,1,3,1)-p_(2,1,2,1),
	p_(3,3,1,1)-p_(2,2,1,1),
	p_(1,2,3,3)-p_(1,4,4,4), 
        p_(1,3,2,3)-p_(1,4,4,4),
	p_(1,3,3,2)-p_(1,4,4,4),
	p_(3,1,2,3)-p_(4,1,4,4),
	p_(3,1,3,2)-p_(4,1,4,4),
	p_(3,3,1,2)-p_(4,4,1,4),
	p_(2,1,3,3)-p_(4,1,4,4),
        p_(2,3,1,3)-p_(4,4,1,4),
	p_(2,3,3,1)-p_(4,4,4,1),
	p_(3,2,1,3)-p_(4,4,1,4),
	p_(3,2,3,1)-p_(4,4,4,1),
	p_(3,3,2,1)-p_(4,4,4,1),
        p_(1,3,3,3)-p_(1,4,4,4),
        p_(3,1,3,3)-p_(4,1,4,4),
	p_(3,3,1,3)-p_(4,4,1,4),
	p_(3,3,3,1)-p_(4,4,4,1),
	p_(4,2,4,4)-p_(2,4,4,4), -- (b)
	p_(4,4,2,4)-p_(2,4,4,4),
	p_(4,4,4,2)-p_(2,4,4,4),
	p_(2,3,3,3)-p_(2,4,4,4), 
	p_(3,2,3,3)-p_(2,4,4,4), 
	p_(3,3,2,3)-p_(2,4,4,4),
	p_(3,3,3,2)-p_(2,4,4,4),
	p_(2,2,3,3)-p_(3,3,2,2), -- (c)
	p_(2,3,2,3)-p_(3,2,3,2),
	p_(2,3,3,2)-p_(3,2,2,3),
	p_(2,2,4,4)-p_(4,4,2,2),
	p_(2,4,2,4)-p_(4,2,4,2),
	p_(2,4,4,2)-p_(4,2,2,4),
	p_(3,3,4,4)-p_(4,4,3,3),
	p_(3,4,3,4)-p_(4,3,4,3),
	p_(3,4,4,3)-p_(4,3,3,4)}

--Sanity check: 60 symmetry equations
numcols mingens F==60

--Sanity check: 
-- the equations in F are symmetry equations of the model
-- note that we are not proving that they are the only symmetry equations

apply(flatten entries gens F,i->f12(i))
apply(flatten entries gens F,i->f13(i))
apply(flatten entries gens F,i->f14(i))

------------------------------------------------------------------------------
------------------------------------------------------------------------------
-- 3. Proposition 5.4: additional binomial equations for a specific tree topology
------------------------------------------------------------------------------
------------------------------------------------------------------------------

G=ideal{p_(3,4,3,4),
        p_(3,4,4,3),
	p_(2,4,2,4)-p_(2,4,4,4),
	p_(2,4,4,2)-p_(2,4,4,4),
	p_(2,3,2,3)-p_(2,4,4,4),
	p_(2,3,3,2)-p_(2,4,4,4)}

--Sanity check
isSubset(F+G,I12) --true
numgens trim(F+G) --66

-- there are 5 more linear equations (see Prop 5.5)

------------------------------------------------------------------------------
------------------------------------------------------------------------------
-- 4. Proposition 5.5: non-binomial equations for a specific tree topology
------------------------------------------------------------------------------
------------------------------------------------------------------------------

--Check the listed polynomials belong to the ideal.
--Note that in the equations of the paper it has been sometimes used to simplify expressions
--that p_1+p_2+p_3+p_4=1

a=(p_3+p_4)^3*p_(2,2,4,4)-((p_1+p_2)^3+(p_3+p_4)^3)*p_(2,2,2,2)+(p_1+p_2)^3*p_(2,2,3,3)
f12(a) --0

b=(p_1+p_2)^2*(p_3^3+p_4^3)*p_(3,3,3,3)-(p_1+p_2)^2*(p_3+p_4)*(p_3-p_4)^2*p_(2,4,4,4)-p_3*p_4*(p_1+p_2)^2*(p_3+p_4)*p_(2,2,3,3)+p_1*p_2*p_3^2*p_4^2*(p_3+p_4)^2*p_(3,3,4,4)
f12(b) --0

c=(p_3+p_4)^2*(p_1^3+p_2^3)*p_(4,4,4,4)-(p_1+p_2)*(p_3+p_4)^2*(p_1-p_2)^2*p_(2,4,4,4)-p_1*p_2*(p_1+p_2)*(p_3+p_4)^2*p_(2,2,4,4)+p_1^2*p_2^2*p_3*p_4*(p_1+p_2)^2*p_(3,3,4,4)
f12(c) --0

d=(p_3+p_4)^2*p_(2,2,4,4)-(p_1+p_2)^2*p_(2,2,3,3)-((p_3+p_4)^2-(p_1+p_2)^2)*p_(2,4,4,4)
f12(d) -- Not 0
sub(f12(d),p_4=>1-p_1-p_2-p_3) --0

e=p_1^2*p_2^2*p_3*p_4*(p_1+p_2)*p_(3,3,4,4)-(p_3+p_4)*(p_1^3+p_2^3)*(p_(4,4,4,4)-p_(2,4,4,4))
f12(e) -- 0

f=(p_3^3+p_4^3)*p_(1,1,1,1)*p_(3,3,3,3)-(p_3^3+p_4^3-p_3*p_4*(p_3+p_4)^2)*p_(1,1,1,1)*p_(2,4,4,4)-p_3*p_4*(p_3+p_4)^2*p_(1,1,4,4)*p_(4,4,1,1)
f12(f) -- 0
sub(f12(f),p_4=>1-p_1-p_2-p_3) --0


------------------------------------------------------------------------------
------------------------------------------------------------------------------
-- 5. Proposition 6.4: space of mixtures for a specific tree topology
------------------------------------------------------------------------------
------------------------------------------------------------------------------

--F defined in Proposition 4.7
--G defined in Proposition 5.4
--a..e defined in Proposition 5.5

E12=trim(F+G+ideal{a,b,c,d,e});
I12linear=ideal select(flatten entries gens I12,i->degree i=={1});
I12linear==E12 --true


------------------------------------------------------------------------------
------------------------------------------------------------------------------
-- 6. Proposition 6.5: non-binomial linear model equations
------------------------------------------------------------------------------
------------------------------------------------------------------------------

a2=(p_3+p_4)^2*p_(2,2,4,4)-(p_1+p_2)^2*p_(2,2,3,3)-((p_3+p_4)^2-(p_1+p_2)^2)*p_(2,4,4,4)
f12(a2)--0
f13(a2)--0
f14(a2)--0
b2=(p_3+p_4)^2*p_(2,4,2,4)-(p_1+p_2)^2*p_(2,3,2,3)-((p_3+p_4)^2-(p_1+p_2)^2)*p_(2,4,4,4)
f12(b2)--0
f13(b2)--0
f14(b2)--0
c2=(p_3+p_4)^2*p_(2,4,4,2)-(p_1+p_2)^2*p_(2,3,3,2)-((p_3+p_4)^2-(p_1+p_2)^2)*p_(2,4,4,4)
f12(c2)--0
f13(c2)--0
f14(c2)--0
d2=p_1*p_2*p_3*p_4*p_(3,3,4,4)-(p_1+p_2)^2*(p_(2,2,3,3)-p_(2,4,4,4))
f12(d2)--not 0
sub(f12(d2),p_4=>1-p_1-p_2-p_3) --0
f13(d2)--0
f14(d2)--0
e2=p_1*p_2*p_3*p_4*p_(3,4,3,4)-(p_1+p_2)^2*(p_(2,3,2,3)-p_(2,4,4,4))
f12(e2)--0
f13(e2)--not 0
sub(f12(e2),p_4=>1-p_1-p_2-p_3) --0
f14(e2)--0
f2=p_1*p_2*p_3*p_4*p_(3,4,4,3)-(p_1+p_2)^2*(p_(2,3,3,2)-p_(2,4,4,4))
f12(f2) --0
f14(f2) --not 0
sub(f14(f2),p_4=>1-p_1-p_2-p_3) --0
g2=p_1*p_2*p_3*p_4*(p_(3,3,4,4)+p_(3,4,3,4)+p_(3,4,4,3))-((p_1+p_2)^3+(p_3+p_4)^3)*(p_(2,2,2,2)-p_(2,4,4,4))
f12(g2) --not 0
sub(f12(g2),p_4=>1-p_1-p_2-p_3) --0
h2=(p_1+p_2)*(p_3^3+p_4^3)*(p_(3,3,3,3)-p_(2,4,4,4))-p_3*p_4*(p_3+p_4)*((p_1+p_2)^3+(p_3+p_4)^3)*(p_(2,2,2,2)-p_(2,4,4,4))
f12(h2) --not 0
sub(f12(h2),p_4=>1-p_1-p_2-p_3) --0
i2=(p_3+p_4)*(p_1^3+p_2^3)*(p_(4,4,4,4)-p_(2,4,4,4))-p_1*p_2*(p_1+p_2)*((p_1+p_2)^3+(p_3+p_4)^3)*(p_(2,2,2,2)-p_(2,4,4,4))
f12(i2) --not 0
sub(f12(i2),p_4=>1-p_1-p_2-p_3) --0

------------------------------------------------------------------------------
------------------------------------------------------------------------------
-- 7. Theorem 6.7: space of mixtures
------------------------------------------------------------------------------
------------------------------------------------------------------------------

-- Union of phylogenetic varieties (intersection of ideals)

-- Note that it takes a few minutes to run
-- Uncomment if you want to rerun
--I=time intersect(I12,I13,I14);
-- used 470.256 seconds
--betti I --69,3,91
--do not rerun unless in need of overwriting
--"I.txt" << toString I << endl << close
I=value get "Auxiliary_files/QuartetIntersectionVanishingIdeals_ptilde.txt";

-- Comparing described equations to linear part of I
E=trim(F+ideal{a2,b2,c2,d2,e2,f2,g2,h2,i2});
Ilinear=ideal select(flatten entries gens I,i->degree i=={1});
Ilinear==E --false
sub(Ilinear,p_4=>1-p_1-p_2-p_3)==sub(E,p_4=>1-p_1-p_2-p_3) --true


------------------------------------------------------------------------------
------------------------------------------------------------------------------
-- 8. Theorem 7.6: complete intersection
------------------------------------------------------------------------------
------------------------------------------------------------------------------

-- Step 1: reproduce the steps we did to obtain the equations

          -- 1.1 cubic equations in the statement have been obtained by adding 1's in tripod cubic equation

xa=p_(1,1,4,4)*p_(1,4,1,4)*p_(1,4,4,1)-p_(1,1,1,1)*p_(1,4,4,4)^2
xb=p_(1,4,4,1)*p_(4,1,4,1)*p_(4,4,1,1)-p_(1,1,1,1)*p_(4,4,4,1)^2

          -- 1.2 quadric equations are obtained from the rank blocks of the flattening matrix	   

-- flattening matrix for tree topology 12|34 for TN93
flat=value get "Auxiliary_files/Flat1234_ptilde.txt"; 
-- adding symmetry constraints to the flattening matrix
flatF81=sub(flat,apply(flatten entries gens F,i->(terms i)_0=>-(terms i)_1));

--build rank blocks (same as for TN93, see Quartet.m2 or "A novel algebraic approach to time-reversible evolutionary models")
s={(1,1),(1,4),(4,1),(2,4),(4,2),(4,4),(2,2),(1,2),(2,1),(3,3),(1,3),(3,1),(2,3),(3,2),(3,4),(4,3)}
--Blocks of rk 1
B11=flatF81_{position(s,i->i==(1,4)),position(s,i->i==(4,1)),position(s,i->i==(2,4)),position(s,i->i==(4,2))}
B12=flatF81_{position(s,i->i==(1,3)),position(s,i->i==(3,1)),position(s,i->i==(3,2)),position(s,i->i==(2,3))}
B13=flatF81_{position(s,i->i==(1,2)),position(s,i->i==(2,1))}
--Block of rk 2
B2=flatF81_{position(s,i->i==(1,1)),position(s,i->i==(1,2)),position(s,i->i==(2,2))}
--Blocks of rk 3
B31=flatF81_{position(s,i->i==(1,1)),position(s,i->i==(1,2)),position(s,i->i==(1,4)),position(s,i->i==(4,4))}
B32=flatF81_{position(s,i->i==(1,1)),position(s,i->i==(1,2)),position(s,i->i==(1,3)),position(s,i->i==(3,3))}

--rk 0 conditions
CI0=ideal{p_(3,4,3,4),p_(4,3,4,3),p_(3,4,4,3),p_(4,3,3,4)}
--rk 1 conditions
CI1F81=trim ((ideal flatten for i to (numcols B11)-1 list (for j to (numrows B11)-1 list det B11_{0,i}^{1,j}))+
(ideal flatten for i to (numcols B12)-1 list (for j to (numrows B12)-1 list det B12_{0,i}^{10,j}))+
(ideal flatten for i to (numcols B13)-1 list (for j to (numrows B13)-1 list det B13_{0,i}^{7,j})));
--rk 2 conditions
CI2F81=trim (ideal flatten for j to (numrows B2)-1 list det B2_{0,1,2}^{0,7,j});
--rk 3 conditions
CI3F81=trim (ideal flatten for j to (numrows B31)-1 list det B31^{0,7,1,j}+ideal flatten for j to (numrows B32)-1 list det B32^{0,7,10,j});
--All ranks
CIF81=trim(CI0+CI1F81+CI2F81+CI3F81);
betti CIF81 --4 linear, 8 quadrics, 3 cubics, 3 quartics
netList CIF81_*

--Remark 5.1: p_(1,1,1,1) and p_(1,4,1,4) can be factored out in each equation 
ciF81=trim sum apply(flatten entries gens CIF81,i->saturate(ideal{i},ideal{p_(1,1,1,1)*p_(1,4,1,4)}));
betti ciF81 --13 linear, 5 quadrics
netList ciF81_*
ciF81linear=ideal select(flatten entries gens ciF81,i->degree i=={1});
isSubset(ciF81linear,E12) --true

xc=ciF81_15
xd=ciF81_17
xe=ciF81_16
xf=ciF81_14
xg=ciF81_13

-- Step 2: check that the equations cut a local complete intersection

X=ideal{xa,xb,xc,xd,xe,xf,xg};
netList X_*
betti X

--Rank of the Jacobian at the no-evolution point
Jac=jacobian X;
noEvol=matrix{toList apply(gens Rp,i->sub(f12(i),apply(gens R,i->i=>1)))}
JacNoEvol=sub(Jac,noEvol)
rank JacNoEvol --7

-- Step 3: check that the saturation of the complete intersection is the whole vanishing ideal
I12==saturate(X+E12,ideal{p_(1,1,1,1)*p_(1,4,1,4)*p_(4,1,4,1)}) --true

------------------------------------------------------------------------------
------------------------------------------------------------------------------
-- 9. Corollary 7.8: rank constraints in flattening matrices
------------------------------------------------------------------------------
------------------------------------------------------------------------------

JA=ideal{xc,xd,xe,xf,xg}+E12;
IA=CIF81+F;

--Saturations are equal
saturate(JA,ideal{p_(1,1,1,1)*p_(1,4,1,4)})==saturate(IA,ideal{p_(1,1,1,1)*p_(1,4,1,4)}) --true

--Saturation contained in the vanishing ideal
JAsaturated=saturate(JA,ideal{p_(1,1,1,1)*p_(1,4,1,4)});
isSubset(JAsaturated,I12) --true

--Saturation equals vanishing ideal up to degree two
I12quadratic=ideal select(flatten entries gens I12,i->degree i=={1} or degree i=={2});
JAsaturated==I12quadratic --true

------------------------------------------------------------------------------
------------------------------------------------------------------------------
-- 10. Corollary 7.10: variety of TN93 intersected with symmetry equations
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
betti CIF -- 62 linear, 8 quadrics, 5 cubics, 3 quartics

-- Check that this ideal is equal to the local complete intersection from Theorem 7.6 up to saturation
CIFsat=trim sum apply(flatten entries gens CIF,i->saturate(ideal{i},ideal{p_(1,1,1,1)*p_(1,4,1,4)}));
betti CIFsat --71 linear, 5 quadrics, 2 cubics
CIFsat==X+E12 --true

-- Prove that CIF is also a local complete intersection
JacCIF=jacobian CIF;
rank JacCIF --78
JacCIFNoEvol=sub(JacCIF,noEvol);
rank JacCIFNoEvol --78
