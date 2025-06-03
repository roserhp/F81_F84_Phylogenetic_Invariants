-----------------------------------------------------------------------------------------------------------
-----------------------------------------------------------------------------------------------------------
--
-- Tripod.m2: this file contains the following computations regarding tripods:
--
-- 1. No-evolution point rho=varphi(Id,Id,Id) in probability coordinates
-- 2. No-evolution point rho=varphi(Id,Id,Id) in ptilde coordinates
-- 3. General point p=varphi(D_1,D_2,D_3) in ptilde coordinates for TN93
-- 4. Vanishing ideal in ptilde coordinates for TN93
-- 5. Complete intersection for TN93
--
-- For computational proofs of Propositions 3.2, 3.3, 3.4, 7.1, 7.4 see Tripod_propositions.m2
----------------------------------------------------------------------------------------------------------- 
----------------------------------------------------------------------------------------------------------- 
restart

-- Ring declaration: 
-- p_1,p_2,p_3,p_4 are parameters representing the root distribution p=(p_1,p_2,p_3,p_4)
-- variable l_(i,j) corresponds to eigenvalue j of transition matrix i

K=frac(QQ[p_1,p_2,p_3,p_4]);
R=K[l_(1,1)..l_(3,4)];

--------------------------------------------------------------------------
-- 1. No-evolution point rho=varphi(Id,Id,Id) in probability coordinates
--------------------------------------------------------------------------

-- Transition matrices of the no-evolution point
M1=id_(R^4)
M2=M1
M3=M1

-- No-evolution point in probability coordinates 
st={1,2,3,4};
S=sort elements (set st)^**3/splice;

rho=mutableMatrix(R,64,1)
for i to 63 do (
	rho_(i,0)=sum apply({1,2,3,4},k->p_k*M1_(position(st,l->l==k),position(st,l->l==(S_i)_0))
	                           *M2_(position(st,l->l==k),position(st,l->l==(S_i)_1))
		                   *M3_(position(st,l->l==k),position(st,l->l==(S_i)_2)))); 
rho=matrix rho; 

--How to check value of specific coordinates?
--rho_(j_1,j_2,j_3)=rho_(position(S,i->i==(j_1,j_2,j_3)),0)
rho_(position(S,i->i==(2,2,2)),0)

---------------------------------------------------------------------------
-- 2. No-evolution point rho=varphi(Id,Id,Id) in ptilde coordinates
---------------------------------------------------------------------------

--Matrix of change of basis
A=matrix{{p_1,p_1*p_3+p_1*p_4,0,(p_1*p_2)/(p_1+p_2)},{p_2,p_2*p_3+p_2*p_4,0,(-p_1*p_2)/(p_1+p_2)},{p_3,-p_1*p_3-p_2*p_3,(p_3*p_4)/(p_3+p_4),0},{p_4,-p_1*p_4-p_2*p_4,(-p_3*p_4)/(p_3+p_4),0}}
Ainv=inverse A

-- No-evolution point in pi-orthogonal coordinates 
rhobar=(Ainv**Ainv**Ainv)*rho;

-------------------------------------------------------------------------------
-- 3. General point p=varphi(D_1,D_2,D_3) in ptilde coordinates for TN93
-------------------------------------------------------------------------------

--Transition matrices of a general point 
D1=diagonalMatrix toList(l_(1,1)..l_(1,4))
D2=diagonalMatrix toList(l_(2,1)..l_(2,4))
D3=diagonalMatrix toList(l_(3,1)..l_(3,4))

-- General point in pi-orthogonal coordinates;
pbar=(D1**D2**D3)*rhobar;

-- General point in ptilde coordinates
ptilde=transpose matrix{toList apply(0..63,i->if(rhobar_(i,0)!=0) then (1/sub(rhobar_(i,0),K))*pbar_(i,0) else pbar_(i,0))};

-- View ptilde coordinates for a general point   
netList apply(S,i->{i,ptilde_(position(S,j->j==i),0)})

-------------------------------------------------------------------
-- 4. Vanishing ideal in ptilde coordinates for TN93
-------------------------------------------------------------------

-- Ring declaration: 
-- variables p_ijk (19 nonvanishing coordinates)

nonZeroEntries=S_(positions(flatten entries ptilde,i->i!=0))
"NonZeroEntries.txt" << toString nonZeroEntries << endl << close

varp=toList apply(nonZeroEntries,i->(symbol p)_i);
Rp=K[varp]; 
gens Rp

--Vanishing ideal as kernel of the map

PTILDE=flatten entries ptilde^(positions(flatten entries ptilde,i->i!=0));
"NonZeroPtilde.txt" << toString PTILDE << endl << close

f=map(R,Rp,PTILDE);
f(p_(1,1,1))
f(p_(4,4,4))
f(p_(3,3,3))
I=time trim kernel f;
-- used 0.531251 seconds
betti I
netList I_*

-------------------------------------------------------------------
-- 5. Complete intersection for TN93
-------------------------------------------------------------------

-- From the paper "A novel algebraic approach to time-reversible evolutionary models" 
-- Proposition 5.3, with ptilde coordinates

CI=ideal{p_(2,2,2)*p_(1,4,4)-p_(1,2,2)*p_(2,4,4), p_(2,2,2)*p_(4,1,4)-p_(2,1,2)*p_(4,2,4),
      p_(2,2,2)*p_(4,4,1)-p_(2,2,1)*p_(4,4,2), p_(1,4,4)*p_(2,3,3)-p_(1,3,3)*p_(2,4,4), p_(4,1,4)*p_(3,2,3)-p_(3,1,3)*p_(4,2,4), p_(4,4,1)*p_(3,3,2)-p_(3,3,1)*p_(4,4,2),
      -p_(1,1,1)*p_(3,3,3)^2+p_(1,3,3)*p_(3,1,3)*p_(3,3,1), -p_(1,1,1)*p_(4,4,4)^2+p_(1,4,4)*p_(4,1,4)*p_(4,4,1), -p_(2,2,2)*p_(3,3,3)^2+p_(2,3,3)*p_(3,2,3)*p_(3,3,2)}
betti CI
netList CI_*


-- Ring declaration:
-- Primary decomposition does not work over a field of fractions
-- Since pi does not appear in the generators of the ideal, we can work over QQ 
T=QQ[varp]
I=sub(I,T)
CI=sub(CI,T)

"TripodVanishingIdeal_ptilde.txt" << toString I << endl << close
"TripodCI_ptilde.txt" << toString CI << endl << close

-- V(I) is an irreducible component of V(CI),
-- saturation of V(CI) w.r.t. p_111*p_222*p_333*p_444 coincides with I
PD=primaryDecomposition CI;
length PD --132
PD_131==I --true
time saturate(CI,sub(ideal{p_(1,1,1)*p_(2,2,2)*p_(3,3,3)*p_(4,4,4)},T))==I --true

