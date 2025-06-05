-----------------------------------------------------------------------------------------------------------
-----------------------------------------------------------------------------------------------------------
--
-- Quartet1234.m2: this file contains the following computations regarding quartets with tree topology 12|34:
--
-- 1. No-evolution point rho=psi(Id,Id,Id,Id,Id) and point q=psi(Id,Id,Id,Id,M) for TN93 in probability coordinates
-- 2. Point q=varphi(Id,Id,Id,Id,D) for TN93 in ptilde coordinates
-- 3. General point p=varphi(D_1,D_2,D_3,D_4,D) in ptilde coordinates for TN93
-- 4. Flattening matrix 12|34 for a general point p=varphi(D_1,D_2,D_3,D_4,D) in ptilde coordinates for TN93
-- 5. Complete intersection in ptilde coordinates for TN93
--
-- For computational proofs see Quartet_propositions.m2
----------------------------------------------------------------------------------------------------------- 
----------------------------------------------------------------------------------------------------------- 
restart

-- Ring declaration: 
-- p_1,p_2,p_3,p_4 are parameters representing the root distribution p=(p_1,p_2,p_3,p_4)
-- variable l_(i,j) corresponds to eigenvalue j of transition matrix i

K=frac(QQ[p_1,p_2,p_3,p_4]);
R=K[l_(1,1)..l_(5,4)];

-----------------------------------------------------------------------------------------------------------------
-- 1. No-evolution point rho=psi(Id,Id,Id,Id,Id) and point q=psi(Id,Id,Id,Id,M) for TN93 in probability coordinates
-----------------------------------------------------------------------------------------------------------------

-- Transition matrices at pendant edges
M1=id_(R^4)
M2=M1
M3=M1
M4=M1

-- Matrix of change of basis
A=sub(matrix{{p_1,p_1*p_3+p_1*p_4,0,(p_1*p_2)/(p_1+p_2)},{p_2,p_2*p_3+p_2*p_4,0,(-p_1*p_2)/(p_1+p_2)},{p_3,-p_1*p_3-p_2*p_3,(p_3*p_4)/(p_3+p_4),0},{p_4,-p_1*p_4-p_2*p_4,(-p_3*p_4)/(p_3+p_4),0}},R);
Ainv=inverse A

-- Transition matrix at the interior edge
M5=transpose(Ainv)*diagonalMatrix{l_(5,1),l_(5,2),l_(5,3),l_(5,4)}*transpose(A)

-- Check matrix M5 belongs to the model TN93
--parameter b
(1/p_2)*M5_(2,1)==(1/p_1)*M5_(2,0) --true
(1/p_2)*M5_(3,1)==(1/p_1)*M5_(2,0) --true
(1/p_1)*M5_(3,0)==(1/p_1)*M5_(2,0) --true
(1/p_3)*M5_(0,2)==(1/p_1)*M5_(2,0) --true
(1/p_3)*M5_(1,2)==(1/p_1)*M5_(2,0) --true
(1/p_4)*M5_(0,3)==(1/p_1)*M5_(2,0) --true
(1/p_4)*M5_(1,3)==(1/p_1)*M5_(2,0) --true
--parameter c
(1/p_2)*M5_(0,1)==(1/p_1)*M5_(1,0)
--parameter d
(1/p_3)*M5_(3,2)==(1/p_4)*M5_(2,3)


-- Point q=psi(Id,Id,Id,Id,M) in probability coordinates
st={1,2,3,4};
S=sort elements (set st)^**4/splice/splice;

q=mutableMatrix(R,256,1)
for i to 255 do (
	q_(i,0)=sum flatten toList apply(st,k->apply(st,kk->p_k*M1_(position(st,l->l==k),position(st,l->l==(S_i)_0))*M2_(position(st,l->l==k),position(st,l->l==(S_i)_1))*M3_(position(st,l->l==kk),position(st,l->l==(S_i)_2))*M4_(position(st,l->l==kk),position(st,l->l==(S_i)_3))*M5_(position(st,l->l==k),position(st,l->l==kk))))
) 
q=matrix q;

--How to check value of specific coordinates?
--q_(j_1,j_2,j_3,j_4)=q_(position(S,i->i==(j_1,j_2,j_3,j_4)),0)
q_(position(S,i->i==(2,2,2,2)),0)

-- No-evolution point rho=psi(Id,Id,Id,Id,Id) in probability coordinates
rho=transpose matrix{toList apply(flatten entries q,i->sub(i,toList apply(1..4,i->l_(5,i)=>1)))};
rho_(position(S,i->i==(2,2,2,2)),0)

---------------------------------------------------------------------------
-- 2. Point q=varphi(Id,Id,Id,Id,D) for TN93 in ptilde coordinates
---------------------------------------------------------------------------

-- No-evolution point in pi-orthogonal coordinates 
rhobar=(Ainv**Ainv**Ainv**Ainv)*rho;

-- Point q for TN93 in pi-orthogonal coordinates 
qbar=(Ainv**Ainv**Ainv**Ainv)*q;

-- Point q for TN93 in ptilde coordinates 
qtilde=transpose matrix{toList apply(0..255,i->if(rhobar_(i,0)!=0) then (1/sub(rhobar_(i,0),K))*qbar_(i,0) else qbar_(i,0))};

qtilde_(position(S,i->i==(2,2,2,2)),0)
qtilde_(position(S,i->i==(1,1,2,2)),0)
qtilde_(position(S,i->i==(2,4,2,4)),0)
qtilde_(position(S,i->i==(3,3,3,3)),0)

-------------------------------------------------------------------------------
-- 3. General point p=varphi(D_1,D_2,D_3,D_4,D) in ptilde coordinates for TN93
-------------------------------------------------------------------------------

--Transition matrices of a general point 
D1=diagonalMatrix toList(l_(1,1)..l_(1,4))
D2=diagonalMatrix toList(l_(2,1)..l_(2,4))
D3=diagonalMatrix toList(l_(3,1)..l_(3,4))
D4=diagonalMatrix toList(l_(4,1)..l_(4,4))

-- General point in pi-orthogonal coordinates;
ptilde=(D1**D2**D3**D4)*qtilde;

-- View ptilde coordinates for a general point   
netList apply(S,i->{i,ptilde_(position(S,j->j==i),0)})

--Indices corresponding to non-zero ptilde coordinates for TN93 for quartet 12|34
nonZeroEntries=S_(positions(flatten entries ptilde,i->i!=0))
"1234NonZeroEntries.txt" << toString nonZeroEntries << endl << close

--Non-zero ptilde coordinates for TN93 for quartet 12|34
PTILDE=flatten entries ptilde^(positions(flatten entries ptilde,i->i!=0));
"1234NonZeroPtilde.txt" << toString PTILDE << endl << close

---------------------------------------------------------------------------------------------------------------
-- 4. Flattening matrix 12|34 for a general point p=varphi(D_1,D_2,D_3,D_4,D) in ptilde coordinates for TN93
---------------------------------------------------------------------------------------------------------------

Rp = K[apply(nonZeroEntries, ind -> p_ind)];

-- artifact to avoid case separation when defining flattening for 12|34
-- note that pbar3344 and pbar4433 are different from zero for a general point for 12|34 
-- but zero for the no-evolution point, so ptilde=pbar
aux=mutableMatrix rhobar;
aux_(position(S,i->i==(3,3,4,4)),0)=1 
aux_(position(S,i->i==(4,4,3,3)),0)=1
rhobar2=sub(matrix aux,Rp);

s={(1,1),(1,4),(4,1),(2,4),(4,2),(4,4),(2,2),(1,2),(2,1),(3,3),(1,3),(3,1),(2,3),(3,2),(3,4),(4,3)}
flatt=mutableMatrix(Rp,16,16);
for i to 15 do for j to 15 do (if member(s_i|s_j,nonZeroEntries)==true then flatt_(i,j)=rhobar2_(position(S,k->k==s_i|s_j),0)*p_(s_i|s_j))
flatt=matrix flatt;    
"Flat1234_ptilde.txt" << toString flatt << endl << close

-- RANK BLOCKS 
-- From the paper "A novel algebraic approach to time-reversible evolutionary models" 
-- Lemma 5.13

--Blocks of rk 1
B11=flatt_{position(s,i->i==(1,4)),position(s,i->i==(4,1)),position(s,i->i==(2,4)),position(s,i->i==(4,2))}
B12=flatt_{position(s,i->i==(1,3)),position(s,i->i==(3,1)),position(s,i->i==(3,2)),position(s,i->i==(2,3))}
B13=flatt_{position(s,i->i==(1,2)),position(s,i->i==(2,1))}

--Block of rk 2
B2=flatt_{position(s,i->i==(1,1)),position(s,i->i==(1,2)),position(s,i->i==(2,2))}

--Blocks of rk 3
B31=flatt_{position(s,i->i==(1,1)),position(s,i->i==(1,2)),position(s,i->i==(1,4)),position(s,i->i==(4,4))}
B32=flatt_{position(s,i->i==(1,1)),position(s,i->i==(1,2)),position(s,i->i==(1,3)),position(s,i->i==(3,3))}

---------------------------------------------------------------------------------------
-- 5. Complete intersection for TN93
---------------------------------------------------------------------------------------

-- Chosing minors with a fixed non-vanishing subminor:

--rk 1 conditions
CI1=trim ((ideal flatten for i to (numcols B11)-1 list (for j to (numrows B11)-1 list det B11_{0,i}^{1,j}))+
(ideal flatten for i to (numcols B12)-1 list (for j to (numrows B12)-1 list det B12_{0,i}^{10,j}))+
(ideal flatten for i to (numcols B13)-1 list (for j to (numrows B13)-1 list det B13_{0,i}^{7,j})));
betti CI1 --28 quadrics: 12+12+4

--rk 2 conditions
CI2=trim (ideal flatten for j to (numrows B2)-1 list det B2_{0,1,2}^{0,7,j});
betti CI2 --4 cubics

--rk 3 conditions
CI3=trim (ideal flatten for j to (numrows B31)-1 list det B31^{0,7,1,j}+ideal flatten for j to (numrows B32)-1 list det B32^{0,7,10,j});
betti CI3 --14 quartics

--All ranks
CI=trim(CI1+CI2+CI3);
betti CI --28 quadrics, 4 cubics, 14 quartics

"Quartet1234CI_ptilde.txt" << toString CI << endl << close
