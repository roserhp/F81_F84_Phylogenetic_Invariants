#define the p-tilde variables and the three constants 
p1111, p1122, p1212, p1221, p1222, p2121, p2211, p2221, p2112, p2122, p2212, p3333, p2232, p2332, p2333 = var('p1111, p1122, p1212, p1221, p1222, p2121, p2211, p2221, p2112, p2122, p2212, p3333, p2232, p2332, p2333')\
c1,c2,c3=var('c1,c2,c3');

#write the 7 equations
f1 = p1122*p1212*p1221 - p1111*p1222;
f2 = p1221*p2121*p2211 - p1111*p2221;
f3 = c1*p1111*p3333 - c2*p1111*p2332 - c3*p1122*p2211;
f4 = p1212*p2121 - p2112*p1221;
f5 = p1212*p2122 - p2112*p1222;
f6 = p1212*p2221 - p2212*p1221;
f7 = p1212*p2333 - p2212*p1222;

#compute the Jacobian for a generic point
J(p1111, p1212, p1221, p1222, p2121, p2211, p2221, p2112, p2122, p3333, p2232, p2333) = jacobian((f1,f2,f3,f4,f5,f6,f7), (p1111, p1212, p1221, p1222, p2121, p2211, p2221, p2112, p2122, p3333, p2232, p2333))\

#define the jacobian at the no-evolution point\
nep_jac = J(1,1,1,1,1,1,1,1,1,1,1,1)\

#compute the rank of the Jacobian at the no evolution point\
nep_jac.rank()}
