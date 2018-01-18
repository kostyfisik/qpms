$Assumptions = k >= 0 && k < k0 && k0 >= 0 && c >= 0 && n >= 0 ;
f = Refine[Integrate[(1 - Exp[-c x])^\[Kappa] (k0 x)^(-q) Exp[
    I k0 x] x BesselJ[n, k x], {x,
       0, \[Infinity]}], {\[Kappa] == kk, q == qq, n == nn}]
CForm[f]
Series[f, {k, \[Infinity], 10}]
Simplify[f]
Quit[ ]
