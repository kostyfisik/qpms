$Assumptions = k >= 0 && k0 >= 0 && c >= 0 && n >= 0 ;
Simplify[Refine[Integrate[(1 - Exp[-c x])^\[Kappa] (k0 x)^(-q) Exp[
    I k0 x] x BesselJ[n, k x], {x, 
       0, \[Infinity]}], {\[Kappa] == kk, q == qq, n == nn}]]
Series[%, {k, \[Infinity], 10}]
Quit[ ]
