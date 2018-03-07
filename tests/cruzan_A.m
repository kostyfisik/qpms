(* Vector translation coefficients as in Journal of Computational Physics 139, 137--165 (1998), eqs. (58), (59), (61) *)
gaunt[m_, n_, mu_, nu_, 
   p_] := (-1)^(m + mu) (2 p + 1) Sqrt[
    Factorial[n + m] Factorial[
      nu + mu] Factorial[p - m - mu]/Factorial[n - m]/
      Factorial[nu - mu] / Factorial[p + m + mu]] ThreeJSymbol[{n, 
     0}, {nu, 0}, {p, 0}] ThreeJSymbol[{n, m}, {nu, 
     mu}, {p, -m - mu}];
bCXcoeff[m_,n_,mu_,nu_,p_]:=(-1)^(mu+m)(2p+3)Sqrt[(n+m)!(nu+mu)!(p-m-mu+1)!/(n-m)!/(nu-mu)!/(p+m+mu+1)!]ThreeJSymbol[{n,m},{nu,mu},{p+1,-m-mu}]ThreeJSymbol[{n,0},{nu,0},{p,0}]
p[q_,n_,nu_]:=n+nu-2q;
ACXcoeff[m_,n_,mu_,nu_,q_]:=(-1)^m (2nu+1)(nu+m)!(nu-mu)!/2/n/(nu+1)/(nu-m)!/(nu+m)!I^p[q,n,nu](n(n+1)+nu(nu+1)-p[q,n,nu](p[q,n,nu]+1))gaunt[-m,n,mu,nu,p[q,n,nu]]
BCXcoeff[m_,n_,mu_,nu_,q_]:=(-1)^(m+1)(2nu+1)(n+m)!(nu-m)!/2/n/(n+1)/(n-m)!(nu+mu)!I^(p[q,n,nu]+1)Sqrt[((p[q,n,nu]+1)^2-(n-nu)^2)((n+nu+1)^2-(p[q,n,nu]+1)^2)]bCXcoeff[-m,n,mu,nu,p[q,n,nu]]`

lMax := 5
For[n = 0, n <= lMax, n++,
 For[nu = 0, nu <= lMax, nu++,
  For[m = -n, m <= n, m++,
   For[mu = -nu, mu <= nu, mu++,
    For[q = 0, q <= Min[n, nu, (n + nu - Abs[m + mu])/2],q++,
     Print[CForm[N[ACXcoeff[m,n,mu,nu,q],16]]]
     ]
    ]
   ]
  ]
 ]

