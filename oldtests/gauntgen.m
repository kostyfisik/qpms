gaunt[m_, n_, mu_, nu_, 
   p_] := (-1)^(m + mu) (2 p + 1) Sqrt[
    Factorial[n + m] Factorial[
      nu + mu] Factorial[p - m - mu]/Factorial[n - m]/
      Factorial[nu - mu] / Factorial[p + m + mu]] ThreeJSymbol[{n, 
     0}, {nu, 0}, {p, 0}] ThreeJSymbol[{n, m}, {nu, 
     mu}, {p, -m - mu}]

lMax := 30
For[n = 0, n <= lMax, n++,
 For[nu = 0, nu <= lMax, nu++,
  For[m = -n, m <= n, m++,
   For[mu = -nu, mu <= nu, mu++,
    For[q = 0, q <= Min[n, nu, (n + nu - Abs[m + mu])/2],q++,
     Print[StringForm["{`1`, `2`, `3`, `4`, `5`, `6`},",m,n,mu,nu,q,CForm[N[gaunt[m, n, mu, nu, n + nu - 2 q],32]]]]
     ]
    ]
   ]
  ]
 ]

