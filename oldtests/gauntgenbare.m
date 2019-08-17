gaunt[m_, n_, mu_, nu_, 
   p_] := (-1)^(m + mu) (2 p + 1) Sqrt[
    Factorial[n + m] Factorial[
      nu + mu] Factorial[p - m - mu]/Factorial[n - m]/
      Factorial[nu - mu] / Factorial[p + m + mu]] ThreeJSymbol[{n, 
     0}, {nu, 0}, {p, 0}] ThreeJSymbol[{n, m}, {nu, 
     mu}, {p, -m - mu}]

lMax := 18
For[n = 0, n <= lMax, n++,
 For[m = -n, m <= n, m++,
  For[nu = 0, nu <= lMax, nu++,
   For[mu = -nu, mu <= nu, mu++,
    For[q = 0, q <= Min[n, nu, (n + nu - Abs[m + mu])/2],q++,
     Print[CForm[N[gaunt[m, n, mu, nu, n + nu - 2 q],16]]]
     ]
    ]
   ]
  ]
 ]

