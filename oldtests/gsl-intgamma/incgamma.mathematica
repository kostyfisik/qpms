l:= ReadList["incgamma.in", Table[Number, {2}]];
For[i = 1, i <= Length[l], i++,
  p := l[[i]];
  j := p[[1]];
  x := p[[2]];
  res := NIntegrate[SigmaIntegrand[n, k, r, x], {x,eta,Infinity}];
  Print[j, " ", CForm[x], " | ", CForm[N[Gamma[1/2-j,x],16]], 
                          " | ", CForm[N[Gamma[1/2-j,x*I],16]], 
                          " | ", CForm[N[Gamma[1/2-j,-x],16]],
                          " | ", CForm[N[Gamma[1/2-j,-x*I],16]]]
];

