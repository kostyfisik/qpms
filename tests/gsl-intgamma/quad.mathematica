SigmaIntegrand[n_, k_, R_, x_] := Exp[-(R x)^2 + (k/(2 x))^2] x^(2 n) ;
l:= ReadList["values.in", Table[Number, {4}]];
For[i = 1, i <= Length[l], i++,
  p := l[[i]];
  n := p[[1]];
  k := p[[2]];
  r := p[[3]];
  eta := p[[4]];
  res := NIntegrate[SigmaIntegrand[n, k, r, x], {x,eta,Infinity}];
  Print[InputForm[n], " ", CForm[k], " ", CForm[r], " ", CForm[res]]
];

