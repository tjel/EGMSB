(* ::Package:: *)

(* ::Input:: *)
fh2[x_]:=1/(72 x^2 (-1+x^2)^2) (24 (-1+x) x^2 (1+x)-24 (-1+x)^2 x^3 (1+x)Zeta[2]-12 (-1+x) x (1+x) (-10+7 x^2) ArcTanh[x]+48 x^3 ArcTanh[x] Log[2]-3 (-1+x) (16+15 x-17 x^2-11 x^3+3 x^4+2 x^5) Log[1-x]^2+48 x^5 ArcTanh[x] Log[2 x]-12 x^4 Log[x] Log[1+x]-3 (1+x)^2 (-16+31 x-10 x^2-5 x^3+4 x^4) Log[1+x]^2+Log[1-x] (12 x^2 (3+x+x^3) Log[1+x]-12 x^4 Log[x (1+x)])-6 (-1+x) (1+x) (-20+19 x^2) Log[1-x^2]+12 x^2 (2+x^4) Log[x] Log[1-x^2]+12 x^2 (1+x)^2 (2-2 x+x^2) PolyLog[2,1-x]-12 (-1+x)^2 x^2 (2+2 x+x^2) PolyLog[2,1/(1+x)]-24 x^3 (1+x^2) PolyLog[2,-1+2/(1+x)])+1/(6 x^2 (-1 + x^2))(-(-2 + x^2) (x ArcTanh[x] + Log[1 - x^2]));
