(* Compton and double Compton scattering in FeynCalc, for the results of J. Chluba's thesis
   "Spectral distortions of the CMB" (2005).  Run:  wolframscript -file compton_dc_feyncalc.wl
   (or the Wolfram kernel with -script) from the repository root.                              *)
$LoadAddOns = {};
Get["FeynCalc`"];
outfile = FileNameJoin[{Directory[], "examples", "chluba", "out", "chluba_feyncalc.txt"}];
str = OpenWrite[outfile];
say[x_] := (WriteString[str, x <> "\n"]; Print[x]);

(* ---------------------------------------------------------------- Compton *)
ClearScalarProducts[];
SP[p, p] = m^2; SP[k, k] = 0; SP[kp, kp] = 0;
SP[p, k] = (s - m^2)/2; SP[p, kp] = (m^2 - u)/2; SP[k, kp] = (s + u - 2 m^2)/2;
chain = GA[nu].(GS[p + k] + m).GA[mu]/(s - m^2) + GA[mu].(GS[p - kp] + m).GA[nu]/(u - m^2);
chainbar = GA[mu].(GS[p + k] + m).GA[nu]/(s - m^2) + GA[nu].(GS[p - kp] + m).GA[mu]/(u - m^2);
{tC, msqC} = AbsoluteTiming[Simplify[ExpandScalarProduct[Contract[TR[(GS[p + k - kp] + m).chain.(GS[p] + m).chainbar]]]]];
textbook = 8 (-(u - m^2)/(s - m^2) - (s - m^2)/(u - m^2) + 4 m^2 (1/(s - m^2) + 1/(u - m^2)) + 4 m^4 (1/(s - m^2) + 1/(u - m^2))^2);
say["Compton: sum |M|^2/e^4 - textbook = " <> ToString[Simplify[msqC - textbook]] <> "   (" <> ToString[Round[tC, 0.01]] <> " s)"];
(* Klein-Nishina in the lab frame: s - m^2 = 2 m w, u - m^2 = -2 m wp *)
wp = w/(1 + w (1 - c)/m);
avg = msqC/4 /. {s -> m^2 + 2 m w, u -> m^2 - 2 m wp};
dsig = (wp/w)^2 avg (4 Pi al)^2/(64 Pi^2 m^2);
kn = (al/m)^2/2 (wp/w)^2 (wp/w + w/wp - (1 - c^2));
say["Klein-Nishina: FeynCalc dsigma/dOmega - KN = " <> ToString[Simplify[dsig - kn]]];
say["Thomson limit: sigma = " <> ToString[InputForm[Integrate[2 Pi Limit[kn, w -> 0], {c, -1, 1}]]] <> "  (8 pi r0^2/3 with r0 = al/m)"];

(* ---------------------------------------------------------------- double Compton *)
mandlSkyrme[pk0_, pk1_, pk2_, ppk0_, ppk1_, ppk2_] := Module[{k0, k1, k2, k0p, k1p, k2p, a, b, cc, x, z, A, B, rho},
  {k0, k1, k2} = {-pk0, pk1, pk2}; {k0p, k1p, k2p} = {ppk0, -ppk1, -ppk2};
  a = 1/k0 + 1/k1 + 1/k2; b = 1/k0p + 1/k1p + 1/k2p; cc = 1/(k0 k0p) + 1/(k1 k1p) + 1/(k2 k2p);
  x = k0 + k1 + k2; z = k0 k0p + k1 k1p + k2 k2p; A = k0 k1 k2; B = k0p k1p k2p;
  rho = k0/k0p + k0p/k0 + k1/k1p + k1p/k1 + k2/k2p + k2p/k2;
  2 (a b - cc) ((a + b) (2 + x) - (a b - cc) - 8) - 2 x (a^2 + b^2) - 2 (a b + cc (1 - x)) rho - 8 cc +
   4 x/(A B) ((A + B) (1 + x) + x^2 (1 - z) + 2 z - (a A + b B) (2 + (1 - x) z/x))];

dcMsq[{pk0_, pk1_, pk2_, k0k1_, k0k2_}] := Module[{k1k2, pp, ph, chains, M, Mbar, den, res},
  ClearScalarProducts[];
  k1k2 = -(pk0 - pk1 - pk2 - k0k1 - k0k2);
  SP[p, p] = 1; SP[k0, k0] = 0; SP[k1, k1] = 0; SP[k2, k2] = 0;
  SP[p, k0] = pk0; SP[p, k1] = pk1; SP[p, k2] = pk2; SP[k0, k1] = k0k1; SP[k0, k2] = k0k2; SP[k1, k2] = k1k2;
  pp = p + k0 - k1 - k2;
  den[q_] := ExpandScalarProduct[SP[q, q]] - 1;
  ph = {{-k0, i1}, {k1, i2}, {k2, i3}};    (* photons as outgoing momenta, with their indices *)
  chains = Map[Function[o, {GA[o[[1, 2]]].(GS[pp + o[[1, 1]]] + 1).GA[o[[2, 2]]].(GS[p - o[[3, 1]]] + 1).GA[o[[3, 2]]],
       GA[o[[3, 2]]].(GS[p - o[[3, 1]]] + 1).GA[o[[2, 2]]].(GS[pp + o[[1, 1]]] + 1).GA[o[[1, 2]]],
       1/(den[pp + o[[1, 1]]] den[p - o[[3, 1]]])}], Permutations[ph]];
  M = Total[#[[1]] #[[3]] & /@ chains]; Mbar = Total[#[[2]] #[[3]] & /@ chains];
  res = -ExpandScalarProduct[Contract[TR[(GS[pp] + 1).M.(GS[p] + 1).Mbar]]];   (* (-g)^3 for the photons *)
  {res, mandlSkyrme[pk0, pk1, pk2, ExpandScalarProduct[SP[pp, k0]], ExpandScalarProduct[SP[pp, k1]], ExpandScalarProduct[SP[pp, k2]]]}];

Do[{t, {msq, X}} = AbsoluteTiming[dcMsq[pt]];
  say["double Compton at " <> ToString[pt] <> ": sum|M|^2/(e^6 X) = " <> ToString[Simplify[msq/X]] <> "   (" <> ToString[Round[t, 0.1]] <> " s)"],
  {pt, {{7/3, 26/25, 4/7, 35/8, 24/23}, {4/21, 7/3, 1/3, 27/7, 16/7}, {2, 4/23, 2/3, 19/3, 37/23}, {13/3, 5/2, 4, 19/18, 5/11}}}];
Close[str];
