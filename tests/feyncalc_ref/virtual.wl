(* FeynCalc + Package-X side of tests/test_feyncalc_loop.sage: sum over spins, polarizations and colours of
   M1 M0^* at numerical points for H -> b b~ (gluon loop), e+e- -> gamma gamma (QED) and Z -> nu nu~ (electroweak).
   Run:  PROC=hbb WolframKernel -script virtual.wl  (PROC = hbb, ee_aa or z_nn).
   e+e- -> gamma gamma is not used in the test: TID meets a zero Gram determinant and Package-X gives Indeterminate. *)
$FeynCalcStartupMessages = False; $LoadAddOns = {"FeynArts", "FeynHelpers"};
Quiet[Get["FeynCalc`"]]; $FCAdvice = False;
which = Environment["PROC"]; t0 = AbsoluteTime[];
If[which === "z_nn", FCSetDiracGammaScheme["NDR-Discard"]];
Switch[which,
 "hbb", proc = {S[1]} -> {F[4, {3}], -F[4, {3}]}; nout = 2; model = SMQCD; sel = !FreeQ[#, V[5, ___]] &; restr = {};
   vecs = {}; pts = {{}}; num = {SMP["m_H"] -> 125, SMP["m_b"] -> 5, SMP["m_W"] -> 80, SMP["sin_W"] -> 1/2, SMP["e"] -> 1/3, SMP["g_s"] -> 1},
 "ee_aa", proc = {F[2, {1}], -F[2, {1}]} -> {V[1], V[1]}; nout = 2; model = SM; sel = True &; restr = {Restrictions -> QEDOnly, ExcludeParticles -> {F[3], F[4]}};
   vecs = {p3, p4}; pts = {{30, -9}, {50, -20}}; num = {SMP["m_e"] -> 1/2, SMP["m_mu"] -> 1, SMP["m_tau"] -> 3/2, SMP["e"] -> 1},
 "z_nn", proc = {V[2]} -> {F[1, {1}], -F[1, {1}]}; nout = 2; model = SM; sel = True &; restr = {};
   vecs = {p1}; pts = {{}}; num = {SMP["m_Z"] -> 91, SMP["m_W"] -> 80, SMP["cos_W"] -> 80/91, SMP["sin_W"] -> Sqrt[1 - (80/91)^2], SMP["m_e"] -> 1/2, SMP["e"] -> 1/3, SMP["m_H"] -> 125}];
n2 = Length[proc[[1]]] == 2;
d0 = InsertFields[CreateTopologies[0, Length[proc[[1]]] -> nout], proc, InsertionLevel -> {Particles}, Model -> model, Sequence @@ restr];
d1 = InsertFields[CreateTopologies[1, Length[proc[[1]]] -> nout, ExcludeTopologies -> {Tadpoles, WFCorrections}], proc,
   InsertionLevel -> {Particles}, Model -> model, Sequence @@ restr];
d1 = DiagramSelect[d1, sel];
moms = If[n2, {{p1, p2}, {p3, p4}}, {{p1}, {p2, p3}}];
opts = {IncomingMomenta -> moms[[1]], OutgoingMomenta -> moms[[2]], UndoChiralSplittings -> True, ChangeDimension -> D,
   SMP -> True, Contract -> True, DropSumOver -> True};
amp0 = FCFAConvert[CreateFeynAmp[d0], List -> False, Sequence @@ opts];
amp1 = FCFAConvert[CreateFeynAmp[d1, PreFactor -> 1], LoopMomenta -> {l}, List -> False, Sequence @@ opts];
FCClearScalarProducts[];
If[n2,
  {m1, m2, m3, m4} = Switch[which, "ee_aa", {SMP["m_e"], SMP["m_e"], 0, 0}];
  SetMandelstam[s, t, u, p1, p2, -p3, -p4, m1, m2, m3, m4],
  {m1, m2, m3} = Switch[which, "hbb", {SMP["m_H"], SMP["m_b"], SMP["m_b"]}, "z_nn", {SMP["m_Z"], 0, 0}];
  SPD[p1, p1] = m1^2; SPD[p2, p2] = m2^2; SPD[p3, p3] = m3^2; SPD[p2, p3] = (m1^2 - m2^2 - m3^2)/2;
  SPD[p1, p2] = (m1^2 + m2^2 - m3^2)/2; SPD[p1, p3] = (m1^2 - m2^2 + m3^2)/2];
x = (amp1 ComplexConjugate[amp0]) // FermionSpinSum // DiracSimplify // SUNSimplify[#, Explicit -> True, SUNNToCACF -> False] &;
Do[x = If[which === "z_nn", DoPolarizationSums[x, v], DoPolarizationSums[x, v, 0]], {v, vecs}];
x = x // DiracSimplify // SUNSimplify[#, SUNNToCACF -> False] & // ExpandScalarProduct;
x = FCReplaceMomenta[x, If[n2, {p4 -> p1 + p2 - p3}, {p3 -> p1 - p2}]] // ExpandScalarProduct;
If[n2, x = x /. u -> 2 SMP["m_e"]^2 - s - t];
Print["algebra done ", Round[AbsoluteTime[] - t0]];
x = TID[x, l, ToPaVe -> True];
x = PaXEvaluate[x, PaXImplicitPrefactor -> 1, PaXC0Expand -> True, PaXD0Expand -> True];
x = FeynAmpDenominatorExplicit[x] // ExpandScalarProduct;
Print["PaX done ", Round[AbsoluteTime[] - t0]];
Do[
  v = x /. If[n2, {s -> pt[[1]], t -> pt[[2]]}, {}] /. num /. {SUNN -> 3, ScaleMu -> 1, D -> 4 - 2 Epsilon};
  ser = Expand[Normal[Series[v, {Epsilon, 0, 0}]]];
  Print["POINT ", pt, " :: ", ToString[N[{Coefficient[ser, Epsilon, -2], Coefficient[ser, Epsilon, -1], Coefficient[ser, Epsilon, 0]}, 30], InputForm]],
  {pt, pts}];
Print["DONE ", Round[AbsoluteTime[] - t0]];
