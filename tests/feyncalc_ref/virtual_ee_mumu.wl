(* FeynCalc + Package-X side of tests/test_feyncalc_loop.sage for e+e- -> mu+mu- in QED, masses kept.
   Run:  WHICH=vp|vertex|box WolframKernel -script virtual_ee_mumu.wl  (vacuum polarization, vertices, boxes). *)
$FeynCalcStartupMessages = False; $LoadAddOns = {"FeynArts", "FeynHelpers"};
Quiet[Get["FeynCalc`"]]; $FCAdvice = False;
t0 = AbsoluteTime[];
proc = {F[2, {1}], -F[2, {1}]} -> {F[2, {2}], -F[2, {2}]};
d0 = InsertFields[CreateTopologies[0, 2 -> 2], proc, InsertionLevel -> {Particles}, Restrictions -> QEDOnly, ExcludeParticles -> {F[3], F[4]}];
d1 = InsertFields[CreateTopologies[1, 2 -> 2, ExcludeTopologies -> {Tadpoles, WFCorrections}], proc,
   InsertionLevel -> {Particles}, Restrictions -> QEDOnly, ExcludeParticles -> {F[3], F[4]}];

opts = {IncomingMomenta -> {p1, p2}, OutgoingMomenta -> {p3, p4}, UndoChiralSplittings -> True, ChangeDimension -> D,
   List -> False, SMP -> True, Contract -> True, DropSumOver -> True};
amp0 = FCFAConvert[CreateFeynAmp[d0], Sequence @@ opts];
amp1all = FCFAConvert[CreateFeynAmp[d1, PreFactor -> 1], LoopMomenta -> {l}, Sequence @@ (opts /. (List -> False) -> (List -> True))];
nl[a_] := Count[a, PropagatorDenominator[m_, ___] /; !FreeQ[m, l], Infinity];
Print["loop propagators per diagram: ", nl /@ amp1all];
which = Environment["WHICH"];
amp1 = Total[Select[amp1all, Switch[which, "box", nl[#] == 4, "vp", nl[#] == 2, "vertex", nl[#] == 3, _, True] &]];
FCClearScalarProducts[];
SPD[p1, p1] = SMP["m_e"]^2; SPD[p2, p2] = SMP["m_e"]^2; SPD[p3, p3] = SMP["m_mu"]^2;
SPD[p1, p2] = (s - 2 SMP["m_e"]^2)/2; SPD[p1, p3] = (SMP["m_e"]^2 + SMP["m_mu"]^2 - t)/2;
(* p2.p3 from (p1 + p2 - p3)^2 = m_mu^2:  p2.p3 = m_e^2 + p1.p2 - p1.p3 *)
SPD[p2, p3] = SMP["m_e"]^2 + (s - 2 SMP["m_e"]^2)/2 - (SMP["m_e"]^2 + SMP["m_mu"]^2 - t)/2;
SPD[p4, p4] = SMP["m_mu"]^2;
x = (amp1 ComplexConjugate[amp0]) // FermionSpinSum // DiracSimplify;
x = FCReplaceMomenta[x, {p4 -> p1 + p2 - p3}] // ExpandScalarProduct;
Print["p4 left: ", !FreeQ[x, p4]];
Print["spin sum done ", Round[AbsoluteTime[] - t0]];
x = TID[x, l, ToPaVe -> True];
Print["TID done ", Round[AbsoluteTime[] - t0]];
x = If[which === "vertex", PaXEvaluate[x, PaXImplicitPrefactor -> 1], PaXEvaluate[x, PaXImplicitPrefactor -> 1, PaXC0Expand -> True, PaXD0Expand -> True]];
x = FeynAmpDenominatorExplicit[x] // ExpandScalarProduct;
Print["PaX done ", Round[AbsoluteTime[] - t0]];
pts = {{30, -9}, {50, -20}};
Do[
  v = x /. {s -> pt[[1]], t -> pt[[2]], SMP["m_e"] -> 1/2, SMP["m_mu"] -> 1, SMP["m_tau"] -> 3/2, SMP["e"] -> 1, ScaleMu -> 1, D -> 4 - 2 Epsilon};
  ser = Expand[Normal[Series[v, {Epsilon, 0, 0}]]];
  Print["POINT ", pt, " :: ", ToString[N[{Coefficient[ser, Epsilon, -2], Coefficient[ser, Epsilon, -1], Coefficient[ser, Epsilon, 0]}, 30], InputForm]],
  {pt, pts}];
Print["A0 :: ", ToString[PaXEvaluate[A0[m^2], PaXImplicitPrefactor -> 1] /. ScaleMu -> 1, InputForm]];
Print["DONE ", Round[AbsoluteTime[] - t0]];
