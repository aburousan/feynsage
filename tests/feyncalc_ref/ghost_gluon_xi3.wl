(* FeynCalc + Package-X side of the ghost-gluon vertex check at xi = 3 in tests/test_feyncalc_loop.sage.
   The projection on p2 keeps a PaVe function that Package-X does not evaluate; only the p3 projection is used. *)
$FeynCalcStartupMessages = False; $LoadAddOns = {"FeynArts", "FeynHelpers"};
Quiet[Get["FeynCalc`"]]; $FCAdvice = False;
d1 = InsertFields[CreateTopologies[1, 1 -> 2, ExcludeTopologies -> {Tadpoles, WFCorrections}], {U[5]} -> {U[5], V[5]},
   InsertionLevel -> {Particles}, Model -> SMQCD];
d1 = DiagramSelect[d1, FreeQ[#, F] && FreeQ[#, S] && FreeQ[#, V[1] | V[2] | V[3]] &];
Print["diagrams: ", Length[CreateFeynAmp[d1]]];
amp = FCFAConvert[CreateFeynAmp[d1, Truncated -> True, GaugeRules -> {}, PreFactor -> 1], IncomingMomenta -> {p1},
   OutgoingMomenta -> {p2, p3}, LorentzIndexNames -> {mu}, SUNIndexNames -> {a, b, c}, LoopMomenta -> {l},
   ChangeDimension -> D, List -> False, SMP -> True, DropSumOver -> True, FinalSubstitutions -> {p1 -> p2 + p3}];
FCClearScalarProducts[]; SPD[p2, p2] = 2; SPD[p3, p3] = 3; SPD[p2, p3] = -1/2;
Do[
  x = Contract[amp FVD[v, mu]] // SUNSimplify[# SUNF[a, b, c], Explicit -> True, SUNNToCACF -> False] & // ExpandScalarProduct;
  x = x /. GaugeXi[_] -> 3 /. SMP["g_s"] -> 1 /. SUNN -> 3;
  x = TID[x, l, ToPaVe -> True];
  x = PaXEvaluate[x, PaXImplicitPrefactor -> 1, PaXC0Expand -> True];
  x = x /. {ScaleMu -> 1, D -> 4 - 2 Epsilon};
  ser = Expand[Normal[Series[x, {Epsilon, 0, 0}]]];
  Print["PROJ ", v, " :: ", ToString[N[{Coefficient[ser, Epsilon, -2], Coefficient[ser, Epsilon, -1], Coefficient[ser, Epsilon, 0]}, 30], InputForm]],
  {v, {p2, p3}}];
