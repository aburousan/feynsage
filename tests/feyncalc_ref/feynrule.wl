$FeynCalcStartupMessages = False; Quiet[Get["FeynCalc`"]]; $FCAdvice = False;
tos[x_] := ToString[FCE[x] //. {Pair[LorentzIndex[a_, ___], LorentzIndex[b_, ___]] :> Metric[a, b], MT[a_, b_] :> Metric[a, b],
   FV[p_, m_] :> Comp[p, m], Pair[LorentzIndex[m_, ___], Momentum[p_, ___]] :> Comp[p, m], SP[p_, q_] :> Dot[p, q],
   SUNF[a_, b_, c_] :> fcolor[a, b, c], SUNDelta[a_, b_] :> Kron[a, b], SMP["g_s"] -> gs}, InputForm];
YM = -(1/4) FieldStrength[al, be, i] . FieldStrength[al, be, i];
Print["G4: ", tos[FeynRule[YM, {QuantumField[GaugeField, {mu}, {a}][p1], QuantumField[GaugeField, {nu}, {b}][p2],
    QuantumField[GaugeField, {rho}, {c}][p3], QuantumField[GaugeField, {sg}, {d}][p4]}]]];
heft = -(1/4) CH FieldStrength[mu1, nu1, a1] . FieldStrength[mu1, nu1, a1] . QuantumField[H];
Print["HGG: ", tos[FeynRule[heft, {QuantumField[GaugeField, {mu}, {a}][p1], QuantumField[GaugeField, {nu}, {b}][p2], QuantumField[H][p3]}]]];
Print["HGGG: ", tos[FeynRule[heft, {QuantumField[GaugeField, {mu}, {a}][p1], QuantumField[GaugeField, {nu}, {b}][p2],
    QuantumField[GaugeField, {rho}, {c}][p3], QuantumField[H][p4]}]]];
