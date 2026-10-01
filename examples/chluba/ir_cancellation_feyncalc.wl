(* The same check in FeynCalc: one-loop QED vertex (Feynman gauge, D = 4 - 2 eps), tensor reduction
   with TID, scalar integrals from Package-X through FeynHelpers (PaXEvaluate).  The 1/eps pole of
   F1 changes with t exactly like the soft-photon angular integral I(t); F2 -> alpha/2pi as t -> 0. *)
$LoadAddOns = {"FeynHelpers"};
Get["FeynCalc`"];
outfile = FileNameJoin[{Directory[], "slides", "out", "chluba_ir_feyncalc.txt"}];
str = OpenWrite[outfile];
say[x_] := (WriteString[str, x <> "\n"]; Print[x]);

ffactors[tv_] := Module[{numG, numP, loopG, loopP, tree, M, AB, F1, F2},
  ClearScalarProducts[];
  SPD[p, p] = 1; SPD[pp, pp] = 1; SPD[p, pp] = 1 - tv/2;
  numG = DiracTrace[(GSD[p] + 1).GAD[mu].(GSD[pp] + 1).GAD[nu].(GSD[pp + l] + 1).GAD[mu].(GSD[p + l] + 1).GAD[nu],
     DiracTraceEvaluate -> True];
  numP = DiracTrace[(GSD[p] + 1).(GSD[pp] + 1).GAD[nu].(GSD[pp + l] + 1).(GSD[p] + GSD[pp]).(GSD[p + l] + 1).GAD[nu],
     DiracTraceEvaluate -> True];
  loopG = PaXEvaluate[ToPaVe[TID[ExpandScalarProduct[Contract[numG]] FAD[{l, 0}, {l + pp, 1}, {l + p, 1}], l], l]/(I Pi^2),
     PaXImplicitPrefactor -> 1];
  loopP = PaXEvaluate[ToPaVe[TID[ExpandScalarProduct[Contract[numP]] FAD[{l, 0}, {l + pp, 1}, {l + p, 1}], l], l]/(I Pi^2),
     PaXImplicitPrefactor -> 1];
  tree[x_] := ExpandScalarProduct[Contract[DiracTrace[x, DiracTraceEvaluate -> True]]];
  M = {{tree[(GSD[p] + 1).GAD[mu].(GSD[pp] + 1).GAD[mu]], tree[(GSD[p] + 1).(GSD[p] + GSD[pp]).(GSD[pp] + 1)]/2},
       {tree[(GSD[p] + 1).(GSD[pp] + 1).(GSD[p] + GSD[pp])], tree[(GSD[p] + 1).(GSD[pp] + 1)] (4 - tv)/2}} /. D -> 4 - 2 Epsilon;
  AB = LinearSolve[M, {loopG, loopP}];
  F1 = Normal[Series[AB[[1]] + AB[[2]] /. ScaleMu -> 1, {Epsilon, 0, 0}]];
  F2 = Normal[Series[-AB[[2]] /. ScaleMu -> 1, {Epsilon, 0, 0}]];
  {Coefficient[Expand[F1], Epsilon, -1], F2 /. Epsilon -> 0}];

Iang[tv_] := 2 (1 - tv/2) NIntegrate[1/(1 - x (1 - x) tv), {x, 0, 1}, WorkingPrecision -> 30] - 2;
{tF, res} = AbsoluteTiming[Table[tv -> ffactors[tv], {tv, {-1/2, -1, -3, -8}}]];
p0 = N[(-1/2 /. res)[[1]], 20];
Do[say["t = " <> ToString[tv] <> ":  pole(t) - pole(-1/2) = " <> ToString[N[(tv /. res)[[1]] - p0, 15]] <>
     "   I(t) - I(-1/2) = " <> ToString[N[Iang[tv] - Iang[-1/2], 15]]], {tv, {-1, -3, -8}}];
say["F2 at t = -1/10^6 in units alpha/(4 pi): " <> ToString[N[ffactors[-1/10^6][[2]], 12]] <> "  (expected 2)"];
say["time " <> ToString[Round[tF, 0.1]] <> " s"];
Export[FileNameJoin[{Directory[], "slides", "out", "chluba_ir_feyncalc.json"}],
  <|"lines" -> Table[{ToString[tv], N[(tv /. res)[[1]] - p0, 15], N[Iang[tv] - Iang[-1/2], 15]}, {tv, {-1, -3, -8}}]|>, "RawJSON"];
Close[str];
