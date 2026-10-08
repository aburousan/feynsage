(* Sum over all spins, polarizations and colours of |M|^2 at tree level in the FeynArts SM (Feynman gauge) for
   2 -> 3 and 1 -> 3 processes, computed by FeynCalc at the phase-space points of points_3body.json: the scalar
   products and, after the traces, the masses and couplings are put in as exact rationals (60 digits).
   Run in this folder:  ORACLE_ONLY=mu_enn WolframKernel -script oracle_3body.wl  -> val_mu_enn.json;
   collect the val_*.json files into values_3body.json (with *^ written as e) for tests/test_feyncalc_3body.sage. *)
$FeynCalcStartupMessages = False; $LoadAddOns = {"FeynArts"};
Quiet[Get["FeynCalc`"]]; $FCAdvice = False;
env = Environment["ORACLE_ONLY"]; only = If[!StringQ[env] || env === "", All, StringSplit[env, ","]];
e[1] = F[2, {1}]; (* helpers *)
f = <|"e-" -> F[2, {1}], "e+" -> -F[2, {1}], "mu-" -> F[2, {2}], "mu+" -> -F[2, {2}],
     "nu_e" -> F[1, {1}], "nu_e~" -> -F[1, {1}], "nu_mu" -> F[1, {2}], "nu_mu~" -> -F[1, {2}],
     "u" -> F[3, {1}], "u~" -> -F[3, {1}], "d" -> F[4, {1}], "d~" -> -F[4, {1}], "c" -> F[3, {2}], "c~" -> -F[3, {2}],
     "t" -> F[3, {3}], "t~" -> -F[3, {3}], "b" -> F[4, {3}], "b~" -> -F[4, {3}],
     "gamma" -> V[1], "Z" -> V[2], "W-" -> V[3], "W+" -> -V[3], "g" -> V[5], "H" -> S[1]|>;
mass = <|"e-" -> SMP["m_e"], "e+" -> SMP["m_e"], "mu-" -> SMP["m_mu"], "mu+" -> SMP["m_mu"], "nu_e" -> 0, "nu_e~" -> 0,
     "nu_mu" -> 0, "nu_mu~" -> 0, "u" -> SMP["m_u"], "u~" -> SMP["m_u"], "d" -> SMP["m_d"], "d~" -> SMP["m_d"],
     "c" -> SMP["m_c"], "c~" -> SMP["m_c"], "t" -> SMP["m_t"], "t~" -> SMP["m_t"], "b" -> SMP["m_b"], "b~" -> SMP["m_b"],
     "gamma" -> 0, "Z" -> SMP["m_Z"], "W-" -> SMP["m_W"], "W+" -> SMP["m_W"], "g" -> 0, "H" -> SMP["m_H"]|>;
vec = {"gamma", "Z", "W-", "W+", "g"};
rules = {SMP["e"] -> EL, SMP["sin_W"] -> SW, SMP["cos_W"] -> CW, SMP["m_W"] -> MW, SMP["m_Z"] -> MZ, SMP["m_H"] -> MH,
   SMP["m_e"] -> ME, SMP["m_mu"] -> MM, SMP["m_u"] -> MU, SMP["m_d"] -> MD, SMP["m_c"] -> MC, SMP["m_t"] -> MT,
   SMP["m_b"] -> MB, SMP["g_s"] -> GS, SUNN -> 3, CA -> 3, CF -> 4/3};

pts = Import["points_3body.json", "RawJSON"];
massless = <|"uu_ggg" -> {SMP["m_u"]}, "ee_wwa" -> {SMP["m_e"]}, "ee_zha" -> {SMP["m_e"]}, "ud_wza" -> {SMP["m_u"], SMP["m_d"]}|>;
num[s_String] := Rationalize[ToExpression[s <> "`60"], 10^-58];
run[name_, in_List, out_List] := Module[{all = Join[in, out], n, mom, tops, diags, amp, sq, t0 = AbsoluteTime[], vals = {}, rec, rl, ms},
  n = Length[all]; mom = Table[ToExpression["p" <> ToString[i]], {i, n}];
  tops = CreateTopologies[0, Length[in] -> Length[out]];
  diags = InsertFields[tops, (f /@ in) -> (f /@ out), InsertionLevel -> {Particles}, Model -> SMQCD];
  amp = FCFAConvert[CreateFeynAmp[diags], IncomingMomenta -> mom[[;; Length[in]]], OutgoingMomenta -> mom[[Length[in] + 1 ;;]],
     UndoChiralSplittings -> True, ChangeDimension -> 4, List -> False, SMP -> True, Contract -> True, DropSumOver -> True,
     TransversePolarizationVectors -> {}];
  amp = amp /. ((# -> 0) & /@ Lookup[massless, name, {}]);
  Do[
   rec = pts[name]["points"][[k]];
   rl = Join[{SMP["e"] -> num[rec["EL"]], SMP["sin_W"] -> num[rec["SW"]], SMP["cos_W"] -> num[rec["CW"]], SMP["m_W"] -> num[rec["MW"]],
      SMP["m_Z"] -> num[rec["MZ"]], SMP["m_H"] -> num[rec["MH"]], SMP["m_e"] -> num[rec["ME"]], SMP["m_mu"] -> num[rec["MM"]],
      SMP["m_u"] -> num[rec["MU"]], SMP["m_d"] -> num[rec["MD"]], SMP["m_c"] -> num[rec["MC"]], SMP["m_t"] -> num[rec["MT"]],
      SMP["m_b"] -> num[rec["MB"]], SMP["g_s"] -> num[rec["GS"]], SUNN -> 3, CA -> 3, CF -> 4/3}];
   ms = (mass /@ all) /. ((# -> 0) & /@ Lookup[massless, name, {}]);
   FCClearScalarProducts[];
   Do[ScalarProduct[mom[[i]], mom[[i]]] = ms[[i]]^2, {i, n}];
   Do[ScalarProduct[mom[[i]], mom[[j]]] = num[rec["d" <> ToString[i] <> ToString[j]]], {i, n}, {j, i + 1, n}];
   sq = (amp ComplexConjugate[amp]) // FeynAmpDenominatorExplicit // FermionSpinSum // DiracSimplify //
        SUNSimplify[#, Explicit -> True, SUNNToCACF -> False] &;
   Do[If[MemberQ[vec, all[[i]]],
      sq = Which[all[[i]] === "gamma", DoPolarizationSums[sq, mom[[i]], 0],
        all[[i]] === "g", DoPolarizationSums[sq, mom[[i]],
           First[Select[Delete[Range[n], i], ms[[#]] === 0 &] /. {} -> {Mod[i, n] + 1}] /. j_Integer :> mom[[j]]],
        True, DoPolarizationSums[sq, mom[[i]]]]], {i, n}];
   sq = sq // FCE // DiracSimplify // SUNSimplify[#, SUNNToCACF -> False] & // ExpandScalarProduct // Expand;
   Module[{x = N[sq /. rl /. {SUNN -> 3}, 50], ep}, ep = Coefficient[Expand[x], Eps[__]]; x = Expand[x] /. Eps[__] -> 0; Print["EPS coefficient size: ", Max[Abs[N[Cases[Expand[N[sq /. rl /. {SUNN -> 3}, 50]], c_ * Eps[__] :> c, {0, 1}]]]] ]; If[!NumericQ[x], Print["LEFT: ", ToString[Union[Cases[x, _Symbol | _SMP | _Pair | _FeynAmpDenominator | _Momentum, Infinity]], InputForm]]];
     AppendTo[vals, ToString[x, InputForm, NumberMarks -> False]]],
   {k, Length[pts[name]["points"]]}];
  Export["val_" <> name <> ".json", vals, "JSON"];
  Print[name, "  ", vals, "  ", Round[AbsoluteTime[] - t0, 0.1], " s"];
];

procs = {
  {"ee_mma", {"e-", "e+"}, {"mu-", "mu+", "gamma"}}, {"ee_uug", {"e-", "e+"}, {"u", "u~", "g"}},
  {"uu_ggg", {"u", "u~"}, {"g", "g", "g"}}, {"ee_wwa", {"e-", "e+"}, {"W+", "W-", "gamma"}},
  {"ee_zha", {"e-", "e+"}, {"Z", "H", "gamma"}}, {"ud_wza", {"u", "d~"}, {"W+", "Z", "gamma"}},
  {"mu_enn", {"mu-"}, {"e-", "nu_e~", "nu_mu"}}, {"t_benu", {"t"}, {"b", "e+", "nu_e"}}, {"h_eea", {"H"}, {"e-", "e+", "gamma"}}};
Do[If[only === All || MemberQ[only, p[[1]]],
   Check[run @@ p, Print[p[[1]], ": FAILED"]]], {p, procs}];
Print["ORACLE DONE"];
