(* Sum over all spins, polarizations and colours of |M|^2 at tree level in the FeynArts SM (Feynman gauge) for
   2 -> 3 and 1 -> 3 processes, by FeynCalc, for tests/test_feyncalc_3body.sage.  Run in this folder:
     ORACLE_ONLY=mu_enn WolframKernel -script oracle_3body.wl     -> val_mu_enn.json, the exact rational value at the
        rational points of points_3body.json (make_points.sage); collect the val_*.json into values_3body.json
     ORACLE_SYMBOLIC=1 ORACLE_ONLY=mu_enn ...                      -> out3_mu_enn.m, the whole formula in the dij
   Each result is {the part without epsilon tensors, the coefficient of Eps(p1, p2, p3, p4)}. *)
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

senv = Environment["ORACLE_SYMBOLIC"]; symbolic = StringQ[senv] && senv === "1";
pts = Import["points_3body.json", "RawJSON"];
massless = <|"uu_ggg" -> {SMP["m_u"]}, "ee_wwa" -> {SMP["m_e"]}, "ee_zha" -> {SMP["m_e"]}, "ud_wza" -> {SMP["m_u"], SMP["m_d"]}|>;
num[s_String] := ToExpression[s];
(* sum over spins, polarizations and colours of |M|^2 with the scalar products p_i.p_j = dots[i, j] (symbols or
   numbers); returns {the part without epsilon tensors, the coefficient of Eps(p1, p2, p3, p4)}.  The epsilon tensors
   are not dropped: the last momentum is put in by momentum conservation and EpsEvaluate writes each as a multiple of
   Eps(p1, p2, p3, p4) (zero for a decay, which has only three independent momenta); at tree level the coefficient
   must vanish, which the test checks *)
square[amp_, all_, mom_, nin_, ms_, dots_] := Module[{n = Length[all], sq, last},
  FCClearScalarProducts[];
  Do[ScalarProduct[mom[[i]], mom[[i]]] = ms[[i]]^2, {i, n}];
  Do[ScalarProduct[mom[[i]], mom[[j]]] = dots[i, j], {i, n}, {j, i + 1, n}];
  sq = (amp ComplexConjugate[amp]) // FeynAmpDenominatorExplicit // FermionSpinSum // DiracSimplify //
       SUNSimplify[#, Explicit -> True, SUNNToCACF -> False] &;
  Do[If[MemberQ[vec, all[[i]]],
     sq = Which[all[[i]] === "gamma", DoPolarizationSums[sq, mom[[i]], 0],
       all[[i]] === "g", DoPolarizationSums[sq, mom[[i]],
          First[Select[Delete[Range[n], i], ms[[#]] === 0 &] /. {} -> {Mod[i, n] + 1}] /. j_Integer :> mom[[j]]],
       True, DoPolarizationSums[sq, mom[[i]]]]], {i, n}];
  sq = sq // FCE // DiracSimplify // SUNSimplify[#, SUNNToCACF -> False] & // ExpandScalarProduct // Expand;
  last = Total[mom[[;; nin]]] - Total[mom[[nin + 1 ;; n - 1]]];
  sq = Expand[sq /. e_Eps :> EpsEvaluate[ExpandScalarProduct[FCReplaceMomenta[e, {mom[[n]] -> last}]]]];
  If[Length[Union[Cases[sq, _Eps, Infinity]]] > 1, Print["more than one epsilon tensor left: ", Union[Cases[sq, _Eps, Infinity]]]];
  {sq /. Eps[__] -> 0, Coefficient[sq /. Eps[__] -> EPSTENSOR, EPSTENSOR]}];

run[name_, in_List, out_List] := Module[{all = Join[in, out], n, mom, amp, ms, t0 = AbsoluteTime[], vals = {}, rec, rl, r},
  n = Length[all]; mom = Table[ToExpression["p" <> ToString[i]], {i, n}];
  amp = FCFAConvert[CreateFeynAmp[InsertFields[CreateTopologies[0, Length[in] -> Length[out]], (f /@ in) -> (f /@ out),
       InsertionLevel -> {Particles}, Model -> SMQCD]],
     IncomingMomenta -> mom[[;; Length[in]]], OutgoingMomenta -> mom[[Length[in] + 1 ;;]],
     UndoChiralSplittings -> True, ChangeDimension -> 4, List -> False, SMP -> True, Contract -> True, DropSumOver -> True,
     TransversePolarizationVectors -> {}];
  amp = amp /. ((# -> 0) & /@ Lookup[massless, name, {}]);
  ms = (mass /@ all) /. ((# -> 0) & /@ Lookup[massless, name, {}]);
  If[symbolic,
   (* ORACLE_SYMBOLIC=1: the whole formula in the symbols dij -> out3_<name>.m, a list {no-epsilon part, epsilon coefficient} *)
   r = square[amp, all, mom, Length[in], ms, ToExpression["d" <> ToString[#1] <> ToString[#2]] &] /. rules /. {SUNN -> 3};
   Print[name, "  left over: ", Union[Cases[r, _SMP | _Pair | _FeynAmpDenominator | _Momentum | _Eps, Infinity]]];
   WriteString["out3_" <> name <> ".m", ToString[r, InputForm]]; Close["out3_" <> name <> ".m"];
   Print[name, "  symbolic done  ", Round[AbsoluteTime[] - t0, 0.1], " s"];
   Return[]];
  Do[
   rec = pts[name]["points"][[k]];
   rl = {SMP["e"] -> num[rec["EL"]], SMP["sin_W"] -> num[rec["SW"]], SMP["cos_W"] -> num[rec["CW"]], SMP["m_W"] -> num[rec["MW"]],
      SMP["m_Z"] -> num[rec["MZ"]], SMP["m_H"] -> num[rec["MH"]], SMP["m_e"] -> num[rec["ME"]], SMP["m_mu"] -> num[rec["MM"]],
      SMP["m_u"] -> num[rec["MU"]], SMP["m_d"] -> num[rec["MD"]], SMP["m_c"] -> num[rec["MC"]], SMP["m_t"] -> num[rec["MT"]],
      SMP["m_b"] -> num[rec["MB"]], SMP["g_s"] -> num[rec["GS"]], SUNN -> 3, CA -> 3, CF -> 4/3};
   r = square[amp, all, mom, Length[in], ms, num[rec["d" <> ToString[#1] <> ToString[#2]]] &] /. rl;
   r = Together /@ r;
   If[!And @@ (NumericQ /@ r), Print["LEFT: ", Union[Cases[r, _Symbol | _SMP | _Pair | _Momentum, Infinity]]]];
   AppendTo[vals, ToString[#, InputForm, NumberMarks -> False] & /@ r],
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
