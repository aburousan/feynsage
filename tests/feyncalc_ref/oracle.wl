(* Sum over all spins, polarizations and colours of |M|^2 at tree level in the FeynArts SM (Feynman gauge),
   computed by FeynCalc, written as InputForm in s, t, u for the comparison with feynsage (tests/test_feyncalc.sage).
   Run in this folder:  WolframKernel -script oracle.wl   (ORACLE_ONLY=ee_ww,gg_tt for a few), then gzip out_*.m.
   Photons are summed with -g, gluons with the physical sum (reference: another massless external momentum),
   massive vectors with -g + k k/M^2, as process() does. *)
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

run[name_, in_List, out_List] := Module[{all = Join[in, out], n, mom, ms, tops, diags, amp, sq, t0 = AbsoluteTime[], res},
  n = Length[all]; mom = Table[ToExpression["p" <> ToString[i]], {i, n}]; ms = mass /@ all;
  tops = CreateTopologies[0, Length[in] -> Length[out]];
  diags = InsertFields[tops, (f /@ in) -> (f /@ out), InsertionLevel -> {Particles}, Model -> SMQCD];
  amp = FCFAConvert[CreateFeynAmp[diags], IncomingMomenta -> mom[[;; Length[in]]], OutgoingMomenta -> mom[[Length[in] + 1 ;;]],
     UndoChiralSplittings -> True, ChangeDimension -> 4, List -> False, SMP -> True, Contract -> True, DropSumOver -> True,
     TransversePolarizationVectors -> {}];
  FCClearScalarProducts[];
  If[n == 4, SetMandelstam[s, t, u, mom[[1]], mom[[2]], -mom[[3]], -mom[[4]], Sequence @@ ms],
     ScalarProduct[mom[[1]], mom[[1]]] = ms[[1]]^2; ScalarProduct[mom[[2]], mom[[2]]] = ms[[2]]^2;
     ScalarProduct[mom[[3]], mom[[3]]] = ms[[3]]^2;
     ScalarProduct[mom[[2]], mom[[3]]] = (ms[[1]]^2 - ms[[2]]^2 - ms[[3]]^2)/2;
     ScalarProduct[mom[[1]], mom[[2]]] = (ms[[1]]^2 + ms[[2]]^2 - ms[[3]]^2)/2;
     ScalarProduct[mom[[1]], mom[[3]]] = (ms[[1]]^2 - ms[[2]]^2 + ms[[3]]^2)/2];
  sq = (amp ComplexConjugate[amp]) // FeynAmpDenominatorExplicit // FermionSpinSum // DiracSimplify //
       SUNSimplify[#, Explicit -> True, SUNNToCACF -> False] &;
  Do[If[MemberQ[vec, all[[i]]],
     sq = Which[all[[i]] === "gamma", DoPolarizationSums[sq, mom[[i]], 0],
       all[[i]] === "g", DoPolarizationSums[sq, mom[[i]],
          First[Select[Delete[Range[n], i], ms[[#]] === 0 &] /. {} -> {Mod[i, n] + 1}] /. j_Integer :> mom[[j]]],
       True, DoPolarizationSums[sq, mom[[i]]]]], {i, n}];
  sq = sq // FCE // DiracSimplify // SUNSimplify[#, SUNNToCACF -> False] & // ExpandScalarProduct;
  res = Together[sq /. rules];
  If[!FreeQ[res, Pair | Momentum | Eps | DiracGamma | Spinor], Print[name, ": UNREDUCED"]];
  Module[{num = Expand[Numerator[res]], den = Denominator[res]},
    Export["out_" <> name <> ".m", StringRiffle[Join[ToString[#, InputForm] & /@ If[Head[num] === Plus, List @@ num, {num}],
       {"DEN", ToString[den, InputForm]}], "\n"], "Text"]];
  Print[name, "  ", Length[diags], " diagrams  ", Round[AbsoluteTime[] - t0, 0.1], " s"];
];

procs = {
  {"ee_mumu", {"e-", "e+"}, {"mu-", "mu+"}}, {"ee_ee", {"e-", "e+"}, {"e-", "e+"}},
  {"ee_ww", {"e-", "e+"}, {"W+", "W-"}}, {"ee_zz", {"e-", "e+"}, {"Z", "Z"}}, {"ee_zh", {"e-", "e+"}, {"Z", "H"}},
  {"ee_az", {"e-", "e+"}, {"gamma", "Z"}}, {"ee_aa", {"e-", "e+"}, {"gamma", "gamma"}}, {"ee_nn", {"e-", "e+"}, {"nu_e", "nu_e~"}},
  {"ee_tt", {"e-", "e+"}, {"t", "t~"}}, {"ee_hh", {"e-", "e+"}, {"H", "H"}},
  {"ud_wa", {"u", "d~"}, {"W+", "gamma"}}, {"ud_wz", {"u", "d~"}, {"W+", "Z"}}, {"ud_wh", {"u", "d~"}, {"W+", "H"}},
  {"aa_ww", {"gamma", "gamma"}, {"W+", "W-"}}, {"ae_nw", {"gamma", "e-"}, {"nu_e", "W-"}},
  {"ww_zz", {"W+", "W-"}, {"Z", "Z"}}, {"ww_hh", {"W+", "W-"}, {"H", "H"}}, {"ww_ww", {"W+", "W-"}, {"W+", "W-"}},
  {"zz_hh", {"Z", "Z"}, {"H", "H"}}, {"hh_hh", {"H", "H"}, {"H", "H"}}, {"zh_zh", {"Z", "H"}, {"Z", "H"}},
  {"tt_hh", {"t", "t~"}, {"H", "H"}}, {"bb_zh", {"b", "b~"}, {"Z", "H"}}, {"en_en", {"e-", "nu_e~"}, {"e-", "nu_e~"}},
  {"uu_tt", {"u", "u~"}, {"t", "t~"}}, {"gg_tt", {"g", "g"}, {"t", "t~"}}, {"ug_uz", {"u", "g"}, {"u", "Z"}},
  {"ug_dw", {"u", "g"}, {"d", "W+"}}, {"uu_ga", {"u", "u~"}, {"g", "gamma"}}, {"gt_wb", {"g", "t"}, {"W+", "b"}},
  {"tb_wh", {"t", "b~"}, {"W+", "H"}}, {"uc_uc", {"u", "c"}, {"u", "c"}},
  {"h_tt", {"H"}, {"t", "t~"}}, {"z_hh", {"Z"}, {"e-", "e+"}}, {"t_bw", {"t"}, {"b", "W+"}}, {"h_ww", {"H"}, {"W+", "W-"}},
  {"w_ud", {"W+"}, {"u", "d~"}}, {"z_tt", {"Z"}, {"t", "t~"}}};
Do[If[only === All || MemberQ[only, p[[1]]],
   Check[run @@ p, Print[p[[1]], ": FAILED"]]], {p, procs}];
Print["ORACLE DONE"];
