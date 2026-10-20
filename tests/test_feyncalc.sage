# process().squared(average=False) against FeynArts + FeynCalc for 38 tree processes of the Standard Model,
# exactly: the two formulas are the same rational function (all masses kept, Feynman gauge, every spin, polarization and colour summed).  The FeynCalc results are
# stored in tests/feyncalc_ref/out_<name>.m.gz; tests/feyncalc_ref/oracle.wl makes them again (Mathematica).
#     sage tests/test_feyncalc.sage [nproc] [names...]
import sys, os, time, gzip
from feynsage import process
from feynsage import models as MD
NP = int(sys.argv[1]) if len(sys.argv) > 1 else 16
ONLY = sys.argv[2:]
PROCS = {
 'ee_mumu': 'e- e+ -> mu- mu+', 'ee_ee': 'e- e+ -> e- e+', 'ee_ww': 'e- e+ -> W+ W-', 'ee_zz': 'e- e+ -> Z Z',
 'ee_zh': 'e- e+ -> Z H', 'ee_az': 'e- e+ -> gamma Z', 'ee_aa': 'e- e+ -> gamma gamma', 'ee_nn': 'e- e+ -> nu_e nu_e~',
 'ee_tt': 'e- e+ -> t t~', 'ee_hh': 'e- e+ -> H H', 'ud_wa': 'u d~ -> W+ gamma', 'ud_wz': 'u d~ -> W+ Z',
 'ud_wh': 'u d~ -> W+ H', 'aa_ww': 'gamma gamma -> W+ W-', 'ae_nw': 'gamma e- -> nu_e W-', 'ww_zz': 'W+ W- -> Z Z',
 'ww_hh': 'W+ W- -> H H', 'ww_ww': 'W+ W- -> W+ W-', 'zz_hh': 'Z Z -> H H', 'hh_hh': 'H H -> H H',
 'zh_zh': 'Z H -> Z H', 'tt_hh': 't t~ -> H H', 'bb_zh': 'b b~ -> Z H', 'en_en': 'e- nu_e~ -> e- nu_e~',
 'uu_tt': 'u u~ -> t t~', 'gg_tt': 'g g -> t t~', 'ug_uz': 'u g -> u Z', 'ug_dw': 'u g -> d W+',
 'uu_ga': 'u u~ -> g gamma', 'gt_wb': 'g t -> W+ b', 'tb_wh': 't b~ -> W+ H', 'uc_uc': 'u c -> u c',
 'h_tt': 'H -> t t~', 'z_hh': 'Z -> e- e+', 't_bw': 't -> b W+', 'h_ww': 'H -> W+ W-', 'w_ud': 'W+ -> u d~', 'z_tt': 'Z -> t t~'}
NAMES = ['EL', 'SW', 'CW', 'MW', 'MZ', 'MH', 'ME', 'MM', 'MU', 'MD', 'MC', 'MT', 'MB', 'GS', 's', 't', 'u']
loc = {n: SR.var(n) for n in NAMES}
sys.path.insert(0, os.path.dirname(os.path.abspath(sys.argv[0])))
from fc_text import evaluate, mathematica
# Exact: both sides are rational functions of s, t, the masses and the couplings.  FeynCalc's text is evaluated
# straight in the fraction field over Q (u = sum m^2 - s - t, M_W = M_Z c_W put in); the numerator of the
# difference is reduced modulo s_W^2 + c_W^2 - 1.  Zero means the two are the same function.
EXACT_NAMES = ['EL', 'SW', 'CW', 'MZ', 'MH', 'ME', 'MM', 'MU', 'MD', 'MC', 'MT', 'MB', 'GS', 's', 't']
RQ = PolynomialRing(QQ, EXACT_NAMES)
FQ = RQ.fraction_field()
IDEAL = RQ.ideal([RQ('SW^2 + CW^2 - 1')])
def in_field(txt, u_val):
    names = {k: FQ(RQ(k)) for k in EXACT_NAMES}
    names['MW'] = names['MZ']*names['CW']
    names['u'] = evaluate(str(u_val).replace('^', '**'), names) if u_val != 0 else FQ(0)
    return evaluate(mathematica(txt), names)
def read_fc(path, u_val):
    lines = gzip.open(path, 'rt').read().strip().split('\n')
    if 'DEN' not in lines:
        return in_field(' '.join(lines), u_val)
    k = lines.index('DEN')
    return sum((in_field(l, u_val) for l in lines[:k]), FQ(0)) / in_field(' '.join(lines[k + 1:]), u_val)

res = []
for name, text in PROCS.items():
    if ONLY and name not in ONLY: continue
    path = os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), 'feyncalc_ref', 'out_%s.m.gz' % name)
    if not os.path.exists(path):
        print('MISS %-8s %-24s (no FeynCalc output)' % (name, text)); continue
    t0 = time.time()
    try:
        P = process(text)
        fs = P.squared(average=False, nproc=NP)
        u_val = sum(m**2 for m in P.masses) - loc['s'] - loc['t'] if len(P.masses) == 4 else 0
        exact = IDEAL.reduce((read_fc(path, u_val) - in_field(str(fs), u_val)).numerator()) == 0
        line = '%-4s %-8s %-24s %3d diagrams  %s  %6.1f s' % ('ok' if exact else 'DIFF', name, text, len(P.diagrams),
                                                             'identical' if exact else 'NOT IDENTICAL', time.time() - t0)
    except Exception as ex:
        line = 'ERR  %-8s %-24s %s' % (name, text, repr(ex)[:200])
    print(line); sys.stdout.flush(); res.append(line)
nok = sum(l.startswith('ok') for l in res)
print('ALL FEYNCALC CHECKS PASS, EXACTLY (%d processes)' % nok if nok == len(res) else 'FAILED: %d of %d' % (len(res) - nok, len(res)))
