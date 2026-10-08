# process().squared(average=False) against FeynArts + FeynCalc for 38 tree processes of the Standard Model
# (all masses kept, Feynman gauge, every spin, polarization and colour summed).  The FeynCalc results are
# stored in tests/feyncalc_ref/out_<name>.m.gz; tests/feyncalc_ref/oracle.wl makes them again (Mathematica).
#     sage tests/test_feyncalc.sage [nproc] [names...]
import sys, os, re, time, random, gzip
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
def _ev(txt):
    txt = re.sub(r'Sqrt\[', 'sqrt(', txt).replace(']', ')').replace('^', '**')
    return SR(sage_eval(txt, locals=loc))
def read_fc(path):
    lines = gzip.open(path, 'rt').read().strip().split('\n')
    if 'DEN' not in lines:
        return _ev(' '.join(lines))
    k = lines.index('DEN')
    num = sum((_ev(l) for l in lines[:k]), SR(0))
    return num / _ev(' '.join(lines[k + 1:]))
RF = RealField(200)
random.seed(int(11))
def point():
    sw = RF(random.uniform(0.3, 0.6)); mz = RF(random.uniform(85, 100))
    v = {loc['SW']: sw, loc['CW']: sqrt(1 - sw**2), loc['MZ']: mz, loc['MW']: mz*sqrt(1 - sw**2), loc['EL']: RF(random.uniform(0.2, 0.4)),
         loc['GS']: RF(random.uniform(0.8, 1.3)), loc['MH']: RF(random.uniform(110, 140))}
    for n in ('ME', 'MM', 'MU', 'MD', 'MC', 'MB'):
        v[loc[n]] = RF(random.uniform(1, 30))
    v[loc['MT']] = RF(random.uniform(200, 300))
    v[loc['s']] = RF(random.uniform(2e5, 9e5)); v[loc['t']] = -RF(random.uniform(1e3, 9e4))
    return v
res = []
for name, text in PROCS.items():
    if ONLY and name not in ONLY: continue
    path = os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), 'feyncalc_ref', 'out_%s.m.gz' % name)
    if not os.path.exists(path):
        print('MISS %-8s %-24s (no FeynCalc output)' % (name, text)); continue
    t0 = time.time()
    try:
        fc = read_fc(path)
        P = process(text)
        fs = P.squared(average=False, nproc=NP)
        m2 = sum(m**2 for m in P.masses)
        worst = 0
        for k in range(3):
            v = point()
            if len(P.masses) == 4:
                v[loc['u']] = m2.subs(v) - v[loc['s']] - v[loc['t']]
            a, b = RF(SR(fs).subs(v)), RF(SR(fc).subs(v))
            worst = max(worst, abs(a/b - 1))
        line = '%-4s %-8s %-24s %3d diagrams  rel.diff %.1e  %6.1f s' % ('ok' if worst < 1e-40 else 'DIFF', name, text, len(P.diagrams), float(worst), time.time() - t0)
    except Exception as ex:
        line = 'ERR  %-8s %-24s %s' % (name, text, repr(ex)[:200])
    print(line); sys.stdout.flush(); res.append(line)
nok = sum(l.startswith('ok') for l in res)
print('ALL FEYNCALC CHECKS PASS (%d processes)' % nok if nok == len(res) else 'FAILED: %d of %d' % (len(res) - nok, len(res)))
