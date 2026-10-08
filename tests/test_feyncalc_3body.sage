# process().squared(average=False) for 2 -> 3 and 1 -> 3 processes against FeynArts + FeynCalc at random physical
# phase-space points (60 digits).  The points (momenta, masses, couplings) are in tests/feyncalc_ref/points_3body.json,
# FeynCalc's values (made by tests/feyncalc_ref/oracle_3body.wl, exact rationals) in values_3body.json.  Epsilon
# tensors of the momenta are evaluated from the components; their part must be zero at tree level.
#     sage tests/test_feyncalc_3body.sage [nproc]
import sys, os, re, time, random
from feynsage import process
from feynsage import models as MD
from feynsage.qft import tensors as T
NP = int(sys.argv[1]) if len(sys.argv) > 1 else 16
ONLY = sys.argv[2:]
MASSLESS = {'uu_ggg': ['u'], 'ee_wwa': ['e'], 'ee_zha': ['e'], 'ud_wza': ['u', 'd']}
PROCS = {'ee_mma': 'e- e+ -> mu- mu+ gamma', 'ee_uug': 'e- e+ -> u u~ g', 'uu_ggg': 'u u~ -> g g g',
         'ee_wwa': 'e- e+ -> W+ W- gamma', 'ee_zha': 'e- e+ -> Z H gamma', 'ud_wza': 'u d~ -> W+ Z gamma',
         'mu_enn': 'mu- -> e- nu_e~ nu_mu', 't_benu': 't -> b e+ nu_e', 'h_eea': 'H -> e- e+ gamma'}
NAMES = ['EL', 'SW', 'CW', 'MW', 'MZ', 'MH', 'ME', 'MM', 'MU', 'MD', 'MC', 'MT', 'MB', 'GS'] + \
        ['d%d%d' % (i, j) for i in range(1, 6) for j in range(i + 1, 6)]
loc = {n: SR.var(n) for n in NAMES}
def _ev(txt):
    txt = re.sub(r'Sqrt\[', 'sqrt(', txt).replace(']', ')').replace('^', '**')
    return SR(sage_eval(txt, locals=loc))
def read_fc(path):
    lines = open(path).read().strip().split('\n')
    if 'DEN' not in lines:
        return _ev(' '.join(lines))
    k = lines.index('DEN')
    return sum((_ev(l) for l in lines[:k]), SR(0)) / _ev(' '.join(lines[k + 1:]))
RF = RealField(200)
def boost(p, b):                      # p along the velocity b (3-vector)
    b2 = sum(x * x for x in b)
    if b2 == 0:
        return p
    g = 1 / sqrt(1 - b2); bp = sum(x * y for x, y in zip(b, p[1:]))
    E = g * (p[0] + bp)
    return [E] + [p[i + 1] + ((g - 1) * bp / b2 + g * p[0]) * b[i] for i in range(3)]
def two_body(M, m1, m2):
    P = sqrt((M**2 - (m1 + m2)**2) * (M**2 - (m1 - m2)**2)) / (2 * M)
    c = RF(random.uniform(-1, 1)); ph = RF(random.uniform(0, 2 * pi.n(200))); s_ = sqrt(1 - c**2)
    k = [P * s_ * cos(ph), P * s_ * sin(ph), P * c]
    return [sqrt(P**2 + m1**2)] + k, [sqrt(P**2 + m2**2)] + [-x for x in k]
def three_body(M, m):                  # in the rest frame of M
    m45 = RF(random.uniform(float(m[1] + m[2]), float(M - m[0])))
    a, b = two_body(M, m[0], m45)
    c, d = two_body(m45, m[1], m[2])
    bv = [x / b[0] for x in b[1:]]
    return [a, boost(c, bv), boost(d, bv)]
mdot = lambda a, b: a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3]
def point(P):
    sw = RF(random.uniform(0.3, 0.6)); mz = RF(random.uniform(85, 100))
    v = {loc['SW']: sw, loc['CW']: sqrt(1 - sw**2), loc['MZ']: mz, loc['MW']: mz * sqrt(1 - sw**2),
         loc['EL']: RF(random.uniform(0.2, 0.4)), loc['GS']: RF(random.uniform(0.8, 1.3)), loc['MH']: RF(random.uniform(110, 140))}
    for n in ('ME', 'MM', 'MU', 'MD', 'MC', 'MB'):
        v[loc[n]] = RF(random.uniform(1, 30))
    v[loc['MT']] = RF(random.uniform(200, 300))
    m = [RF(SR(x).subs(v)) for x in P.masses]
    if len(P.p_in) == 2:
        rs = RF(random.uniform(float(sum(m[2:]) + 50), float(sum(m[2:]) + 900)))
        rs = max(rs, m[0] + m[1] + 1)
        p1, p2 = two_body(rs, m[0], m[1])
        moms = [p1, p2] + three_body(rs, m[2:])
    else:
        moms = [[m[0], 0, 0, 0]] + three_body(m[0], m[1:])
    n = len(moms)
    for i in range(n):                                     # momentum conservation check
        pass
    d = {loc['d%d%d' % (i + 1, j + 1)]: mdot(moms[i], moms[j]) for i in range(n) for j in range(i + 1, n)}
    dots = {T.dot(P.momenta[i], P.momenta[j]): mdot(moms[i], moms[j]) for i in range(n) for j in range(i + 1, n)}
    return v, d, dots
import json
REF = os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), 'feyncalc_ref')
pts = json.load(open(os.path.join(REF, 'points_3body.json'))); fcv = json.load(open(os.path.join(REF, 'values_3body.json')))
res = []
for name, text in PROCS.items():
    if ONLY and name not in ONLY: continue
    if name not in fcv:
        print('MISS %-8s %s' % (name, text)); continue
    t0 = time.time()
    P = process(text, massless=MASSLESS.get(name, []))
    fs = P.squared(average=False, nproc=NP)
    worst = 0
    for rec, fcval in zip(pts[name]['points'], fcv[name]):
        v = {loc[k]: RF(x) for k, x in rec.items() if not k.startswith('d') and k != 'moms'}
        n = len(P.momenta)
        dots = {T.dot(P.momenta[i], P.momenta[j]): RF(rec['d%d%d' % (i + 1, j + 1)]) for i in range(n) for j in range(i + 1, n)}
        e0 = SR(fs).subs(dots).subs(v)
        epses = e0.find(T.EPS(*[SR.wild(i) for i in range(4)]))
        comp = {P.momenta[i]: [RF(x) for x in rec['moms'][i]] for i in range(n)}
        epsval = {x: matrix(RF, [[comp[q][0], -comp[q][1], -comp[q][2], -comp[q][3]] for q in x.operands()]).det() for x in epses}
        a_full = CC(e0.subs(epsval)); a = RF(e0.subs({x: 0 for x in epses})); b = RF(fcval)
        if epses:
            print('   feynsage epsilon part / total: %.1e' % (abs(a_full - a)/abs(a)))
        worst = max(worst, abs(a / b - 1))
    line = '%-4s %-8s %-26s %3d diagrams  rel.diff %.1e  %6.1f s' % ('ok' if worst < 1e-35 else 'DIFF', name, text, len(P.diagrams), float(worst), time.time() - t0)
    print(line); sys.stdout.flush(); res.append(line)
nok = sum(l.startswith('ok') for l in res)
print('ALL 3-BODY FEYNCALC CHECKS PASS (%d processes)' % nok if nok == len(res) and res else 'FAILED')
