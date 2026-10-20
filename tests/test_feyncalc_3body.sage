# process().squared(average=False) for 2 -> 3 and 1 -> 3 processes against FeynArts + FeynCalc, exactly, two ways:
#  formula: FeynCalc's whole result in the scalar products (feyncalc_ref/oracle_3body.wl with ORACLE_SYMBOLIC=1,
#           out3_<name>.m.gz) and feynsage's, as rational functions of the same independent invariants;
#  points:  rational phase-space points (points_3body.json from make_points.sage: every momentum, mass and coupling
#           rational, exactly on shell), where both sides are rational numbers (values_3body.json from oracle_3body.wl).
# Epsilon tensors are kept on both sides; at tree level their coefficient must vanish identically.
#     sage tests/test_feyncalc_3body.sage [nproc] [names]
import sys, os, re, time, json, gzip
from feynsage import process
from feynsage.qft import tensors as T
NP = int(sys.argv[1]) if len(sys.argv) > 1 else 16
ONLY = sys.argv[2:]
MASSLESS = {'uu_ggg': ['u'], 'ee_wwa': ['e'], 'ee_zha': ['e'], 'ud_wza': ['u', 'd']}
PROCS = {'ee_mma': 'e- e+ -> mu- mu+ gamma', 'ee_uug': 'e- e+ -> u u~ g', 'uu_ggg': 'u u~ -> g g g',
         'ee_wwa': 'e- e+ -> W+ W- gamma', 'ee_zha': 'e- e+ -> Z H gamma', 'ud_wza': 'u d~ -> W+ Z gamma',
         'mu_enn': 'mu- -> e- nu_e~ nu_mu', 't_benu': 't -> b e+ nu_e', 'h_eea': 'H -> e- e+ gamma'}
NAMES = ['EL', 'SW', 'CW', 'MZ', 'MH', 'ME', 'MM', 'MU', 'MD', 'MC', 'MT', 'MB', 'GS', 'E_eps'] + \
        ['d%d%d' % (i, j) for i in range(1, 6) for j in range(i + 1, 6)]
loc = {n: SR.var(n) for n in NAMES + ['MW']}
R = PolynomialRing(QQ, NAMES + ['fsI'])
IDEAL = R.ideal([R('SW^2 + CW^2 - 1'), R('fsI^2 + 1')])
REF = os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), 'feyncalc_ref')
pts = json.load(open(os.path.join(REF, 'points_3body.json')))
fcv = json.load(open(os.path.join(REF, 'values_3body.json')))

FQ = R.fraction_field()

def independent(P):
    """The products with the last leg and d_(n-2)(n-1) in terms of the others, from p_i.(sum p_in - sum p_out) = 0,
    and every EPS of the momenta as a multiple of E_eps = EPS(p1, ..., p4) (zero for a decay: three independent
    momenta).  Returns the scalar-product symbols for Dot(p_i, p_j), the eliminated ones and the EPS rule."""
    n, nin = len(P.momenta), len(P.p_in)
    sg = [1 if k < nin else -1 for k in range(n)]
    d = lambda i, j: P.masses[i]^2 if i == j else loc['d%d%d' % (min(i, j) + 1, max(i, j) + 1)]
    sol = solve([sum(sg[k]*d(i, k) for k in range(n)) == 0 for i in range(n)],
                [d(i, n - 1) for i in range(n - 1)] + [d(n - 3, n - 2)], solution_dict=True)[0]
    last = [-sg[k]*sg[n - 1] for k in range(n - 1)]             # p_n in terms of p_1 ... p_(n-1)
    vec = lambda q: last if str(q) == str(P.momenta[n - 1]) else [int(str(q) == str(b)) for b in P.momenta[:n - 1]]
    eps = lambda x: matrix(QQ, [vec(q) for q in x.operands()]).det()*loc['E_eps'] if n == 5 else 0
    dots = {T.dot(P.momenta[i], P.momenta[j]): d(i, j) for i in range(n) for j in range(i + 1, n)}
    return dots, sol, eps

sys.path.insert(0, os.path.dirname(os.path.abspath(sys.argv[0])))
from fc_text import evaluate, mathematica, items

def in_field(txt, sol):
    """A formula (Python syntax) evaluated in the fraction field over Q: the eliminated products, M_W = M_Z c_W
    and i put in."""
    names = {k: FQ(R(k)) for k in R.variable_names()}
    names.update({'MW': names['MZ']*names['CW'], 'I': names['fsI']})
    names.update({str(k): evaluate(str(v).replace('^', '**'), names) for k, v in sol.items()})
    return evaluate(txt, names)

def zero(f):
    """f = 0 in the fraction field, with s_W^2 + c_W^2 = 1 and i^2 = -1."""
    return IDEAL.reduce(f.numerator()) == 0

def formula_check(P, fs):
    dots, sol, eps = independent(P)
    a = SR(fs).subs(dots)
    a = a.subs({x: eps(x) for x in a.find(T.EPS(*[SR.wild(i) for i in range(4)]))})
    fc0, fce = [in_field(mathematica(x), sol) for x in items(gzip.open(os.path.join(REF, 'out3_%s.m.gz' % name), 'rt').read())]
    A = in_field(str(a).replace('^', '**'), sol)
    E = R('E_eps')
    A0, Ae = A.numerator().subs({E: 0})/A.denominator(), A.numerator().derivative(E)
    return zero(A0 - fc0) and zero(FQ(Ae)) and zero(fce)

def point_check(P, fs):
    n = len(P.momenta)
    for rec, (fc0, fce) in zip(pts[name]['points'], fcv[name]):
        v = {loc[k]: QQ(x) for k, x in rec.items() if not k.startswith('d') and k != 'moms'}
        dots = {T.dot(P.momenta[i], P.momenta[j]): QQ(rec['d%d%d' % (i + 1, j + 1)]) for i in range(n) for j in range(i + 1, n)}
        e = SR(fs).subs(dots).subs(v)
        comp = {P.momenta[i]: [QQ(x) for x in rec['moms'][i]] for i in range(n)}
        epses = e.find(T.EPS(*[SR.wild(i) for i in range(4)]))
        epsval = {x: matrix(QQ, [[comp[q][0], -comp[q][1], -comp[q][2], -comp[q][3]] for q in x.operands()]).det() for x in epses}
        a, a0 = QQ(e.subs(epsval)), QQ(e.subs({x: 0 for x in epses}))
        if not (a0 == QQ(fc0) and a == a0 and QQ(fce) == 0):
            return False
    return True

res = []
for name, text in PROCS.items():
    if ONLY and name not in ONLY:
        continue
    has_f, has_p = os.path.exists(os.path.join(REF, 'out3_%s.m.gz' % name)), name in fcv
    if not (has_f or has_p):
        print('MISS %-8s %-26s (no FeynCalc result yet)' % (name, text)); continue
    t0 = time.time()
    P = process(text, massless=MASSLESS.get(name, []))
    fs = P.squared(average=False, nproc=NP)
    f_ok = formula_check(P, fs) if has_f else None
    p_ok = point_check(P, fs) if has_p else None
    word = lambda ok, yes: 'not made yet' if ok is None else (yes if ok else 'DIFFERENT')
    line = '%-4s %-8s %-26s %3d diagrams  formula: %-12s  rational points: %-12s %6.1f s' % (
        'ok' if False not in (f_ok, p_ok) else 'DIFF', name, text, len(P.diagrams),
        word(f_ok, 'identical'), word(p_ok, 'equal'), time.time() - t0)
    print(line); sys.stdout.flush(); res.append(line)
nok = sum(l.startswith('ok') for l in res)
print('ALL 3-BODY FEYNCALC CHECKS PASS (%d processes)' % nok if res and nok == len(res) else 'FAILED')
