# Compton and double Compton scattering from the Feynman diagrams, with feynsage + FORM.
# Checks the results of J. Chluba's thesis "Spectral distortions of the CMB" (2005):
#   Compton: Klein-Nishina; double Compton: Mandl & Skyrme's X (thesis eq. D.1, |M|^2 = e^6 X);
#   the soft limit of double Compton (the 1/nu_2 infrared divergence).
# Kinematics: e(P) + gamma(K0) -> e(P') + gamma(K1) + gamma(K2), m = 1, exact rationals.
import sys, json, time, itertools, random; sys.path.insert(0, '.')
from feynsage import form
random.seed(7r)
O = 'examples/chluba/out/'
dot_ = function('dot')
m = var('m')

def evaluate(expr, table):
    return expr.substitute_function(dot_, lambda a, b: table[tuple(sorted((str(a), str(b))))]).subs(m=1).expand()

# ------------------------------------------------------------------ Compton
def compton_msq(table, pk, pkp):
    """sum over spins and polarisations of |M|^2 / e^4 for e(p) gamma(k) -> e(pp) gamma(kp)."""
    s_m, u_m = 2 * pk, -2 * pkp             # s - m^2 = 2 p.k,  u - m^2 = -2 p.kp
    terms = [(['nu', 'p + k + m', 'mu'], ['mu', 'p + k + m', 'nu'], s_m),
             (['mu', 'p - kp + m', 'nu'], ['nu', 'p - kp + m', 'mu'], u_m)]
    tot = 0
    for a, an, da in terms:
        for b, bn, db in terms:
            tr = form.dirac_trace(['pp + m'] + a + ['p + m'] + bn, vectors=['p', 'pp', 'k', 'kp'])
            tot += evaluate(tr, table) / (da * db)
    return tot                              # (-g)^2 = +1 for two photon polarisation sums

def table_compton(pk, pkp):
    kkp = pk - pkp                           # from (p + k - kp)^2 = m^2
    t = {('p', 'p'): 1, ('k', 'k'): 0, ('kp', 'kp'): 0, ('k', 'p'): pk, ('kp', 'p'): pkp, ('k', 'kp'): kkp}
    # pp = p + k - kp
    t[('p', 'pp')] = 1 + pk - pkp; t[('k', 'pp')] = pk - kkp; t[('kp', 'pp')] = pkp + kkp - 0
    t[('pp', 'pp')] = 1
    return {tuple(sorted(k)): v for k, v in t.items()}

res = {}
t0 = time.time()
pk, pkp = QQ(3) / 7, QQ(2) / 9
M2 = compton_msq(table_compton(pk, pkp), pk, pkp)
s_m, u_m = 2 * pk, -2 * pkp
peskin = 8 * (-u_m / s_m - s_m / u_m + 4 * (1 / s_m + 1 / u_m) + 4 * (1 / s_m + 1 / u_m) ** 2)   # 4 x Peskin (5.87), m = 1
res['compton_sum_vs_textbook'] = str(M2 - peskin)
print("Compton: sum |M|^2/e^4 =", M2, "  textbook:", peskin, "  difference:", M2 - peskin)

# Klein-Nishina: lab frame, omega, theta; dsigma/dOmega = (omega'/omega)^2 <|M|^2>/(64 pi^2 m^2)
w, c = var('w c')                            # omega/m and cos(theta)
wp = w / (1 + w * (1 - c))
avg = (8 * (-(-2 * wp) / (2 * w) - (2 * w) / (-2 * wp) + 4 * (1 / (2 * w) + 1 / (-2 * wp)) + 4 * (1 / (2 * w) + 1 / (-2 * wp)) ** 2)) / 4
e2 = 4 * pi * var('alpha')
dsig = (wp / w) ** 2 * avg * e2 ** 2 / (64 * pi ** 2)
r0 = var('alpha')                            # r0 = alpha/m with m = 1
kn = r0 ** 2 / 2 * (wp / w) ** 2 * (wp / w + w / wp - (1 - c ** 2))
res['klein_nishina'] = bool((dsig - kn).simplify_full() == 0)
sigT = integrate(kn.subs(w=0) * 2 * pi, c, -1, 1)
res['thomson'] = str(sigT)
print("Klein-Nishina reproduced:", res['klein_nishina'], "   Thomson limit sigma =", sigT, "(= 8 pi r0^2/3)")

# ------------------------------------------------------------------ double Compton
V = ['p', 'pp', 'k0', 'k1', 'k2']
def dc_table(pk0, pk1, pk2, k0k1, k0k2):
    # on-shell P' = P + K0 - K1 - K2 fixes k1.k2
    k1k2 = -(pk0 - pk1 - pk2 - k0k1 - k0k2)
    t = {('p', 'p'): 1, ('k0', 'k0'): 0, ('k1', 'k1'): 0, ('k2', 'k2'): 0,
         ('k0', 'p'): pk0, ('k1', 'p'): pk1, ('k2', 'p'): pk2, ('k0', 'k1'): k0k1, ('k0', 'k2'): k0k2, ('k1', 'k2'): k1k2}
    d = lambda a, b: t[tuple(sorted((a, b)))]
    pp = {'p': 1, 'k0': 1, 'k1': -1, 'k2': -1}
    for v in ['p', 'k0', 'k1', 'k2']:
        t[tuple(sorted(('pp', v)))] = sum(cc * d(u, v) for u, cc in pp.items())
    t[('pp', 'pp')] = sum(c1 * c2 * d(u1, u2) for u1, c1 in pp.items() for u2, c2 in pp.items())
    return {tuple(sorted(k)): v for k, v in t.items()}

def dc_msq(table):
    """sum over spins and polarisations of |M|^2 / e^6, all 6 diagrams (the 3! orders of the photons)."""
    d = lambda a, b: table[tuple(sorted((a, b)))]
    q = {'a': ('-k0', 'mu1'), 'b': ('k1', 'mu2'), 'c': ('k2', 'mu3')}  # photons as outgoing momenta
    chains = []
    for o in itertools.permutations('abc'):
        (q1, i1), (q2, i2), (q3, i3) = q[o[0]], q[o[1]], q[o[2]]
        # u(P) first emits o[2], then o[1], then o[0]
        prop1 = 'pp + %s' % q1                              # P' + q_sigma1
        prop2 = 'p - (%s)' % q3                             # P - q_sigma3
        chains.append((o, [i1, prop1 + ' + m', i2, prop2 + ' + m', i3], prop1, prop2))
    def denom(expr):
        # momentum as linear combination over V
        coeffs = {}
        e = SR(expr.replace('pp', 'PPx'))
        for v in V:
            coeffs[v] = e.coefficient(SR('PPx' if v == 'pp' else v))
        val = sum(coeffs[a] * coeffs[b] * d(a, b) for a in V for b in V)
        return val - 1
    for v in V:
        SR.var(v)
    SR.var('PPx')
    tot = 0
    for o1, ch1, a1, b1 in chains:
        for o2, ch2, a2, b2 in chains:
            fac = ['pp + m'] + ch1 + ['p + m'] + list(reversed(ch2))
            tr = form.dirac_trace(fac, vectors=V)
            tot += evaluate(tr, table) / (denom(a1) * denom(b1) * denom(a2) * denom(b2))
    return -tot                              # (-g)^3 for three photon polarisation sums

def mandl_skyrme_X(table):
    d = lambda a, b: table[tuple(sorted((a, b)))]
    k0, k1, k2 = -d('p', 'k0'), d('p', 'k1'), d('p', 'k2')
    k0p, k1p, k2p = d('pp', 'k0'), -d('pp', 'k1'), -d('pp', 'k2')
    a = 1 / k0 + 1 / k1 + 1 / k2; b = 1 / k0p + 1 / k1p + 1 / k2p
    c = 1 / (k0 * k0p) + 1 / (k1 * k1p) + 1 / (k2 * k2p)
    x = k0 + k1 + k2; z = k0 * k0p + k1 * k1p + k2 * k2p
    A = k0 * k1 * k2; B = k0p * k1p * k2p
    rho = k0 / k0p + k0p / k0 + k1 / k1p + k1p / k1 + k2 / k2p + k2p / k2
    return (2 * (a * b - c) * ((a + b) * (2 + x) - (a * b - c) - 8) - 2 * x * (a ** 2 + b ** 2)
            - 2 * (a * b + c * (1 - x)) * rho - 8 * c
            + 4 * x / (A * B) * ((A + B) * (1 + x) + x ** 2 * (1 - z) + 2 * z - (a * A + b * B) * (2 + (1 - x) * z / x)))

ratios, points = [], []
for trial in range(4):
    r = lambda: QQ(random.randint(1r, 40r)) / random.randint(5r, 30r)
    pt = (r(), r(), r(), r(), r())
    points.append([str(v) for v in pt])
    tab = dc_table(*pt)
    t1 = time.time()
    m2 = dc_msq(tab)
    X = mandl_skyrme_X(tab)
    ratios.append(m2 / X)
    print("double Compton point %d: sum|M|^2/e^6 = %s,  X = %s,  ratio %s  (%.1f s)" % (trial, m2, X, m2 / X, time.time() - t1))
res['dc_ratios'] = [str(x) for x in ratios]
res['dc_points'] = points

# ------------------------------------------------------------------ soft limit
# scale K2 -> lam K2 (keeping P' on shell): sum|M_DC|^2 -> e^2 S(K2) sum|M_C|^2, S the eikonal factor
base = (QQ(3) / 4, QQ(1) / 3, QQ(1) / 2, QQ(3) / 4 - QQ(1) / 3, QQ(1) / 3)   # p.k0, p.k1, p.k2, k0.k1 = p.k0 - p.k1, k0.k2
soft = []
for lam in (QQ(1) / 10 ** 4, QQ(1) / 10 ** 6, QQ(1) / 10 ** 8):
    pk0, pk1, pk2, k0k1, k0k2 = base
    tab = dc_table(pk0, pk1, lam * pk2, k0k1 - lam * QQ(1) / 5, lam * k0k2)   # P'.K2 = lam/5
    d = lambda a, b: tab[tuple(sorted((a, b)))]
    S = 2 * d('p', 'pp') / (d('p', 'k2') * d('pp', 'k2')) - 1 / d('p', 'k2') ** 2 - 1 / d('pp', 'k2') ** 2
    # Compton part: e(P) gamma(K0) -> e(P') gamma(K1) with P' = P + K0 - K1 - K2 (K2 -> 0)
    tc = table_compton(pk0, pk1)
    mc = compton_msq(tc, pk0, pk1)
    ratio = dc_msq(tab) / (S * mc)
    soft.append((float(lam), float(ratio)))
    print("soft limit, K2 scaled by %g: sum|M_DC|^2 / (S sum|M_C|^2) = %.10f" % (lam, ratio))
res['soft'] = soft
res['time'] = time.time() - t0
json.dump(res, open(O + 'chluba_feynsage.json', 'w'), indent=1)
