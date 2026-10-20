# Exact rational phase-space points for the 2 -> 3 and 1 -> 3 checks against FeynCalc.
# Every momentum component, mass and coupling is a rational number, every particle is exactly on its mass
# shell, momentum is conserved exactly, M_W = M_Z c_W and s_W^2 + c_W^2 = 1 hold exactly.  So feynsage and
# FeynCalc must give the same rational number: equality at a few random rational points of two rational
# functions is an exact identity test (Schwartz-Zippel), not a numerical one.
#     sage make_points.sage  ->  points_3body.json
# How the points are made:
#   unit vectors  n(a, b) = (2a, 2b, 1 - a^2 - b^2)/(1 + a^2 + b^2)
#   mass m, speed parameter t:  E = m (1 + t^2)/(1 - t^2),  |k| = 2 m t/(1 - t^2)
#   the last two outgoing particles: q = w n massless and f with mass m_f closing the momentum,
#   E_f^2 = m_f^2 + |S + w n|^2 is the hyperbola (E_f - w - b)(E_f + w + b) = m_f^2 + |S|^2 - b^2 (b = S.n),
#   rational in r = E_f - w - b.
import json, random, os
from feynsage import process
random.seed(int(20261009))
R = lambda lo, hi, den=60: QQ(random.randint(int(lo * den), int(hi * den))) / den

def unit():
    a, b = R(-2, 2), R(-2, 2)
    d = 1 + a^2 + b^2
    return vector(QQ, [2*a/d, 2*b/d, (1 - a^2 - b^2)/d])

def onshell(m, scale):
    """A rational momentum (E, k) with E^2 - k^2 = m^2, |k| of order scale."""
    n = unit()
    if m == 0:
        k = R(scale/4, scale)
        return [k] + list(k*n)
    t = R(QQ(1)/10, QQ(8)/10)
    return [m*(1 + t^2)/(1 - t^2)] + list(2*m*t/(1 - t^2)*n)

def close_up(S, mf):
    """q = w n (massless) and f (mass mf) with q + f = -S (spatial); returns q, f with rational energies."""
    while True:
        n = unit()
        b = S*n
        D = mf^2 + S*S - b^2
        if D <= 0:
            continue
        r = QQ(max(1, floor(RR(D).sqrt()*random.uniform(0.05, 0.9)*50)))/50      # 0 < r < sqrt(D)
        w = (D/r - r)/2 - b
        if w <= 0:
            continue
        Ef = (D/r + r)/2
        q = [w] + list(w*n)
        f = [Ef] + list(-(S + w*n))
        return q, f

MASSLESS = {'uu_ggg': ['u'], 'ee_wwa': ['e'], 'ee_zha': ['e'], 'ud_wza': ['u', 'd']}
PROCS = {'ee_mma': 'e- e+ -> mu- mu+ gamma', 'ee_uug': 'e- e+ -> u u~ g', 'uu_ggg': 'u u~ -> g g g',
         'ee_wwa': 'e- e+ -> W+ W- gamma', 'ee_zha': 'e- e+ -> Z H gamma', 'ud_wza': 'u d~ -> W+ Z gamma',
         'mu_enn': 'mu- -> e- nu_e~ nu_mu', 't_benu': 't -> b e+ nu_e', 'h_eea': 'H -> e- e+ gamma'}
mdot = lambda a, b: a[0]*b[0] - a[1]*b[1] - a[2]*b[2] - a[3]*b[3]

def point(P):
    t = R(QQ(16)/100, QQ(32)/100)                       # s_W between about 0.3 and 0.6
    v = {'SW': 2*t/(1 + t^2), 'CW': (1 - t^2)/(1 + t^2), 'MZ': R(85, 100, 4), 'EL': R(QQ(2)/10, QQ(4)/10),
         'GS': R(QQ(8)/10, QQ(13)/10), 'MH': R(110, 140, 4), 'MT': R(200, 300, 4)}
    v['MW'] = v['MZ']*v['CW']
    for nm in ('ME', 'MM', 'MU', 'MD', 'MC', 'MB'):
        v[nm] = R(1, 30, 4)
    names = [str(m) for m in P.masses]
    nin = len(P.p_in)
    out = list(range(nin, len(names)))
    # the massless outgoing leg taken as q, the last other leg as f, the rest free
    qi = next(i for i in reversed(out) if names[i] == '0')
    fi = next(i for i in reversed(out) if i != qi)
    free = [i for i in out if i not in (qi, fi)]
    # an incoming or parent mass is an output of the construction (fixed by the outgoing energies)
    fixed = set(names[i] for i in range(nin))
    val = lambda nm: QQ(0) if nm == '0' else v[nm]
    moms = {}
    for i in free:
        moms[i] = onshell(val(names[i]), 150)
    S = sum((vector(QQ, moms[i][1:]) for i in free), vector(QQ, [0, 0, 0]))
    moms[qi], moms[fi] = close_up(S, val(names[fi]))
    Et = sum(moms[i][0] for i in out)
    if nin == 1:
        v[names[0]] = Et
        moms[0] = [Et, 0, 0, 0]
    else:
        assert names[0] == names[1], "unequal incoming masses not handled"
        if names[0] == '0':
            moms[0], moms[1] = [Et/2, 0, 0, Et/2], [Et/2, 0, 0, -Et/2]
        else:
            u_ = R(QQ(3)/10, QQ(9)/10)                      # m = (E/2)(1 - u^2)/(1 + u^2), p_z = (E/2) 2u/(1 + u^2)
            m = Et/2*(1 - u_^2)/(1 + u_^2); pz = Et/2*2*u_/(1 + u_^2)
            v[names[0]] = m
            moms[0], moms[1] = [Et/2, 0, 0, pz], [Et/2, 0, 0, -pz]
    n = len(names)
    p = [moms[i] for i in range(n)]
    # exact checks of the kinematics
    for i in range(n):
        assert mdot(p[i], p[i]) == val(names[i])^2, (i, names[i])
    for c in range(4):
        assert sum(p[i][c] for i in range(nin)) == sum(p[i][c] for i in out)
    assert v['SW']^2 + v['CW']^2 == 1 and v['MW'] == v['MZ']*v['CW']
    rec = {k: str(x) for k, x in v.items()}
    for i in range(n):
        for j in range(i + 1, n):
            rec['d%d%d' % (i + 1, j + 1)] = str(mdot(p[i], p[j]))
    rec['moms'] = [[str(x) for x in q] for q in p]
    return rec

res = {}
for name, text in PROCS.items():
    P = process(text, massless=MASSLESS.get(name, []))
    res[name] = {'text': text, 'masses': [str(m) for m in P.masses], 'points': [point(P) for k in range(2)]}
    print(name, 'ok', max(len(x) for x in res[name]['points'][0].values() if isinstance(x, str)), 'characters at most')
here = os.path.dirname(os.path.abspath(sys.argv[0]))
json.dump(res, open(os.path.join(here, 'points_3body.json'), 'w'), indent=0)
print('written points_3body.json')
