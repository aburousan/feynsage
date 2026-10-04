# One-loop features that match or extend Package-X 2.1.1: tensor reduction by the Passarino-Veltman
# recursion, raised propagator powers, exact derivatives, Taylor series, discontinuities and double
# spectral functions, pentagons, branch-aware logarithms and dilogarithms.  Reference numbers are from
# Package-X 2.1.1, LoopTools 2.16 (E0) or finite differences, as marked.
#     sage tests/test_px_parity.sage
import os, sys, time
sys.path.insert(0, os.path.abspath('.'))
from feynsage import *
from feynsage.pv import Kin, _Engine, _finalize, _C0f, _D0f, DiscB, mu

t_start = time.time()
fails = []


def close(name, a, b, tol=1e-12):
    a, b = CC(a), CC(b)
    err = abs(a - b) / max(abs(b), 1e-300)
    ok = err < tol
    print("%-46s %-4s rel %.1e" % (name, "ok" if ok else "FAIL", err))
    if not ok:
        fails.append(name)


def same(name, flag):
    print("%-46s %s" % (name, "ok" if flag else "FAIL"))
    if not flag:
        fails.append(name)


# 1. tensor coefficients: PV recursion = projection, exactly
def tensors_agree(vectors, kin, spec, r, symbols=()):
    out = []
    for use in (True, False):
        K = Kin(vectors, kin, symbols=symbols)
        eng = _Engine(K); eng.use_pv = use
        props = tuple((q, K.to_K(m)) for q, m in spec)
        T = eng.tensor(props, r)
        out.append({st: _finalize(c, K, False) for st, c in T.items()})
    A, B = out
    return all(bool((A.get(k, 0) - B.get(k, 0)).simplify_full() == 0) for k in set(A) | set(B))


same("PV recursion = projection: bubble rank 4",
     tensors_agree(["p"], {"p^2": "s"}, [({}, "m1"), ({"p": 1}, "m2")], 4, ["m1", "m2"]))
same("PV recursion = projection: triangle rank 3",
     tensors_agree(["p1", "p2"], {"p1^2": "s1", "p2^2": "s2", "p1.p2": "x"},
                   [({}, "m0"), ({"p1": 1}, "m1"), ({"p2": 1}, "m2")], 3, ["m0", "m1", "m2"]))
same("PV recursion = projection: box rank 3",
     tensors_agree(["p1", "p2", "p4"], {"p1^2": 0, "p2^2": 0, "p4^2": 0, "p1.p2": "s/2", "p1.p4": "t/2", "p2.p4": "-(s+t)/2"},
                   [({}, "m"), ({"p1": 1}, "m"), ({"p1": 1, "p2": 1}, "m"), ({"p4": -1}, "m")], 3, ["m"]))
t0 = time.time()
loop("l^a l^b l^c l^e", ["l", "m"], ["l+p1", "m"], ["l+p1+p2", "m"], ["l-p4", "m"],
     kin={"p1^2": 0, "p2^2": 0, "p4^2": 0, "p1.p2": "s/2", "p1.p4": "t/2", "p2.p4": "-(s+t)/2"})
same("rank-4 box in under 5 s (was > 200 s)", time.time() - t0 < 5)

# 2. raised powers (Package-X: Weights)
r = loop("1", ["l", "0"], ["l", "0"], ["l + p", "0"], kin={"p^2": 3}).coefficients()['1']
close("doubled massless bubble, finite part", finite_part(r, 1).n(80), CC(0.36620409622270323046, -1.0471975511965977462))
r = loop("1", ["l", "2"], ["l", "2"], ["l", "2"]).coefficients()['1']
close("tadpole cubed (Package-X -1/8)", r.n(), -0.125)
r = loop("1", ["l", "1"], ["l", "1"], ["l + p", "2"], kin={"p^2": 3}).coefficients()['1']
close("doubled massive bubble (finite difference)", finite_part(r, 1).n(80), -0.37355072789095, 1e-11)   # Package-X: Indeterminate

# 3. derivatives and series
x = var('x')
h = QQ(1) / 10**4
def fd(F, x0):
    g = lambda hh: (F(x0 + hh) - F(x0 - hh)) / (2 * hh)
    return (4 * g(h / 2) - g(h)) / 3
close("dDiscB/dm1 (finite difference)", loop_diff(DiscB(3, x, 2), x).subs(x=QQ(3)/2).n(80),
      fd(lambda v: DiscB(3, v, 2).n(digits=40), QQ(3)/2), 1e-10)
for i in range(6):
    a = [QQ(2), QQ(5), QQ(-3), QQ(1), QQ(2), QQ(3)]; x0 = a[i]; a[i] = x
    close("dC0/d(argument %d) (finite difference)" % i, loop_diff(_C0f(*a), x).subs(x=x0).n(80),
          fd(lambda v: _C0f(*[v if j == i else a[j] for j in range(6)]).n(digits=40), x0), 1e-11)
a = [QQ(7)/3, QQ(13)/5, QQ(3), QQ(17)/4, QQ(-5), QQ(-6), QQ(1), QQ(2), QQ(3), QQ(4)]
for i in (0, 5, 9):
    b = list(a); x0 = b[i]; b[i] = x
    close("dD0/d(argument %d) (finite difference)" % i, loop_diff(_D0f(*b), x).subs(x=x0).n(80),
          fd(lambda v: _D0f(*[v if j == i else b[j] for j in range(10)]).n(digits=40), x0), 1e-11)
s, m, m1, m2 = var('s m m1 m2')
from feynsage.pv import eps
same("B0(s,m,m) about s = 0 (Package-X: 1/6, 1/60)",
     bool((PVB(0, 0, s, m, m, series=(s, 0, 2)) - (1/eps + log(mu**2) - log(m**2) + s/(6*m**2) + s**2/(60*m**4))).subs(mu=2, m=3).simplify_full() == 0))
same("C0(s,s,0;m,m,m) about s = 0 (Package-X)",
     bool((PVC(0, 0, 0, s, s, 0, m, m, m, series=(s, 0, 1)) - (-1/(2*m**2) - s/(12*m**4))).simplify_full() == 0))
a_, b_, c_, e_, f_, g_ = var('a_ b_ c_ e_ f_ g_')
same("D0 about zero momenta (Package-X)",
     bool((PVD(0, 0, 0, 0, a_*x, b_*x, c_*x, e_*x, f_*x, g_*x, m, m, m, m, series=(x, 0, 1))
           - (1/(6*m**4) + (a_ + b_ + c_ + e_ + f_ + g_)*x/(60*m**6))).simplify_full() == 0))
B1s = PVB(0, 1, s, m1, m2, series=(s, 0, 1)).coefficient(s, 1)
pxB1 = -(2*m1**6 + 3*m1**4*m2**2 - 6*m1**2*m2**4 + m2**6 - 6*m1**4*m2**2*log(m1**2/m2**2)) / (6*(m1**2 - m2**2)**4)
same("B1 slope at s = 0 (Package-X)", bool((B1s - pxB1).expand_log().simplify_full() == 0))

# 4. discontinuities (Package-X Part -> Discontinuity)
t = var('t')
close("Disc_s B0(s;1,2) at s = 10", PVB(0, 0, s, 1, 2, disc=s).subs(s=10).n(), CC(0, 1.8849555921538759431))
close("Disc_s B1(s;1,2) at s = 10", PVB(0, 1, s, 1, 2, disc=s).subs(s=10).n(), CC(0, -0.65973445725385658008))
close("Disc_s C0(1,1,s;2,1,1) at s = 10", PVC(0, 0, 0, 1, 1, s, 2, 1, 1, disc=s).subs(s=10).n(), CC(0, -0.81661184686571930982))
close("Disc_s C1(3,-2,s;1,2,3) at s = 30", PVC(0, 1, 0, 3, -2, s, 1, 2, 3, disc=s).subs(s=30).n(), CC(0, 0.075036691516836366257))
close("Disc_s D0(1,2,3,4,s,-5;1,2,3,4) at s = 50", PVD(0, 0, 0, 0, 1, 2, 3, 4, s, -5, 1, 2, 3, 4, disc=s).subs(s=50).n(),
      CC(0, 0.012340404646077019714))
close("double spectral, equal masses", PVD(0, 0, 0, 0, 0, 0, 0, 0, s, t, 1, 1, 1, 1, disc=(s, t)).subs(s=10, t=8).n(),
      CC(0, 0.24836470664490253086))
close("double spectral, unequal masses", PVD(0, 0, 0, 0, 1, 2, 3, 4, s, t, 1, 2, 3, 4, disc=(s, t)).subs(s=50, t=60).n(),
      CC(0, 0.0065955577399639492024))

# 5. pentagons (LoopTools 2.16 E0 and Feynman parameters)
mink = lambda u, v: u[0]*v[0] - u[1]*v[1] - u[2]*v[2] - u[3]*v[3]
Q = [vector(QQ, v) for v in [(1/10, 1, 0, 0), (1/5, 3/10, 11/10, 0), (3/20, -2/5, 1/2, 9/10), (1/20, 3/5, -7/10, 2/5)]]
ms = [QQ(1), QQ(6)/5, QQ(3)/2, QQ(4)/5, QQ(11)/10]
nm = ["p1", "p2", "p3", "p4"]
kin5 = {(("%s^2" % nm[i]) if i == j else ("%s.%s" % (nm[i], nm[j]))): mink(Q[i], Q[j]) for i in range(4) for j in range(i, 4)}
pr5 = [["l", str(ms[0])]] + [["l + %s" % nm[i], str(ms[i + 1])] for i in range(4)]
close("pentagon E0 (LoopTools)", loop("1", *pr5, kin=kin5).coefficients()['1'].n(), -0.013835163216527903)
r1 = loop("l^mu", *pr5, kin=kin5).coefficients()
close("pentagon E_1 (Feynman parameters)", r1['p1^mu'].n(), 0.002610329175592079)
close("pentagon E_3 (Feynman parameters)", r1['p3^mu'].n(), 0.003182886783294718)
r2 = loop("l^mu l^nu", *pr5, kin=kin5).coefficients()
close("pentagon E_12, 4D basis (Feynman parameters)", r2['p1^mu p2^nu'].n(), -2.174507583563007, 1e-11)

# 6. helpers
close("Ln(-2, +)", Ln(-2, 1).n(80), CC(0.69314718055994530942, 3.1415926535897932385))
close("DiLog(3, -)", DiLog(3, -1).n(80), CC(2.3201804233130983964, -3.4513922952232026614))
close("ContinuedDiLog(-2-, -3-) (Package-X)", continued_dilog(-2, -1, -3, -1), CC(16.989929676117908948, 10.112396644223725364))
close("ContinuedDiLog(-2+, 3-) (Package-X)", continued_dilog(-2, 1, 3, -1), CC(1.2482731820994244579, -6.1132570288179918185))
close("DiscB(3,1,2) to 30 digits (Package-X)", DiscB(3, 1, 2).n(digits=32), -1.2091995761561452337293855051, 1e-28)
s1, s2, s3, s4, s12, s23 = var('s1 s2 s3 s4 s12 s23')
same("Kibble polynomial has 22 terms", kibble(s1, s2, s3, s4, s12, s23).expand().nops() == 22)

# 7. strings are read the same whatever the session defines
l = 1; p = 3
r = loop("1", ["l", "m"], ["l + p", "m"], kin={"p^2": "s"}).coefficients()['1']
same("a session variable l = 1 does not leak into loop()", r.has(SR.var('s')) and r.has(SR.var('m')))
del l, p

print("time %.1f s" % (time.time() - t_start))
print("ALL PACKAGE-X PARITY CHECKS PASS" if not fails else "FAILED: %s" % fails)
