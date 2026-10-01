# The code cells of the feynsage tutorial, in order.  Each is (id, code, options); the cells
# share one namespace, like a notebook.  Options: fig=<name> saves that matplotlib figure,
# text=True shows the value as plain text.
CELLS = [

# ------------------------------------------------------------------ 1. first steps
("import", """
from feynsage import *
""", {}),

("info", """
info(A0)
""", {}),

("a0", """
m = var('m')
a = A0(m)
a
""", {}),

("a0parts", """
print("pole:           ", uv_part(a))
print("finite part:    ", finite_part(a))
print("finite, mu = m: ", finite_part(a, m))
""", {}),

("b0", """
s = var('s')
B0(s, m, m)
""", {}),

("b0numbers", """
for sv in [1, 3, 5, 10]:
    val = finite_part(B0(sv, 1, 1), 1).n(digits=12)
    beta = sqrt(1 - 4/sv) if sv > 4 else 0
    print("s = %2d   B0 finite = %-38s pi*beta = %.12f" % (sv, val, (pi*beta).n()))
""", {}),

("b0plot", """
from feynsage.plotting import quick_plot, set_theme
set_theme()
fig = quick_plot([B0(s, 1, 1), B0(s, 1, 2)], (s, -2, 14),
                 labels=["$m_1 = m_2 = 1$", "$m_1 = 1,\\\\ m_2 = 2$"])
""", {"fig": "fig"}),

("a0b0", """
(A0(m) - m^2*(B0(0, m, m) + 1)).simplify_full()
""", {}),

("norm", """
from feynsage.oneloop import tadpole, expand_eps, D
note = expand_eps(tadpole(1, 1), 0)        # the note: e^(eps gamma_E) Int [dl] 1/(l^2 + 1), Euclidean
px = A0(1).subs(mu=1)                       # Package-X normalisation, Minkowski, m = mu = 1
print("note's tadpole:", note)
print("A0(1):         ", px)
print("sum:           ", (note + px).simplify_full())
""", {}),

# ------------------------------------------------------------------ 2. one loop with numerators
("loop1", """
loop("1", ["l", "m1"], ["l + p", "m2"], kin={"p^2": "s"})
""", {}),

("loopmu", """
loop("l^mu", ["l", "m1"], ["l + p", "m2"], kin={"p^2": "s"})
""", {}),

("b1check", """
sv, a1, a2 = 10, 1, 2
hand = (A0(a1) - A0(a2) - (sv + a1^2 - a2^2)*B0(sv, a1, a2))/(2*sv)
pv = PVB(0, 1, sv, a1, a2)                      # Package-X's PVB[0, 1, ...] = B1
print("poles:        ", uv_part(hand), uv_part(pv))
print("finite parts: ", finite_part(hand, 1).n(), "   ", finite_part(pv, 1).n())
""", {}),

("loopmunu", """
loop("l^mu l^nu", ["l", "m"], ["l + p", "m"], kin={"p^2": "s"})
""", {}),

("gram", """
s1, s2, s12 = var('s1 s2 s12')
p1p2 = (s1 + s2 - s12)/2                    # p1.p2 from s12 = (p1 - p2)^2
Gram = matrix([[s1, p1p2], [p1p2, s2]])
print("det Gram =", Gram.det().expand())
print("at the g-2 point s1 = s2 = m^2, q^2 = s12 = 0:", Gram.det().subs(s1=m^2, s2=m^2, s12=0))
""", {}),

("c0", """
c = C0(1, 2, 3, 1, 2, 3)
print(c, "=", c.n(digits=30))
print("direct numerical integration:", c0_numeric(1, 2, 3, 1, 2, 3))
print("as logs and dilogs:", len(str(explicit(c))), "characters, value", explicit(c).n(digits=20).real())
""", {}),

("c0hand", """
from scipy.integrate import dblquad
s1v, s12v, s2v, m0, m1, m2 = -1.0, -2.0, -3.0, 1.0, 2.0, 3.0      # all invariants spacelike: F > 0
def Fpar(x1, x2):
    x0 = 1 - x1 - x2                         # x0, x1, x2 on the lines (l, m0), (l + p1, m1), (l + p2, m2)
    return x0*m0**2 + x1*m1**2 + x2*m2**2 - x0*x1*s1v - x0*x2*s2v - x1*x2*s12v
by_hand = -dblquad(lambda x2, x1: 1/Fpar(x1, x2), 0, 1, 0, lambda x1: 1 - x1, epsabs=1e-13)[0]
print("by hand:  ", by_hand)
print("feynsage: ", C0(-1, -2, -3, 1, 2, 3).n(digits=16))
""", {}),

("d0", """
d = D0(-1, -2, -3, -4, -5, -6, 1, 2, 3, 4)
print(d.n(digits=25))
print(d0_numeric(-1, -2, -3, -4, -5, -6, 1, 2, 3, 4))
""", {}),

("cir", """
pole_parts(C0(0, s, 0, 0, 0, 0))
""", {}),

("irgraph", """
tri = graph("A-B, B-C, C-A", {"A": "p1 + p2", "B": "-p1", "C": "-p2"},
            kin={"p1^2": 0, "p2^2": 0, "p1.p2": "Q2/2"}, euclidean=True)
print("U =", tri.U(), "     F =", tri.F())
ep = var('eps')
dirichlet = gamma(1 + ep)*gamma(-ep)^2/gamma(1 - 2*ep)
(exp(ep*euler_gamma)*dirichlet).series(ep, 1).truncate()
""", {}),

("gramc", """
print("C1  =", PVC(0, 1, 0, m^2, 0, m^2, 0, m, m))
print("C00 =", PVC(1, 0, 0, m^2, 0, m^2, 0, m, m).log_expand().expand())
""", {}),

# ------------------------------------------------------------------ 3. graphs
("graphbub", """
bub = graph("A-B:ma, A-B:mb", {"A": "p", "B": "-p"}, kin={"p^2": "pp"}, euclidean=True)
print("lines N =", bub.N, "  vertices V =", bub.V, "  loops L = N - V + 1 =", bub.L)
fig = bub.plot(figsize=(2.6, 2.6))
""", {"fig": "fig"}),

("kiteplot", """
from feynsage.plotting import draw_panels, draw_sectors, draw_graph
kite = diagram("kite")
print(kite.lines)
fig = kite.plot(figsize=(2.8, 2.8))
""", {"fig": "fig"}),

("kitetrees", """
xs = lambda rem: "$" + "".join("x_{%d}" % (i + 1) for i in rem) + "$"
trees = kite.spanning_trees()
print(len(trees), "spanning trees; removed lines (from 0):", trees)
fig = draw_panels(kite, trees, titles=[xs(r) for r in trees], momenta=False, ncols=8, size=1.5)
""", {"fig": "fig"}),

("kiteU", """
kite.U()
""", {}),

("kiteforests", """
def entering(comp):
    tot = {}
    for v in comp:
        for k, c in kite.ext.get(v, {}).items():
            tot[k] = tot.get(k, 0) + c
    return {k: c for k, c in tot.items() if c}

forests = kite.two_forests()
carry = [(r, c) for r, c in forests if entering(c)]
print(len(forests), "2-forests,", len(carry), "of them carry the momentum p")
fig = draw_panels(kite, [r for r, c in carry], titles=[xs(r) + "$\\\\,p^2$" for r, c in carry],
                  momenta=False, ncols=8, size=1.5)
""", {"fig": "fig"}),

("kiteF", """
kite.F()
""", {}),

("threeways", """
import time
def three_ways(g, name):
    t0 = time.time()
    U, F = g.U(), g.F()
    t1 = time.time()
    props, loops = g.family()
    Um, Fm = IntegralFamily(name, loops, g.kin, props).UF()
    same = g.U_kirchhoff() == U and g.R(str(Um)) == U and g.R(str(Fm)) == F
    print("%-10s L = %d, %2d lines, %3d trees, %3d 2-forests (%.3f s); trees = Kirchhoff = det M: %s"
          % (name, g.L, g.N, len(g.spanning_trees()), len(g.two_forests()), t1 - t0, same))

for name in ["bubble_mass", "kite", "vertex2", "banana3", "ladder3", "triplebox"]:
    three_ways(diagram(name), name)
""", {}),

("kirchhoff", """
X = [SR.var('x%d' % (i + 1)) for i in range(kite.N)]
idx = {v: n for n, v in enumerate(kite.vertices)}
Lap = matrix(SR, kite.V, kite.V)
for i, (u, v, _) in enumerate(kite.lines):     # line i between u and v, weight 1/x_i
    a, b = idx[u], idx[v]
    Lap[a, a] += 1/X[i]; Lap[b, b] += 1/X[i]; Lap[a, b] -= 1/X[i]; Lap[b, a] -= 1/X[i]
print("vertices in this order:", kite.vertices)
Lap
""", {}),

("kirchhoff2", """
reduced = Lap[1:, 1:]                           # strike out one row and column
U_from_det = (prod(X) * reduced.det()).expand()
print("prod(x) * det(reduced Laplacian) =", U_from_det)
print("equal to U from the trees:", bool(U_from_det == SR(str(kite.U())).expand()))
""", {}),

("vertex2", """
v2 = diagram("vertex2")
fig = v2.plot(figsize=(3.0, 3.0))
""", {"fig": "fig"}),

("vertex2UF", """
R = v2.R
x1, x2, x3, x4, x5, x6 = R.gens()
q2 = R.base_ring()('q2')
U_note = (x2 + x3)*(x1 + x4 + x5) + (x1 + x2 + x3 + x4 + x5)*x6
F_note = -q2*(x1*x3*x4 + x1*x2*(x3 + x4) + x2*x3*(x4 + x5) + (x1 + x2)*(x3 + x4)*x6)
print("U as in the note:", v2.U() == U_note, "   F as in the note:", v2.F() == F_note)
""", {}),

# ------------------------------------------------------------------ 4. sectors
("kitefam", """
props, loops = kite.family()
kfam = IntegralFamily("kite", loops, kite.kin, props)
for i, (q, m2) in enumerate(props):
    print("line %d (x%d):  D%d = (%s)^2" % (i + 1, i + 1, i + 1,
          " + ".join(("" if c == 1 else "-" if c == -1 else str(c)) + k for k, c in q.items()).replace("+ -", "- ")))
""", {}),

("sectorof", """
def sector(a):
    return tuple(1 if x > 0 else 0 for x in a)

for a in [(1, 1, 1, 1, 1), (2, 1, 1, 1, 3), (1, 1, 1, 1, -2), (1, 1, 1, 1, 0), (0, 1, 1, 0, 1), (0, 2, 1, -1, 1)]:
    print("F%s  is in sector %s" % (a, sector(a)))
""", {}),

("contract", """
import matplotlib.pyplot as plt
h, kept = kite.contract([5])           # the kite without line 5
print("lines kept:", kept, "   vertices:", h.vertices)
fig, axes = plt.subplots(1, 2, figsize=(5.2, 2.5))
draw_graph(kite, ax=axes[0], momenta=False, removed=[5], title="line 5 removed (dashed)")
draw_graph(h, ax=axes[1], names=kept, momenta=False, title="line 5 shrunk to a point")
""", {"fig": "fig"}),

("contractUF", """
to_kite = h.R.hom([kite.R.gen(k - 1) for k in kept], kite.R)    # x_k of h -> x_kept[k] of the kite
x5 = kite.R.gen(4)
print("U(kite) at x5 = 0  =  U(kite with line 5 shrunk):", kite.U().subs({x5: 0}) == to_kite(h.U()))
print("F(kite) at x5 = 0  =  F(kite with line 5 shrunk):", kite.F().subs({x5: 0}) == to_kite(h.F()))
print("the same graph with line 5 deleted instead has", len(kite._graph([0, 1, 2, 3]).connected_components()),
      "component(s) and only", kite.L - 1, "loop: not the same integral")
""", {}),

("allsectors", """
from itertools import product as cartesian
secs = [s for s in cartesian([1, 0], repeat=5) if sum(s) in (4, 3)]
titles = [str(s).replace(" ", "") + ("  zero" if kfam.is_zero_sector(s) else "") for s in secs]
print(sum(kfam.is_zero_sector(s) for s in secs), "of", len(secs), "are zero")
fig = draw_sectors(kite, secs, titles=titles, momenta=False, ncols=5, size=2.0)
""", {"fig": "fig"}),

("zerosector", """
s0 = (1, 1, 1, 0, 0)
fig = draw_sectors(kite, [s0, (1, 0, 1, 0, 1)], titles=[str(s0), "(1,0,1,0,1)"], momenta=True, ncols=2, size=2.4)
U0, F0 = kfam.UF(s0)
print("U =", U0, "    F =", F0)
print("zero sector (Lee's criterion):", kfam.is_zero_sector(s0))
""", {"fig": "fig"}),

("leecrit", """
for s in [(1, 1, 1, 0, 0), (0, 1, 1, 0, 1), (1, 1, 1, 1, 0), (1, 1, 1, 1, 1)]:
    U0, F0 = kfam.UF(s)
    print(s, " U + F =", U0 + F0, "   zero:", kfam.is_zero_sector(s))
""", {}),

("count", """
from collections import Counter
allsec = [s for s in cartesian([0, 1], repeat=5) if sum(s) > 0]
nz = Counter(sum(s) for s in allsec if not kfam.is_zero_sector(s))
tot = Counter(sum(s) for s in allsec)
for k in sorted(tot):
    print("%d lines: %2d sectors, %2d non-zero" % (k, tot[k], nz.get(k, 0)))
""", {}),

("kitemasters", """
red = ibp_reduce(kfam, ["F(1,1,1,1,1)"])
print("masters:", red.masters)
fig = draw_sectors(kite, [tuple(1 if x > 0 else 0 for x in m) for m in red.masters],
                   titles=["master " + str(m).replace(" ", "") for m in red.masters], momenta=False, ncols=2, size=2.4)
""", {"fig": "fig"}),

("kitered", """
red
""", {}),

("sectorsym", """
from feynsage.laporta import Reducer
r = Reducer(kfam, symmetries=symmetries(kfam))
for a in [(1, 0, 0, 1, 1), (2, 0, 0, 1, 1), (1, 0, 0, 1, 3), (0, 1, 1, 0, 1)]:
    print("F%s  ->  F%s" % (a, r.canon(a)))
fig = draw_sectors(kite, [(0, 1, 1, 0, 1), (1, 0, 0, 1, 1)], titles=["(0,1,1,0,1)", "(1,0,0,1,1)"],
                   momenta=True, ncols=2, size=2.4)
""", {"fig": "fig"}),

("pak", """
for sec in [(0, 1, 1, 0, 1), (1, 0, 0, 1, 1)]:
    U0, F0 = kfam.UF(sec)
    print(sec, "  U + F =", U0 + F0)
""", {}),

("sunset63", """
from feynsage.ff import reduce_ff
sun = family([("k", "1"), ("k - l", "1"), ("l + p", "1"), "k + p", "l"], kin={"p^2": "pp"})
targets = [(2, 1, 1, 0, 0), (2, 2, 1, 0, 0), (1, 1, 2, 0, 0)]
for ss in [False, True]:
    r = Reducer(sun, symmetries=symmetries(sun), sector_symmetries=ss)
    tab = reduce_ff(r, targets, rmax=2, smax=1)
    ms = sorted({m for row in tab.values() for m in row})
    print("sector symmetries %-5s: %d masters %s" % (ss, len(ms), ms))
""", {}),

("numrel", """
r = Reducer(sun, symmetries=symmetries(sun), numerator_relations=True)
for rel in r.symmetry_relations((1, 1, 1, -1, 0)):
    print(" + ".join("(%s) F%s" % (c, a) for a, c in rel.items()).replace("+ (-", "- ("), "= 0")
""", {}),

# ------------------------------------------------------------------ 5. IBP and masters
("bubfam", """
bubf = family([("l", "1"), ("l - p", "1")], kin={"p^2": "pp"}, euclidean=True)
bubf.info()
""", {}),

("bubibp", """
for k, ident in enumerate(bubf.ibp((1, 1))):
    print("identity %d:" % (k + 1), " + ".join("(%s) F%s" % (c, a) for a, c in ident.items()), "= 0")
""", {}),

("bubred", """
rb = ibp_reduce(bubf, ["F(2,1)", "F(2,2)"])
rb
""", {}),

("bubcheck", """
row = rb[(2, 1)]
dd, ppp = row[(1, 1)].parent().gens()
print(row[(1, 1)] == -(dd - 3)/(ppp + 4), row[(0, 1)] == -(dd - 2)/(2*(ppp + 4)))
""", {}),

("kitevalue", """
from feynsage.oneloop import G, D, expand_eps
kq = family(["l1", "l1 + q", "l1 + l2", "l1 + l2 + q", "l2"], kin={"q^2": 1}, euclidean=True)
rk = ibp_reduce(kq, ["F(1,1,1,1,1)"])
value = {(0, 1, 1, 0, 1): G(1, 1)*G(1, 2 - D/2), (1, 1, 1, 1, 0): G(1, 1)^2}
kite_value = sum(SR(str(c)).subs({SR.var('d'): D}) * value[mm] for mm, c in rk[(1, 1, 1, 1, 1)].items())
expand_eps(kite_value, 0, loops=2)
""", {}),

("methods", """
res = {meth: ibp_reduce(kq, ["F(2,2,1,2,2)"], method=meth) for meth in ["ff", "trimmed", "exact"]}
for meth, r in res.items():
    print("%-8s %.2f s" % (meth, r.seconds))
print("the same coefficients:", res["ff"].table == res["trimmed"].table == res["exact"].table)
""", {}),

("ffverbose", """
r = ibp_reduce(kq, ["F(2,2,1,2,2)"], method="ff", verbose=True)
""", {}),

("lp", """
def critical_points(fam, sector, numbers):
    U, F = fam.UF()
    n = len(sector)
    R = PolynomialRing(QQ, ['x%d' % (i + 1) for i in range(n)] + ['t'])
    x, t = R.gens()[:n], R.gens()[n]
    G = R(str(SR(str(U + F)).subs(numbers)))
    live = [x[i] for i in range(n) if sector[i]]
    I = R.ideal([G.derivative(v) for v in live] + [x[i] for i in range(n) if not sector[i]] + [t*prod(live) - 1])
    return 0 if I.is_one() else I.vector_space_dimension()

pp = var('pp')
for sec in [(1, 1), (1, 0), (0, 1)]:
    print(sec, critical_points(bubf, sec, {pp: 7/3}))
""", {}),

# ------------------------------------------------------------------ 6. FORM
("tr1", """
form.trace(['mu', 'nu', 'rho', 'sigma'])
""", {}),

("tr2", """
print(form.dirac_trace(['mu', 'mu'], vectors=[], dim='D'))
print(form.dirac_trace(['g5', 'mu', 'nu', 'rho', 'sigma'], vectors=[]))
print(form.dirac_trace(['p + m', 'mu', 'k + m', 'nu'], vectors=['p', 'k']))
""", {}),

("ee", """
code = '''
Symbols s, t, u;
Vectors p, pp, k, kp;
Indices mu, nu;
Local M = g_(1, pp, mu, p, nu) * g_(2, k, mu, kp, nu);
trace4, 1;
trace4, 2;
id p.k = -t/2;   id pp.kp = -t/2;
id p.kp = -u/2;  id pp.k = -u/2;
Print;
.end
'''
traces = SR(" ".join(form.run_form(code).split("M =")[1].split(";")[0].split()))
traces
""", {}),

("eexs", """
s, t, u, e, th, alpha = var('s t u e theta alpha')
msq = (e^4/s^2 * traces/4).subs(t=-s/2*(1 - cos(th)), u=-s/2*(1 + cos(th)))
dsig = (msq/(64*pi^2*s)).subs(e=sqrt(4*pi*alpha))
print("dsigma/dOmega =", dsig.trig_simplify().factor())
print("sigma =", integrate(integrate(dsig*sin(th), th, 0, pi), var('phi'), 0, 2*pi))
""", {}),

]
