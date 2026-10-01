# Figures and numbers for the slides on the graph features of feynsage.
# Run from the repository root:  sage slides/graphs.sage
import sys, time, json; sys.path.insert(0, '.')
from itertools import product as cartesian
from feynsage import *
from feynsage.plotting import draw_graph, draw_panels, draw_sectors, STYLE
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib as mpl
O = 'slides/out/'
xs = lambda rem: "$" + "".join("x_{%d}" % (i + 1) for i in rem) + "$"

# 1. a gallery: every picture drawn by g.plot() / draw_graph from the graph alone
names = ["bubble_mass", "triangle", "box", "sunset", "kite", "phi4_two", "vertex2", "triplebox"]
with mpl.rc_context(STYLE):
    fig, axes = plt.subplots(2, 4, figsize=(9.6, 4.8))
    for ax, n in zip(axes.flat, names):
        g = diagram(n)
        draw_graph(g, ax=ax, momenta=False, title="%s   (L = %d)" % (n, g.L))
    fig.tight_layout()
fig.savefig(O + 'gallery.svg', bbox_inches='tight')

# 2. the two-loop vertex of the lecture note: spanning trees and the 2-forests that carry q^2
v2 = diagram("vertex2")
trees = v2.spanning_trees()
fig = draw_panels(v2, trees, titles=[xs(r) for r in trees], momenta=False, ncols=6, size=2.1)
fig.savefig(O + 'v2_trees.svg', bbox_inches='tight')
carry = []
for rem, comp in v2.two_forests():
    P = {}
    for vtx in comp:
        for k, c in v2.ext.get(vtx, {}).items():
            P[k] = P.get(k, 0) + c
    P = {k: c for k, c in P.items() if c}
    if P and v2.kin.square_external(P) != 0:
        carry.append(rem)
fig = draw_panels(v2, carry, titles=[xs(r) + r"$\,q^2$" for r in carry], momenta=False, ncols=5, size=2.1)
fig.savefig(O + 'v2_forests.svg', bbox_inches='tight')
json.dump({"trees": len(trees), "forests": len(v2.two_forests()), "carry": len(carry),
           "U": latex(SR(str(v2.U()))), "F": latex(SR(str(v2.F())).factor())}, open(O + 'v2.json', 'w'))
print("vertex2:", len(trees), "trees,", len(v2.two_forests()), "2-forests,", len(carry), "carry q^2")

# 3. sectors of the kite: lines with index 0 shrunk to points, zero sectors from Lee's criterion
kite = diagram("kite")
props, loops = kite.family()
kfam = IntegralFamily('kite', loops, kite.kin, props)
secs = [s for s in cartesian([1, 0], repeat=5) if sum(s) in (4, 3)]
titles = [str(s).replace(" ", "") + ("  zero" if kfam.is_zero_sector(s) else "") for s in secs]
fig = draw_sectors(kite, secs, titles=titles, momenta=False, ncols=5, size=2.1)
fig.savefig(O + 'kite_sectors.svg', bbox_inches='tight')
print("kite sectors:", sum(kfam.is_zero_sector(s) for s in secs), "of", len(secs), "zero")

# 4. momentum routing (cycle basis) and symmetries (automorphisms of U + F)
def momtext(q):
    t = " + ".join(("" if c == 1 else "-" if c == -1 else str(c) + "*") + k for k, c in q.items())
    return t.replace("+ -", "- ")
open(O + 'kite_routing.txt', 'w').write("\n".join("D%d = (%s)^2" % (i + 1, momtext(q)) for i, (q, m2) in enumerate(props)))
open(O + 'kite_syms.txt', 'w').write(str(symmetries(kfam)))
print(open(O + 'kite_routing.txt').read()); print(symmetries(kfam))

# 5. how many masters? Lee-Pomeransky critical points of G = U + F, sector by sector, against IBP
def group_closure(gens, t):
    G = {tuple(range(t))}
    frontier = list(G)
    while frontier:
        new = []
        for a in frontier:
            for g in gens:
                b = tuple(g[a[i]] for i in range(t))
                if b not in G:
                    G.add(b); new.append(b)
        frontier = new
    return G

def critical_points(fam, sector, numbers):
    U, F = fam.UF()
    n = len(sector)
    R = PolynomialRing(QQ, ['x%d' % (i + 1) for i in range(n)] + ['t'])
    x, t = R.gens()[:n], R.gens()[n]
    G = R(str(SR(str(U + F)).subs(numbers)))
    live = [x[i] for i in range(n) if sector[i]]
    I = R.ideal([G.derivative(v) for v in live] + [x[i] for i in range(n) if not sector[i]] + [t * prod(live) - 1])
    return 0 if I.is_one() else I.vector_space_dimension()

def count_masters(fam, numbers):
    t = fam.t
    grp = group_closure(symmetries(fam), t)
    seen, total = set(), 0
    for s in cartesian([0, 1], repeat=t):
        if sum(s) == 0 or s in seen:
            continue
        orbit = {tuple(s[p.index(i)] for i in range(t)) for p in grp}
        seen |= orbit
        if fam.is_zero_sector(s):          # scaleless: no master (and G = U has a whole curve of critical points)
            continue
        c = critical_points(fam, s, numbers)
        assert c < Infinity, "non-isolated critical points in the non-zero sector %s" % (s,)
        total += c
    return total

# the sunset needs two numerator lines (k + p, l) for IBP; the count uses its three propagators
sun3 = [("k", "1"), ("k - l", "2"), ("l + p", "3")]
cases = [
    ("bubble, equal masses", family([("l", "1"), ("l - p", "1")], kin={"p^2": "pp"}, euclidean=True), None, ["F(2,1)", "F(2,2)"]),
    ("bubble, two masses", family([("l", "1"), ("l - p", "2")], kin={"p^2": "pp"}, euclidean=True), None, ["F(2,1)", "F(1,2)"]),
    ("massless box", family(["k", "k + p1", "k + p1 + p2", "k - p3"],
                            kin={"p1^2": 0, "p2^2": 0, "p3^2": 0, "p1.p2": "s/2", "p1.p3": "1/2", "p2.p3": "-(s+1)/2"}), None,
     ["F(2,1,1,1)", "F(1,2,1,1)"]),
    ("sunset, three masses", family(sun3, kin={"p^2": "pp"}, euclidean=True),
     family(sun3 + ["k + p", "l"], kin={"p^2": "pp"}, euclidean=True),
     ["F(2,1,1,0,0)", "F(1,2,1,0,0)", "F(1,1,2,0,0)", "F(2,2,1,0,0)", "F(1,1,1,-1,0)"]),
]
rows = []
for name, fam, ibpfam, targets in cases:
    nums = {SR.var(str(v)): QQ(7) / 3 + k for k, v in enumerate(fam.kin.R.gens()[1:])}
    t0 = time.time(); lp = count_masters(fam, nums); tlp = time.time() - t0
    nibp = len(ibp_reduce(ibpfam or fam, targets).masters)
    rows.append({"name": name, "lp": int(lp), "ibp": int(nibp), "t": tlp})
    print(name, "critical points:", lp, " IBP masters:", nibp, " (%.2f s)" % tlp)
json.dump(rows, open(O + 'lp_rows.json', 'w'), indent=1)
