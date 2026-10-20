# Loop integrals as graphs (feynsage.topology) against FeynCalc 10.2: the examples of its documentation
# (FCLoopIntegralToGraph, FCLoopFindTopologyMappings, FCLoopFindSubtopologies, FCLoopBasisFindCompletion,
# FCGraphCuttableQ) and its output for FCLoopGLILowerDimension / FCLoopGLIRaiseDimension (run on the laptop,
# printed in the comments), plus checks that do not need FeynCalc: closed forms, IBP, isomorphism.
#     sage tests/test_topology.sage
import time, itertools, os
from feynsage import *
from feynsage.easy import _LIBRARY
from feynsage.graph import FeynmanGraph
from sage.all import Graph
t_start = time.time()
FAIL = []
def check(name, ok):
    print(("ok   " if ok else "FAIL ") + name)
    if not ok:
        FAIL.append(name)
def raises(f, err):
    try:
        f()
    except err:
        return True
    return False
def sage_graph(g):
    G = Graph(multiedges=True, loops=True)
    for u, v, m in g.lines:
        G.add_edge(u, v, str(m))
    return G

# ---------------------------------------------------------------------------- integral_to_graph
f = family(["q - k1", "k1", "q - k2", "k2", ("k2 - k3", "mb"), ("k1 - k3", "mb")], kin={"q^2": "s"}, loops="k1 k2 k3")
fc = FeynmanGraph([(1, 3, 0), (1, 4, 0), (2, 3, 0), (2, 4, 0), (3, 4, "mb^2"), (3, 4, "mb^2")], {1: {"q": 1}, 2: {"q": -1}}, f.kin)
check("FCLoopIntegralToGraph documentation example: the same graph as FeynCalc's",
      sage_graph(integral_to_graph(f)).is_isomorphic(sage_graph(fc), edge_labels=True))
g2 = integral_to_graph(family(["p1", "Q - p1", "p2", "Q - p2"], kin={"Q^2": "s"}, loops="p1 p2"))
check("a factorizing integral (two bubbles in a chain, where FeynCalc aborts) is a graph with a cut vertex",
      len(g2.vertices) == 3 and len(g2.ext) == 2)
ok = True
for name in _LIBRARY:
    G0 = diagram(name)
    props, loops = G0.family()
    G1 = integral_to_graph(IntegralFamily(name, loops, G0.kin, props))      # U and F are compared inside
    ok = ok and sage_graph(G0).is_isomorphic(sage_graph(G1), edge_labels=True)
check("graph -> momenta -> graph gives back an isomorphic graph for all %d ready-made diagrams" % len(_LIBRARY), ok)
check("a tree line (no loop momentum) is refused",
      raises(lambda: integral_to_graph(family(["l", "l + p", "p"], kin={"p^2": "s"}, loops="l")), NotAGraphError))

# ---------------------------------------------------------------------------- topology mappings
k3 = {"Q^2": "s"}
T1 = [family(p, kin=k3, loops="p1 p2 p3", name=n) for n, p in [
    ("fctopology1", ["p3", "p2", "p1", "p2 + p3", "p2 - Q", "p1 - Q", "p2 + p3 - Q", "p1 + p3 - Q", "p1 + p2 + p3 - Q"]),
    ("fctopology2", ["p3", "p2", "p1", "p2 + p3", "p2 - Q", "p1 - Q", "p2 + p3 - Q", "p1 + p2 - Q", "p1 + p2 + p3 - Q"]),
    ("fctopology3", ["p3", "p2", "p1", "p2 + p3", "p1 + p3", "p2 - Q", "p2 + p3 - Q", "p1 + p3 - Q", "p1 + p2 + p3 - Q"]),
    ("fctopology4", ["p3", "p2", "p1", "p2 + p3", "p1 + p3", "p2 - Q", "p1 - Q", "p1 + p3 - Q", "p1 + p2 + p3 - Q"]),
    ("fctopology5", ["p3", "p2", "p1", "p1 + p3", "p2 - Q", "p1 - Q", "p1 + p3 - Q", "p1 + p2 - Q", "p1 + p2 + p3 - Q"])]]
m1, ind1 = topology_mappings(T1)
check("FCLoopFindTopologyMappings, 3-loop example: 3 mappings, 2 independent (fctopology1, fctopology2)",
      len(m1) == 3 and [t.name for t in ind1] == ["fctopology1", "fctopology2"] and all(x["shift"] for x in m1))
k2 = {"p^2": "pp"}
T2 = [family(p, kin=k2, loops="q1 q2", name=n) for n, p in [
    ("fctopology1", ["q2", "q1", "q1 + q2", "p + q1", "p - q2"]), ("fctopology2", ["q2", "q1", "p + q2", "p - q1"]),
    ("fctopology3", ["q2", "q1", "p - q1", "p - q1 + q2"])]]
P2 = [family(p, kin=k2, loops="q1 q2", name=n) for n, p in [
    ("prop2L", ["q1", "q2", "q1 - q2", "q1 - p", "q2 - p"]), ("prop2LX1", ["q2", "q1 - q2", "q1 - p", "q2 - p"]),
    ("prop2LX3", ["q1", "q2", "q1 - p", "q2 - p"]), ("prop2LX15", ["q2", "q1 - q2", "q1 - p"])]]
m2, ind2 = topology_mappings(T2, preferred=P2)
want = {"fctopology1": ("prop2L", {"q1": {"q2": -1}, "q2": {"q1": 1}}, {0: 0, 1: 1, 2: 2, 3: 4, 4: 3}),
        "fctopology2": ("prop2LX3", {"q1": {"q2": 1}, "q2": {"q1": -1}}, {0: 0, 1: 1, 2: 2, 3: 3}),
        "fctopology3": ("prop2LX1", {"q1": {"q2": 1}, "q2": {"q1": -1, "q2": 1}}, {0: 1, 1: 0, 2: 3, 3: 2})}
check("FCLoopFindTopologyMappings with PreferredTopologies: the same targets, line maps and shifts as FeynCalc",
      all((x["target"].name, x["shift"], x["lines"]) == want[x["family"].name] for x in m2)
      and [t.name for t in ind2] == ["prop2L", "prop2LX1", "prop2LX3"])
r = map_integral(T2[0], (1, 2, 1, 1, 3), P2)
check("map_integral: fctopology1(1,2,1,1,3) = prop2L(1,2,1,3,1)", r is not None and r[0].name == "prop2L" and r[1] == (1, 2, 1, 3, 1))

# ---------------------------------------------------------------------------- subtopologies
tri = family(["p1", "p2", "p1 + Q1", "p1 + p2 + Q1", "-p1 + Q2", "-p1 - p2 + Q2"], loops="p1 p2",
             kin={"Q1^2": "a", "Q2^2": "b", "Q1.Q2": "c"})
check("FCLoopFindSubtopologies: 19 for the 2-loop example, 36 for the 3-loop one", len(subtopologies(tri)) == 19
      and len(subtopologies(T1[0])) == 36)
check("scaleless: a massless tadpole and a bubble without its external momentum",
      scaleless(family(["l"])) and scaleless(family(["l", "l"])) and not scaleless(family([("l", "m")])))

# ---------------------------------------------------------------------------- basis and partial fractions
f1 = family(["q", ("q - p + l", "m")], loops="q", kin={"p^2": "pp", "l^2": "ll", "p.l": "pl"})
f2 = family([("q1", "m1"), ("q2", "m2")], loops="q1 q2")
check("FCLoopBasisFindCompletion examples: one propagator added to each, now a basis",
      basis_status(f1) == (2, 3, 2) and basis_status(complete_basis(f1)) == (3, 3, 3) and complete_basis(f1).complete
      and basis_status(f2) == (2, 3, 2) and complete_basis(f2).complete)
f3 = family(["l", "l + p", "l - p"], kin={"p^2": "pp"})
K3 = f3.kin.K; pp = K3(f3.kin.R.gens()[1])
check("apart: 1/(l^2 (l+p)^2 (l-p)^2) = [J(1,0,1) + J(1,1,0) - 2 J(0,1,1)]/(2 p^2) (checked at random points inside)",
      apart(f3, (1, 1, 1)) == {(1, 0, 1): 1/(2*pp), (1, 1, 0): 1/(2*pp), (0, 1, 1): -1/pp})
f4 = family([("l", "m"), ("l + p", "m"), ("l + 2*p", "m"), ("l - p", "m")], kin={"p^2": "pp"})
r4 = apart(f4, (2, 1, 1, 1))
check("apart: four massive lines of one loop momentum become pairs (checked at random points inside)",
      all(sum(1 for x in k if x > 0) <= 2 for k in r4))

# ---------------------------------------------------------------------------- dimension shifts
kite = family(["p1", "p2", "Q - p1 - p2", "Q - p2", "Q - p1"], loops="p1 p2", kin={"Q^2": "qq"})
fc_lo = [(1, 1, 1, 2, 2), (1, 1, 2, 1, 2), (1, 1, 2, 2, 1), (1, 2, 1, 1, 2), (1, 2, 2, 1, 1), (2, 1, 1, 2, 1), (2, 1, 2, 1, 1), (2, 2, 1, 1, 1)]
check("FCLoopGLILowerDimension, kite: the same 8 integrals as FeynCalc", dimension_lower(kite, (1,) * 5) == {k: 1 for k in fc_lo})
K = kite.kin.K; d, qq = [K(g) for g in kite.kin.R.gens()[:2]]
num = {(-1, 0, 1, 1, 1): -1, (0, -1, 1, 1, 1): -1, (0, 0, 0, 1, 1): 1, (0, 0, 1, 0, 1): 1, (0, 0, 1, 1, 0): 1, (0, 0, 1, 1, 1): qq,
       (0, 1, 0, 0, 1): -1, (0, 1, 0, 1, 1): qq, (0, 1, 1, 0, 0): 1, (0, 1, 1, 1, 0): -qq, (1, 0, 0, 1, 0): -1, (1, 0, 0, 1, 1): qq,
       (1, 0, 1, 0, 0): 1, (1, 0, 1, 0, 1): -qq, (1, 1, -1, 1, 1): -qq, (1, 1, 0, 0, 0): 1, (1, 1, 0, 0, 1): qq, (1, 1, 0, 1, 0): qq,
       (1, 1, 0, 1, 1): -qq^2, (1, 1, 1, -1, 0): -1, (1, 1, 1, 0, -1): -1, (1, 1, 1, 0, 0): qq}
check("FCLoopGLIRaiseDimension, kite: the same 22 integrals and coefficients as FeynCalc",
      dimension_raise(kite, (1,) * 5) == {k: v/((d - 2)*(d - 1)*qq) for k, v in num.items()})
# one loop: the massive tadpole in closed form, J(a, d) = c (-1)^a Gamma(a - d/2)/Gamma(a) (m^2)^(d/2 - a), c = I (Minkowski)
dd, m2 = var('dd m2')
ok = True
for eu in (True, False):
    T = lambda a, dim: (1 if eu else (-1)^a*I) * gamma(a - dim/2)/gamma(a) * m2^(dim/2 - a)
    tad = family([("l", "m")], euclidean=eu)
    ev = lambda rel: sum(SR(str(c)).subs(d=dd, m=sqrt(m2)) * T(k[0], dd) for k, c in rel.items())
    ok = ok and bool((T(1, dd - 2)/ev(dimension_lower(tad, (1,)))).simplify_full() == 1)
    ok = ok and bool((T(2, dd + 2)/ev(dimension_raise(tad, (2,)))).simplify_full() == 1)
check("one loop: lowering and raising reproduce the closed-form massive tadpole (Euclidean and Minkowski)", ok)
# FeynCalc gives FCLoopGLILowerDimension[bub(1,1)] = J(1,2) + J(2,1) and FCLoopGLIRaiseDimension[tad(2)] = 2/D (J(1) + mm J(2)):
# the Euclidean relations; for Minkowski propagators 1/(p^2 - m^2) the closed form above needs (-1)^L, as feynsage has
bub = family([("l", "m"), ("l - Q", "m")], kin={"Q^2": "qq"})
back = {}
for k, c in dimension_lower(bub, (1, 1)).items():
    for k2, c2 in dimension_raise(bub, k).items():
        back[k2] = back.get(k2, 0) + c * c2
red = ibp_reduce(bub, list(back), method="exact")
Kb = bub.kin.K
tot = {}
for k, c in back.items():
    for mst, cm in red.table[k].items():
        tot[mst] = tot.get(mst, Kb(0)) + Kb(c) * Kb(cm)
check("raise(lower(J)) = J after IBP reduction, massive bubble", {m: v for m, v in tot.items() if v != 0} == {(1, 1): 1})

# ---------------------------------------------------------------------------- cuts
Kc = kinematics(["q1"], {"q1^2": "s"}, extra=[SR.var("m1"), SR.var("m2")])
g1 = FeynmanGraph([(1, 3, "m1^2"), (1, 4, "m1^2"), (2, 3, "m1^2"), (2, 4, "m2^2"), (2, 4, "m1^2"), (3, 4, 0)], {1: {"q1": 1}, 2: {"q1": -1}}, Kc)
g2c = FeynmanGraph([(1, 3, 0), (1, 4, "m1^2"), (2, 4, 0), (2, 5, "m1^2"), (3, 5, "m1^2"), (3, 6, "m1^2"), (4, 6, "m2^2"), (5, 6, "m1^2")],
                   {1: {"q1": 1}, 2: {"q1": -1}}, Kc)
check("FCGraphCuttableQ examples: graph1 False, graph2 True (avoiding m1)", not cuttable(g1, ["m1"]) and cuttable(g2c, ["m1"]))
bubg = integral_to_graph(family([("l", "m"), ("l - p", "m")], kin={"p^2": "s"}))
check("the equal-mass bubble has one cut, threshold s = 4 m^2", [(c[0], str(c[2])) for c in cuts(bubg)] == [([0, 1], "4*m^2")])

# ---------------------------------------------------------------------------- topologies of an amplitude
exec(open(os.path.join("tests", "feyncalc_ref", "gh_gh_2l_integrals.py")).read())
fc_tops = [family(p, kin={"p^2": "pp"}, loops="q1 q2") for p in (["q2", "q1", "q1 + q2", "p + q1", "p - q2"],
           ["q2", "q1", "p + q2", "p - q1"], ["q2", "q1", "p - q1", "p - q1 + q2"])]
t3, r3_ = find_topologies(ints, kin={"p^2": "pp"}, loops="q1 q2", merge_sectors=False, complete=False)
check("FCLoopFindTopologies, 2-loop ghost self-energy (%d integrals): the same 3 topologies as FeynCalc" % len(ints),
      sorted(str(canonical_form(t)[0]) for t in t3) == sorted(str(canonical_form(f)[0]) for f in fc_tops))
t1, r1_ = find_topologies(ints, kin={"p^2": "pp"}, loops="q1 q2")
check("... and with sectors merged (shifts of the loop momenta): one complete topology, every integral with its shift",
      len(t1) == 1 and t1[0].complete and all(r is not None and r[2] is not None for r in r1_))

# ---------------------------------------------------------------------------- tensor decomposition (Tdec)
q, p, q1, q2, q3, l, p1, p2, k1, k2 = momenta("q p q1 q2 q3 l p1 p2 k1 k2")
mu, nu, rho, sg = lorentz_indices("mu nu rho sg")
dd_ = spacetime_dimension()
check("Tdec documentation: q^mu -> p^mu (p.q)/p^2", bool((tdec([(q, mu)], [p]) - comp(p, mu)*dot(p, q)/dot(p, p)).simplify_full() == 0))
r2 = tdec([(q1, mu), (q2, nu), (q3, rho)], [p])
check("Tdec documentation: q1^mu q2^nu q3^rho with p, contracted with p p p = (p.q1)(p.q2)(p.q3)",
      bool((contract(r2*comp(p, mu)*comp(p, nu)*comp(p, rho)) - dot(p, q1)*dot(p, q2)*dot(p, q3)).simplify_full() == 0))
r3 = tdec([(l, mu), (l, nu)], [p1, p2], rules={dot(p1, p1): 0, dot(p2, p2): 0})
fc_g = -(dot(l, l)*dot(p1, p2)^2 - 2*dot(l, p1)*dot(l, p2)*dot(p1, p2))/((2 - dd_)*dot(p1, p2)^2)
fc_22 = (2*dot(l, p1)^2*dot(p1, p2)^2 - dd_*dot(l, p1)^2*dot(p1, p2)^2)/((2 - dd_)*dot(p1, p2)^4)
check("Tdec documentation: l^mu l^nu with p1^2 = p2^2 = 0, the g^{mu nu} and p2^mu p2^nu coefficients as FeynCalc's",
      bool((r3.coefficient(metric(mu, nu)) - fc_g).simplify_full() == 0)
      and bool((r3.coefficient(comp(p2, mu)*comp(p2, nu)) - fc_22).simplify_full() == 0))
check("Tdec: l^4 with two momenta (Pak symmetrisation: 9 tensors) and q1 q1 q2 at two loops (checked by FORM inside)",
      len(tdec([(l, mu), (l, nu), (l, rho), (l, sg)], [k1, k2]).operands()) == 9
      and len(tdec([(q1, mu), (q1, nu), (q2, rho)], [k1, k2], rules={"k1^2": 0, "k2^2": 0, "k1.k2": "s/2"}).operands()) == 10)

print("time %.0f s" % (time.time() - t_start))
print("ALL TOPOLOGY CHECKS PASS" if not FAIL else "FAILED: %s" % FAIL)
