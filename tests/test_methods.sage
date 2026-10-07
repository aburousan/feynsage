# Tests of the method modules added for the lecture note: expansions, de, counting, sector, gpl,
# parametric, calculus.  Every check has a known answer from the note (checked there with Mathematica).
#     sage tests/test_methods.sage
from feynsage import *
from feynsage import gpl
ok = []
def check(name, cond):
    ok.append(bool(cond)); print(("PASS " if cond else "FAIL ") + name)

# expansions: Gamma(-eps) and psi at a pole
L = laurent(gamma(-eps)*exp(-eps*euler_gamma), eps, 1)
check("laurent Gamma(-eps)", (L - (-1/eps - pi^2/12*eps)).simplify_full() == 0)
check("laurent psi(-1+eps)", (laurent(psi(-1 + eps), eps, 0) - (-1/eps + 1 - euler_gamma)).simplify_full() == 0)
z = var('z')
check("residue Gamma(z) at -2", residue(gamma(z), z, -2) == 1/2)
check("residue Gamma(-1-eps-z) at -1-eps", (residue(gamma(-1 - eps - z), z, -1 - eps) + 1).simplify_full() == 0)

# de: the massless bubble J(1,1) ~ pp^(d/2 - 2), so dJ/dpp = (d - 4)/(2 pp) J; the equal-mass bubble has
#     two masters J(1,1), J(0,1) and dJ(0,1)/dpp = 0
J = family(["l", "l - p"], kin={"p^2": "pp"}, euclidean=True, name="J")
A = differential_equation(J, ["J(1,1)"], "pp")
d, pp = var('d pp')
check("massless bubble DE", (A[0, 0] - (d - 4)/(2*pp)).simplify_full() == 0)
Jm = family([["l", "m2"], ["l - p", "m2"]], kin={"p^2": "pp"}, euclidean=True, name="Jm")
Am = differential_equation(Jm, ["J(1,1)", "J(0,1)"], "pp")
check("equal-mass bubble DE", Am.nrows() == 2 and Am[1, 0] == 0 and Am[1, 1] == 0)

# counting: the equal-mass bubble sector has one master (critical point of G = U + F), the family has
#     two (the tadpole sector adds one); the massless box family has three
x1, x2 = var('x1 x2')
check("critical points bubble sector", critical_points(x1 + x2 + (x1 + x2)^2 + x1*x2, [x1, x2]) == 1)
eq = family([["l", "m^2"], ["l - p", "m^2"]], kin={"p^2": "pp"}, euclidean=True)
check("master_count equal-mass bubble", sum(master_count(eq, [(1, 1), (1, 0)], {"pp": 11/5, "m": 1}).values()) == 2)

# sector decomposition: massless on-shell triangle, s = -1: -1/eps^2 + pi^2/12 (with -Gamma(1+eps) e^(eps gE))
s = var('s')
g = graph("A-B, B-C, C-A", {"A": "p1 + p2", "B": "-p1", "C": "-p2"}, kin={"p1^2": 0, "p2^2": 0, "p1.p2": "s/2"})
xs = [SR.var('x%d' % i) for i in (1, 2, 3)]
secs = sector_decompose([(SR(g.U()), -1 + 2*eps), (SR(g.F()), -1 - eps)], xs)
T = (laurent(-gamma(1 + eps)*exp(eps*euler_gamma), eps, 2)*integrate_sectors(secs, eps, 0, {s: -1})).expand()
check("sector triangle", abs(T.coefficient(eps, -2) + 1) < 1e-10 and abs(T.coefficient(eps, -1)) < 1e-10
      and abs(T.coefficient(eps, 0) - (pi^2/12).n()) < 1e-10)

# gpl: G(1,0;1) = pi^2/6, G(-1,0,0;1) = 3 zeta(3)/4, and one order of a 1x1 canonical DE
check("gpl table", gpl.at_one(gpl.GPL.word(("1", "0"))) == pi^2/6
      and gpl.at_one(gpl.GPL.word(("-1", "0", "0"))) == 3*zeta(3)/4)
f1 = gpl.solve_canonical([gpl.GPL.const(1)], [(matrix([[2]]), 0)], [], [0])
check("gpl canonical order", str(f1[0].to_sr(lambda w, v: gpl.to_log(w, v) if w else 1)) == "2*log(x)")

# parametric: the massless bubble, U = x1 + x2, F = pp x1 x2
pref, integrand, xv = feynman_parametrize(J)
check("feynman_parametrize bubble", integrand.has(pp) and len(xv) == 2)

# calculus
x = var('x')
check("delta_integrate", delta_integrate(x^2, x - 1/2, x, 0, 1) == 1/4)
check("principal_value", (principal_value(1/(x - 1/2), x, 0, 1, 1/2)).simplify_full() == 0)
check("integrate_termwise", (integrate_termwise(1/x - 1/(1 + x), x) - (log(x) - log(x + 1))).simplify_full() == 0)

# edge cases found in review
from feynsage.expansions import coefficient
check("laurent removable singularity in Gamma", (laurent(gamma((exp(eps) - 1)/eps - 1), eps, -1) - 2/eps).simplify_full() == 0)
check("coefficient at a Gamma pole", coefficient(gamma(-eps), eps, -1) == -1)
yv = SR.var('y')
check("gpl first leg at y = 0", str(gpl.solve_canonical([gpl.GPL.const(1)], [(matrix(SR, [[1]]), yv)], [], [0])[0]) == "(1)*G(0; x)")
check("critical points, empty sector", critical_points(1, []) == 1)
check("delta root on an end point", delta_integrate(1, x, x, 0, 1) == 0)
import tempfile, os
tf = tempfile.NamedTemporaryFile('w', suffix='.m', delete=False); tf.write("{F[1,1] -> F[1,0], F[2,1] -> 2*F[1,1]}"); tf.close()
kt = read_kira(tf.name); os.unlink(tf.name)
check("read_kira one-line list", len(kt) == 2 and str(kt[(2, 1)]) == "2*F(1, 1)")

# C0 where the Gram determinant vanishes (vertex at q^2 = 0, the g - 2 kinematics), large mass ratio:
# C0(m^2, m^2, 0; m, M, m) = -Int_0^1 dx (1 - x)/((1 - x)^2 m^2 + x M^2), exact by partial fractions
from feynsage.pv import C0, DiscB
a_ = var('a_'); assume(a_ > 1)
mm, MM = QQ(511)/1000, 60000
R0 = MM^2/mm^2; av = ((R0 - 2) + sqrt((R0 - 2)^2 - 4))/2
exact_c0 = (-integrate(((1 - x)/((x + a_)*(x + 1/a_))).partial_fraction(x), x, 0, 1).subs(a_=av)/mm^2).n(digits=150)
got = C0(mm^2, mm^2, 0, mm, MM, mm).n(digits=60)
check("C0 at vanishing Gram determinant, M/m = 1.2e5, 40 digits", abs((got.real() - exact_c0)/exact_c0) < 1e-40)
check("DiscB(0; m, m) = -2 exactly", DiscB(0, 1, 1).n(digits=60) == -2)

print("%d of %d passed" % (sum(ok), len(ok)))
