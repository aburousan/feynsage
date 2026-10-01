# The code cells of the feynsage tutorial, in order.  They follow sir's two lectures "Feynman
# integral calculus" (NISER, 2026), part 1 and part 2.  Each is (id, code, options); the cells
# share one namespace, like a notebook.  Options: fig=<name> saves that matplotlib figure,
# text=True shows the value as plain text.
CELLS = [

# ------------------------------------------------------------------ the tadpole
("import", """
from feynsage import *
from feynsage.oneloop import tadpole, bubble_equal_mass, expand_eps, D
""", {}),

("solidangle", """
n = var('n')
Omega = 2*pi^(n/2)/gamma(n/2)            # from the Gaussian integral
[Omega.subs(n=k) for k in [2, 3, 4]]
""", {}),

("In", """
m2 = var('m2')
tadpole(n, m2)                           # Int d^Dk/pi^(D/2) 1/(k^2 + m^2)^n
""", {}),

("derivrel", """
for k in [2, 3, 4]:
    lhs = tadpole(k, m2)
    rhs = (-1)^(k - 1)/factorial(k - 1) * diff(tadpole(1, m2), m2, k - 1)
    print("n = %d:  I_n - (-1)^(n-1)/(n-1)! (d/dm^2)^(n-1) I_1 =" % k, (lhs - rhs).simplify_full())
""", {}),

("pole", """
expand_eps(tadpole(1, m2), 0)
""", {}),

("phi4", """
lam = var('lam', latex_name=r'\\lambda'); mu2 = var('mu2', latex_name=r'\\mu^2'); w, e = var('omega e')
# (1/2)(-lambda)(mu^2)^(2-w) Int d^(2w)l/(2 pi)^(2w) 1/(l^2 + m^2), the integral done as above
T = 1/2*(-lam)*mu2^(2 - w) * pi^w*gamma(1 - w)*m2^(w - 1)/(2*pi)^(2*w)
Texp = T.subs({w: 2 - e}).series(e, 1).truncate()                          # e = 2 - omega
lecture = lam*m2/(32*pi^2)*(1/e + psi(2) + log(4*pi*mu2/m2))
print("minus the lecture's result:", (Texp - lecture).canonicalize_radical().simplify_full())
Texp.collect(e)
""", {}),

("psi2", """
psi(2), psi(2).n()
""", {}),

# ------------------------------------------------------------------ IBP for the tadpole
("tadfam", """
tad = family([("k", "m")], euclidean=True, name="T")
tad.info()
""", {}),

("tadibp", """
for k in [1, 2, 3]:
    print("seed T(%d):" % k, tad.ibp((k,)))      # {(n,): c} means c T(n)
""", {}),

("tadred", """
r = ibp_reduce(tad, ["T(2)", "T(3)", "T(4)"])
r
""", {}),

("taddraw", """
fig = r.draw(graph("A-A:m", {}), size=1.0)      # the tadpole: one line from A back to A
""", {"fig": "fig"}),

("tadcheck", """
mm = var('m')
for k in [2, 3, 4]:
    c = SR(str(r[(k,)][(1,)])).subs({SR.var('d'): D})
    print("T(%d)/T(1) from IBP equals the Gamma-function ratio:" % k,
          bool((c - tadpole(k, mm^2)/tadpole(1, mm^2)).simplify_full() == 0))
""", {}),

("diffeq", """
Uformula = lambda n: gamma(n - D/2)/(gamma(1 - D/2)*gamma(n))      # times U(1)
(-(n - D/2)*Uformula(n) + n*Uformula(n + 1)).simplify_full()
""", {}),

("factorial", """
t = var('t')
v = function('v')(t)
sol = desolve(diff(v, t) == (D/2 - t)/(t*(t - 1))*v, v, ivar=t)
print("v(t) =", sol)
Dv, nv = 2.6, 3
num = numerical_integral(t^(nv - 1) * t^(-Dv/2)*(1 - t)^(Dv/2 - 1) / gamma(Dv/2), 0, 1)[0]
print("Int_0^1 t^(n-1) v(t) dt with v0 = 1/Gamma(D/2):", num)
print("U(n) = Gamma(n - D/2)/Gamma(n):                ", tadpole(nv, 1).subs(D=Dv).n())
""", {}),

# ------------------------------------------------------------------ the basis
("basis", """
print("A0(1)                     =", A0(1))
print("B0(5; 1, 1), finite part  =", finite_part(B0(5, 1, 1), 1).n(digits=12))
print("C0(1, 2, 3; 1, 2, 3)      =", C0(1, 2, 3, 1, 2, 3).n(digits=12))
print("D0(1,2,3,4; -5,-6; 1,...) =", D0(1, 2, 3, 4, -5, -6, 1, 1, 1, 1).n(digits=12))
""", {}),

# ------------------------------------------------------------------ the bubble
("bubfam", """
bub = family([("p", "m"), ("p - k", "m")], kin={"k^2": "kk"}, loops=["p"], euclidean=True, name="J")
bub.info()
""", {}),

("bubibp", """
for ident in bub.ibp((1, 1)):
    print(" + ".join("(%s) J(%s)" % (c, ",".join(map(str, a))) for a, c in ident.items()), "= 0")
""", {}),

("bubred", """
rb = ibp_reduce(bub, ["J(2,1)", "J(1,2)", "J(0,2)"])
rb
""", {}),

("bubcheck", """
d_, kk_, m_ = rb[(2, 1)][(1, 1)].parent().gens()
print("dotted bubble:  ", rb[(2, 1)][(1, 1)] == -(d_ - 3)/(kk_ + 4*m_^2),
      "  tadpole part:", rb[(2, 1)][(0, 1)] == 1/(kk_ + 4*m_^2) * (-(d_ - 2)/(2*m_^2)))
print("dotted tadpole: ", rb[(0, 2)][(0, 1)] == -(d_ - 2)/(2*m_^2))
""", {}),

("bubdraw", """
gb = graph("A-B:m, A-B:m", {"A": "k", "B": "-k"}, kin={"k^2": "kk"}, euclidean=True)
fig = rb.draw(gb, targets=[(2, 1), (0, 2)], rename={"kk": "k^2"}, size=1.0)
""", {"fig": "fig"}),

("bubde", """
# k.dJ/dk = J(0,2) - J(1,1) - k^2 J(1,2) and dJ/dk^2 = (k.dJ/dk)/(2 k^2); reduce the right side
kk = SR.var('kk'); dd = SR.var('d'); m = SR.var('m')
c = lambda a, mm_: SR(str(rb[a].get(mm_, 0)))
A = (c((0, 2), (1, 1)) - 1 - kk*c((1, 2), (1, 1)))/(2*kk)        # coefficient of J(1,1)
B = (c((0, 2), (0, 1)) - kk*c((1, 2), (0, 1)))/(2*kk)            # coefficient of T(1)
print("dJ/dk^2 = A J + B T(1) with")
print("  A =", A.factor())
print("  B =", B.factor())
print("as in the notes:", bool((A + (1/kk - (dd - 3)/(kk + 4*m^2))/2).simplify_full() == 0),
      bool((B + (dd - 2)/(4*m^2)*(1/kk - 1/(kk + 4*m^2))).simplify_full() == 0))
""", {}),

("dex", """
x = var('x')
# write k^2 = 4 m^2 x; then dJ/dx = 4 m^2 dJ/dk^2
Ax = (4*m^2*A).subs(kk=4*m^2*x).subs({dd: D}).simplify_full()
Bx = (4*m^2*B).subs(kk=4*m^2*x).subs({dd: D}).simplify_full()
print("dJ/dx = (%s) J + (%s) T(1)" % (Ax.partial_fraction(x), Bx.partial_fraction(x)))
""", {}),

("dehom", """
J0 = function('J0')(x)
hom = desolve(diff(J0, x) == Ax*J0, J0, ivar=x)       # the homogeneous equation
hom.canonicalize_radical()                             # write e^(a log u) as u^a
""", {}),

("deseries", """
N = 6
cs = [var('c%d' % i) for i in range(N + 1)]
T1 = var('T1')                                                  # the master T(1)
Js = sum(cs[i]*x^i for i in range(N + 1))                       # J as a power series in x
eq = (x*(1 + x)*(diff(Js, x) - Ax*Js - Bx*T1)).simplify_rational().expand()
eqs = [eq.coefficient(x, k) == 0 for k in range(N)]
print("the x^0 equation:", eqs[0])
sol = solve(eqs, cs[:N], solution_dict=True)[0]
is_T2 = bool((sol[cs[0]] + (D - 2)/(2*m^2)*T1).simplify_full() == 0)
print("so J(D, 0) = c0 =", sol[cs[0]].factor(), "   equal to T(2):", is_T2)
for n in range(4):
    ratio = (sol[cs[n + 1]]/sol[cs[n]]).factor()
    hyp = -(2 - D/2 + n)*(1 + n)/((3/2 + n)*(1 + n))   # term ratio of the 2F1
    print("c%d/c%d =" % (n + 1, n), ratio, "   same as 2F1:", bool((ratio - hyp).simplify_full() == 0))
""", {}),

("denum", """
from scipy.integrate import solve_ivp
Dv, mv = 3.3, 1
T1v = float(tadpole(1, mv^2).subs(D=Dv))
fA = fast_callable(Ax.subs(m=mv, D=Dv), vars=[x])
fB = fast_callable(Bx.subs(m=mv, D=Dv), vars=[x])
x0 = 1e-6                                                       # start next to the boundary x = 0
J_start = float(sum(sol[cs[i]].subs(D=Dv, m=mv, T1=T1v)*x0^i for i in range(N)))
num = solve_ivp(lambda t, y: [fA(t)*y[0] + fB(t)*T1v], [x0, 2.0], [J_start],
                rtol=1e-11, atol=1e-13, dense_output=True)
for xv in [0.125, 0.5, 2.0]:
    exact = bubble_equal_mass(4*mv^2*xv, mv^2).subs(D=Dv).n()     # Gamma(2-D/2) m^(D-4) 2F1
    print("x = %.3f   solved numerically: %.10f    closed form: %.10f" % (xv, num.sol(xv)[0], exact))
""", {}),

("hyp", """
a, b, c, z = 7/20, 1, 3/2, -1/2           # 2F1(2 - D/2, 1; 3/2; -x) at D = 3.3, x = 1/2
h = hypergeometric([a, b], [c], z)
print("Sage:          ", h.n(digits=25))
term = lambda n: rising_factorial(a, n)*rising_factorial(b, n)/rising_factorial(c, n)*z^n/factorial(n)
print("the series:    ", sum(term(n) for n in range(80)).n(digits=25))
import mpmath; mpmath.mp.dps = 25
print("mpmath.hyp2f1: ", mpmath.hyp2f1(a, b, c, z))
assume(x > 0)
h3 = hypergeometric([1/2, 1], [3/2], -x)                # D = 3
print("D = 3 in closed form:", h3.simplify_hypergeometric(algorithm='maxima'))
""", {}),

# ------------------------------------------------------------------ the fish and Feynman parameters
("feyntrick", """
A_, B_, x = var('A B x')
assume(A_ > 0, B_ > 0, B_ - A_ > 0)          # any order of A and B gives the same
integrate(1/(x*A_ + (1 - x)*B_)^2, x, 0, 1).factor()
""", {}),

("fishfp", """
Dv, mv, pv = 3.3, 1, 3
# after the shift l -> l + p(1 - x) the loop integral is a tadpole with mass^2 = m^2 + p^2 x(1 - x)
fp = numerical_integral(lambda xx: tadpole(2, mv^2 + pv*xx*(1 - xx)).subs(D=Dv).n(), 0, 1)[0]
print("Feynman parameter integral:", fp)
print("bubble from the DE (2F1):  ", bubble_equal_mass(pv, mv^2).subs(D=Dv).n())
""", {}),

("ramond", """
p2v = 3
finite = -numerical_integral(lambda xx: log(1 + xx*(1 - xx)*p2v), 0, 1)[0]
b = sqrt(1 + 4/p2v)
print("finite part, integral over x:    ", finite)
print("Ramond: 2 - b ln((b + 1)/(b - 1)):", (2 - b*log((b + 1)/(b - 1))).n())
""", {}),

("continuation", """
print("ln(-2) =", log(-2).n(), "   ln 2 + i pi =", (log(2) + I*pi).n())
for sv in [3, 5, 10]:
    val = finite_part(B0(sv, 1, 1), 1).n(digits=12)
    beta = sqrt(1 - 4/sv) if sv > 4 else 0
    print("s = %2d (m = 1):  Im of the bubble = %.10f    pi*beta = %.10f" % (sv, imag(val), (pi*beta).n()))
""", {}),

("b0plot", """
from feynsage.plotting import quick_plot, set_theme
set_theme()
s = var('s')
fig = quick_plot([B0(s, 1, 1)], (s, -2, 14), labels=["bubble, $m = 1$"])
""", {"fig": "fig"}),

# ------------------------------------------------------------------ part 2: graph polynomials
("bub2", """
from feynsage.plotting import draw_panels
xs = lambda rem: "$" + "".join("\\\\alpha_{%d}" % (i + 1) for i in rem) + "$"
g2 = graph("A-B:m1, A-B:m2", {"A": "p", "B": "-p"}, kin={"p^2": "pp"}, euclidean=True)
trees = g2.spanning_trees(); forests = [r for r, c in g2.two_forests()]
fig = draw_panels(g2, trees + forests, titles=["1-tree: " + xs(r) for r in trees] + ["2-tree: " + xs(r) for r in forests],
                  momenta=False, ncols=3, size=1.9)
print("U =", g2.U())
print("V =", g2.F0())
print("F = V + (m1^2 x1 + m2^2 x2) U =", g2.F())
""", {"fig": "fig"}),

("box", """
# Smirnov Fig. 3.6 as drawn in the lecture: line 1 on top, 2 on the left, 3 on the right, 4 at the
# bottom; p1 and p2 enter on the left, p3 and p4 on the right, s = (p1 + p2)^2, t = (p1 + p3)^2
box = graph("TL-TR, TL-BL, TR-BR, BL-BR", {"TL": "p1", "BL": "p2", "TR": "p3", "BR": "-p1-p2-p3"},
            kin={"p1^2": 0, "p2^2": 0, "p3^2": 0, "p1.p2": "s/2", "p1.p3": "t/2", "p2.p3": "-(s+t)/2"},
            euclidean=True)
fig = box.plot(figsize=(2.8, 2.8))
print("U =", box.U(), "     V =", box.F0())
""", {"fig": "fig"}),

("box2trees", """
def momentum_in(g, comp):
    tot = {}
    for v in comp:
        for k, c in g.ext.get(v, {}).items():
            tot[k] = tot.get(k, 0) + c
    return {k: c for k, c in tot.items() if c}

ft = box.two_forests()
titles = []
for rem, comp in ft:
    P2 = box.kin.square_external(momentum_in(box, comp))
    titles.append(xs(rem) + (":  0" if P2 == 0 else "  $%s$" % latex(P2)))
fig = draw_panels(box, [r for r, c in ft], titles=titles, momenta=False, ncols=3, size=2.0)
""", {"fig": "fig"}),

("kite", """
kite = graph("L-B, L-T, B-R, T-R, T-B", {"L": "p", "R": "-p"}, kin={"p^2": "pp"}, euclidean=True)
fig = kite.plot(figsize=(2.8, 2.8))
""", {"fig": "fig"}),

("kitetrees", """
trees = kite.spanning_trees()
fig = draw_panels(kite, trees, titles=[xs(r) for r in trees], momenta=False, ncols=8, size=1.5)
print(len(trees), "1-trees:  U =", kite.U())
""", {"fig": "fig"}),

("kite2trees", """
carry = [r for r, c in kite.two_forests() if momentum_in(kite, c)]
fig = draw_panels(kite, carry, titles=[xs(r) + "$\\\\,p^2$" for r in carry], momenta=False, ncols=8, size=1.5)
print(len(carry), "2-trees with p^2:  F =", kite.F())
""", {"fig": "fig"}),

("fpformula", """
# the massless bubble from U and F with the general formula: N = 2 lines, L = 1 loop, powers 1
from feynsage.oneloop import bubble_massless
g0 = graph("A-B, A-B", {"A": "p", "B": "-p"}, kin={"p^2": 1}, euclidean=True)
U0, F0 = g0.U(), g0.F()
Dv, N, L = 3.3, 2, 1
xx = var('xx')
u = SR(str(U0)).subs(x1=xx, x2=1 - xx); f = SR(str(F0)).subs(x1=xx, x2=1 - xx)
val = gamma(N - L*Dv/2) * numerical_integral(u^(N - (L + 1)*Dv/2)/f^(N - L*Dv/2), 0, 1)[0]
print("U =", U0, "   F =", F0)
print("from U and F:", val, "    closed form:", bubble_massless(1, 1, 1).subs(D=Dv).n())
""", {}),

("square", """
# the two-mass bubble, Euclidean: D1 = k^2 + m1^2, D2 = (k + p)^2 + m2^2
# write the scalar products as symbols: kk = k.k, kp = k.p, pp = p.p
x1, x2, m1, m2, kk, kp, pp = var('x1 x2 m1 m2 kk kp pp')
S = (x1*(kk + m1^2) + x2*(kk + 2*kp + pp + m2^2)).expand()       # sum_i x_i D_i
M = S.coefficient(kk)                                           # the coefficient of k.k
Qc = -S.coefficient(kp)/2                                        # -2 Q.k: Q = Qc * p
J = S.subs(kk=0, kp=0)                                          # what has no loop momentum
print("sum x_i D_i =", S)
print("M =", M, "    Q = (%s) p" % Qc, "    J =", J)
U_sq = M
F_sq = (M*(J - Qc^2*pp/M)).simplify_rational().expand()          # det M (J - Q.M^-1.Q), with Q.Q = Qc^2 p.p
print("U = det M =", U_sq)
print("F = det M (J - Q M^-1 Q) =", F_sq.collect(pp))
""", {}),

("squarecheck", """
fam2 = family([("k", "m1"), ("k + p", "m2")], kin={"p^2": "pp"}, euclidean=True)
U_f, F_f = fam2.UF()
print("fam.UF() gives the same:      ", bool(SR(str(U_f)) == U_sq), bool((SR(str(F_f)) - F_sq).expand() == 0))
print("and the 2-trees gave the same:", bool((SR(str(g2.F())) - F_sq).expand() == 0))
""", {}),

("sirnb", """
# sir's notebook completes the square, sum x_i D_i = l.M.l - 2 Q.l + J: U = det M, F = det M (J - Q.M^-1.Q)
tri = family(["k", "k + p1", "k + p1 + p2"], kin={"p1^2": "P1Sq", "p2^2": "P2Sq", "p1.p2": "(QSq - P1Sq - P2Sq)/2"})
bb = family([("k", "m1"), ("k + p1", "m2")], kin={"p1^2": "p1sq"})
bx = family(["k", "k + p1", "k + p1 + p2", "k - p3"],
            kin={"p1^2": 0, "p2^2": 0, "p3^2": 0, "p1.p2": "s/2", "p1.p3": "t/2", "p2.p3": "-(s+t)/2"})
for name, fam in [("triangle", tri), ("bubble", bb), ("box", bx)]:
    U_, F_ = fam.UF()
    print("%-9s U = %-18s F = %s" % (name, U_, F_))
""", {}),

]

# ------------------------------------------------------------------ a real calculation: pi0 -> gamma gamma
CELLS += [

("piongraph", """
# the quark triangle: P = pion vertex, A and B = photon vertices; all three lines are quarks of mass mq
# line 1 = P-A (Feynman parameter x), 2 = B-P (y), 3 = A-B (1 - x - y); Minkowski signs
tri = graph("P-A:mq, B-P:mq, A-B:mq", {"P": "k1 + k2", "A": "-k1", "B": "-k2"},
            kin={"k1^2": 0, "k2^2": 0, "k1.k2": "mpi2/2"})
fig = tri.plot(figsize=(2.6, 2.6))
print("U =", tri.U())
print("F =", tri.F())
""", {"fig": "fig"}),

("piontrees", """
trees = tri.spanning_trees(); two = tri.two_forests()
lab = lambda rem: "$" + "".join("x_{%d}" % (i + 1) for i in rem) + "$"
def term(rem, comp):
    P2 = tri.kin.square_external(momentum_in(tri, comp))
    return lab(rem) + ("$\\;m_\\pi^2$" if P2 != 0 else ":  0")
fig = draw_panels(tri, trees + [r for r, c in two],
                  titles=["1-tree " + lab(r) for r in trees] + ["2-tree " + term(r, c) for r, c in two],
                  momenta=False, ncols=3, size=2.0)
""", {"fig": "fig"}),

("piondelta", """
x, y = var('x y')
m, m_pi = var('m m_pi')
Delta = SR(str(tri.F())).subs(x1=x, x2=y, x3=1 - x - y).expand()      # U = 1 on the simplex
Delta = Delta.subs({SR.var('mq'): m, SR.var('mpi2'): m_pi^2})
Delta
""", {}),

("pionUF", """
# the general formula with N = 3 lines (powers 1), L = 1 loop, D = 4, Minkowski:
# Int d^4l/(i pi^2) 1/(D1 D2 D3) = (-1)^N Gamma(N - L D/2) Int dx delta(1 - sum x) U^(N-(L+1)D/2) / F^(N-L D/2)
N, L, Dv = 3, 1, 4
print("prefactor (-1)^N Gamma(N - L D/2) =", (-1)^N*gamma(N - L*Dv/2),
      "   power of U:", N - (L + 1)*Dv/2, "   power of F:", -(N - L*Dv/2))
Ux, Fx = SR(str(tri.U())), SR(str(tri.F()))
integrand = ((-1)^N*gamma(N - L*Dv/2) * Ux^(N - (L + 1)*Dv/2) / Fx^(N - L*Dv/2)).subs(x3=1 - x1 - x2)
integrand = integrand.subs({SR.var('mq'): m, SR.var('mpi2'): m_pi^2}).simplify_full()
print("integrand on the simplex x1 + x2 + x3 = 1:", integrand)
""", {}),

("pionnum", """
from scipy.integrate import dblquad
mv, rv = 1.0, 0.5                                 # m = 1 and m_pi^2 = r m^2 with r = 0.5
f_num = fast_callable(integrand.subs(m=mv, m_pi=sqrt(rv)), vars=[SR.var('x1'), SR.var('x2')])
val = dblquad(lambda b, a: f_num(a, b), 0, 1, 0, lambda a: 1 - a)[0]
print("Int d^4l/(i pi^2) 1/(D1 D2 D3) from U and F:  ", val)
print("-I(r)/m^2 with I(r) = (2/r) arcsin^2(sqrt(r)/2):", -(2/rv*arcsin(sqrt(rv)/2)^2)/mv^2)
""", {}),

("pionnorm", """
Ir = var('I_r')
# the loop measure of the amplitude is d^4l/(2 pi)^4 = (i pi^2/(2 pi)^4) d^4l/(i pi^2)
(I*pi^2/(2*pi)^4 * (-Ir/m^2)).simplify_full()
""", {}),

("pionIr", """
r = var('r')
series_I = sum(r^k * integrate(integrate((x*y)^k, y, 0, 1 - x), x, 0, 1) for k in range(10))
closed = 2/r * arcsin(sqrt(r)/2)^2
print("I(r) term by term:  ", series_I.series(r, 4).truncate())
print("closed form, series:", closed.taylor(r, 0, 3))
print("heavy quark, r -> 0:", limit(closed, r=0))
print("m = 330 MeV:  2 I(r) =", (2*closed.subs(r=(134.9768/330)^2)).n(digits=5))
""", {}),

("piontrace", """
# Tr[ g5 (q - x k1 - (1-y) k2 + m) g^nu (q - x k1 + y k2 + m) g^mu (q + (1-x) k1 + y k2 + m) ],
# each momentum combination standing for gamma.(that momentum); answer with the usual epsilon
form.dirac_trace(['g5', 'q - x*k1 - (1-y)*k2 + m', 'nu', 'q - x*k1 + y*k2 + m', 'mu',
                  'q + (1-x)*k1 + y*k2 + m'], vectors=['q', 'k1', 'k2'], levi_civita="usual")
""", {}),

("pionamp", """
Nc, e, Q, g, f, mq_ = var('N_c e Q g f_pi m_q')
A_one = Nc * e^2 * Q^2 * g / (4*pi^2*mq_) * 2*limit(closed, r=0)    # one quark, both orderings, heavy quark
print("one quark:             ", A_one)
print("with g = m_q/f_pi:     ", A_one.subs(g=mq_/f))
alpha = var('alpha')
A = (A_one.subs(g=mq_/f, Q=2/3) - A_one.subs(g=mq_/f, Q=-1/3)).subs(e=sqrt(4*pi*alpha))   # u minus d
print("u and d, e^2 = 4 pi alpha:", A.simplify_full())
""", {}),

("pioneps", """
# Sage-style input; feynsage writes the FORM program, runs it and reads the answer back
form.compute("eps(mu,nu,rho,sigma)*eps(mu,nu,al,be)*k1(rho)*k2(sigma)*k1(al)*k2(be)",
             vectors=["k1", "k2"], rules={"k1.k1": 0, "k2.k2": 0}, levi_civita="usual", show_code=True)
""", {}),

("pionrate", """
mpi = var('m_pi')
Gamma = (A^2 * mpi^3 / (64*pi)).simplify_full()
print("Gamma =", Gamma)
vals = {alpha: 1/137.035999, mpi: 134.9768, f: 130.2/sqrt(2)}          # MeV (PDG values)
hbar, BR = 6.582119569e-22, 0.98823                                   # MeV s, BR(pi0 -> gamma gamma)
G3 = (Gamma.subs(vals).subs({Nc: 3}) * 1e6).n()                        # eV
G1 = (Gamma.subs(vals).subs({Nc: 1}) * 1e6).n()
tau = lambda G: hbar*1e6*BR/G                                         # s, G in eV
print("width,    feynsage: %.2f eV          PrimEx-II (2020): 7.80 +- 0.12 eV   -> %.1f sigma" % (G3, abs(G3 - 7.80)/0.12))
print("lifetime, feynsage: %.3g s     from the PrimEx-II width:  %.3g s" % (tau(G3), tau(7.80)))
print("                                        PDG world average: (8.43 +- 0.13)e-17 s -> %.1f sigma"
      % (abs(tau(G3) - 8.43e-17)/0.13e-17))
print("the PDG lifetime as a width: %.2f eV (it includes older, lower measurements)" % (hbar*1e6*BR/8.43e-17))
print("with N_c = 1: %.2f eV, nine times too small" % G1)
vals_f = dict(vals); vals_f[f] = 1.01*130.2/sqrt(2)
print("f_pi 1%% larger: Gamma = %.2f eV (Gamma ~ 1/f_pi^2: 2%% lower)" % (Gamma.subs(vals_f).subs({Nc: 3})*1e6).n())
""", {}),

]
