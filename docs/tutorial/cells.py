# The code cells of the feynsage tutorial, in order.  They follow sir's two lectures "An
# Introduction to Feynman Integrals" (NISER, 2026), part 1 and part 2.  Each is (id, code, options); the cells
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
    print("n = %d:  I_n - (-1)^(n-1)/(n-1)! (d/dm^2)^(n-1) I_1 =" % k,
          (lhs - rhs).simplify_full())
""", {}),

("pole", """
expand_eps(tadpole(1, m2), 0)
""", {}),

("phi4", """
lam = var('lam', latex_name=r'\\lambda'); mu2 = var('mu2', latex_name=r'\\mu^2')
w, e = var('omega e')
# (1/2)(-lambda)(mu^2)^(2-w) Int d^(2w)l/(2 pi)^(2w) 1/(l^2 + m^2),
# the integral done as above
T = 1/2*(-lam)*mu2^(2 - w) * pi^w*gamma(1 - w)*m2^(w - 1)/(2*pi)^(2*w)
Texp = T.subs({w: 2 - e}).series(e, 1).truncate()                          # e = 2 - omega
lecture = lam*m2/(32*pi^2)*(1/e + psi(2) + log(4*pi*mu2/m2))
gap = (Texp - lecture).canonicalize_radical().simplify_full()
print("minus the lecture's result:", gap)
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
bub = family([("p", "m"), ("p - k", "m")], kin={"k^2": "kk"}, loops=["p"], euclidean=True,
             name="J")
bub.info()
""", {}),

("bubibp", """
for ident in bub.ibp((1, 1)):
    terms = ["(%s) J(%s)" % (c, ",".join(map(str, a))) for a, c in ident.items()]
    print(" + ".join(terms), "= 0")
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
# k.dJ/dk = J(0,2) - J(1,1) - k^2 J(1,2) and dJ/dk^2 = (k.dJ/dk)/(2 k^2);
# reduce the right side
kk = SR.var('kk'); dd = SR.var('d'); m = SR.var('m')
c = lambda a, mm_: SR(str(rb[a].get(mm_, 0)))
A = (c((0, 2), (1, 1)) - 1 - kk*c((1, 2), (1, 1)))/(2*kk)        # coefficient of J(1,1)
B = (c((0, 2), (0, 1)) - kk*c((1, 2), (0, 1)))/(2*kk)            # coefficient of T(1)
print("dJ/dk^2 = A J + B T(1) with")
print("  A =", A.factor())
print("  B =", B.factor())
print("as in the notes:",
      bool((A + (1/kk - (dd - 3)/(kk + 4*m^2))/2).simplify_full() == 0),
      bool((B + (dd - 2)/(4*m^2)*(1/kk - 1/(kk + 4*m^2))).simplify_full() == 0))
""", {}),

("dex", """
x = var('x')
# write k^2 = 4 m^2 x; then dJ/dx = 4 m^2 dJ/dk^2
Ax = (4*m^2*A).subs(kk=4*m^2*x).subs({dd: D}).simplify_full()
Bx = (4*m^2*B).subs(kk=4*m^2*x).subs({dd: D}).simplify_full()
print("dJ/dx = (%s) J" % Ax.partial_fraction(x))
print("        + (%s) T(1)" % Bx.partial_fraction(x))
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
    same = bool((ratio - hyp).simplify_full() == 0)
    print("c%d/c%d =" % (n + 1, n), ratio, "   same as 2F1:", same)
""", {}),

("denum", """
from scipy.integrate import solve_ivp
Dv, mv = 3.3, 1
T1v = float(tadpole(1, mv^2).subs(D=Dv))
fA = fast_callable(Ax.subs(m=mv, D=Dv), vars=[x])
fB = fast_callable(Bx.subs(m=mv, D=Dv), vars=[x])
x0 = 1e-6                                          # start next to the boundary x = 0
J_start = float(sum(sol[cs[i]].subs(D=Dv, m=mv, T1=T1v)*x0^i for i in range(N)))
num = solve_ivp(lambda t, y: [fA(t)*y[0] + fB(t)*T1v], [x0, 2.0], [J_start],
                rtol=1e-11, atol=1e-13, dense_output=True)
for xv in [0.125, 0.5, 2.0]:
    exact = bubble_equal_mass(4*mv^2*xv, mv^2).subs(D=Dv).n()   # Gamma(2-D/2) m^(D-4) 2F1
    print("x = %.3f   solved numerically: %.10f    closed form: %.10f"
          % (xv, num.sol(xv)[0], exact))
""", {}),

("hyp", """
a, b, c, z = 7/20, 1, 3/2, -1/2           # 2F1(2 - D/2, 1; 3/2; -x) at D = 3.3, x = 1/2
h = hypergeometric([a, b], [c], z)
print("Sage:          ", h.n(digits=25))
poch = rising_factorial                                   # (a)_n = a (a+1) ... (a+n-1)
term = lambda n: poch(a, n)*poch(b, n)/poch(c, n)*z^n/factorial(n)
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
# after the shift l -> l + p(1 - x) the loop integral is a tadpole
# with mass^2 = m^2 + p^2 x(1 - x)
tad_x = lambda xx: tadpole(2, mv^2 + pv*xx*(1 - xx)).subs(D=Dv).n()
fp = numerical_integral(tad_x, 0, 1)[0]
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
print("ln(-2)      =", log(-2).n())
print("ln 2 + i pi =", (log(2) + I*pi).n())
for sv in [3, 5, 10]:
    val = finite_part(B0(sv, 1, 1), 1).n(digits=12)
    beta = sqrt(1 - 4/sv) if sv > 4 else 0
    print("s = %2d (m = 1):  Im of the bubble = %.10f    pi*beta = %.10f"
          % (sv, imag(val), (pi*beta).n()))
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
titles = ["1-tree: " + xs(r) for r in trees] + ["2-tree: " + xs(r) for r in forests]
fig = draw_panels(g2, trees + forests, titles=titles, momenta=False, ncols=3, size=1.9)
print("U =", g2.U())
print("V =", g2.F0())
print("F = V + (m1^2 x1 + m2^2 x2) U =", g2.F())
""", {"fig": "fig"}),

("box", """
# Smirnov (2006) Fig. 3.6: line 1 on top, 2 on the left, 3 at the bottom, 4 on the right;
# p1 and p2 enter on the left, p3 and p4 on the right, s = (p1 + p2)^2, t = (p1 + p3)^2
box = graph("TL-TR, TL-BL, BL-BR, TR-BR",
            {"TL": "p1", "BL": "p2", "TR": "p3", "BR": "-p1-p2-p3"},
            kin={"p1^2": 0, "p2^2": 0, "p3^2": 0,
                 "p1.p2": "s/2", "p1.p3": "t/2", "p2.p3": "-(s+t)/2"},
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
fig = draw_panels(box, [r for r, c in ft], titles=titles, momenta=False, ncols=3,
                  size=2.0)
""", {"fig": "fig"}),

("kite", """
kite = graph("L-B, L-T, B-R, T-R, T-B", {"L": "p", "R": "-p"}, kin={"p^2": "pp"},
             euclidean=True)
fig = kite.plot(figsize=(2.8, 2.8))
""", {"fig": "fig"}),

("kitetrees", """
trees = kite.spanning_trees()
fig = draw_panels(kite, trees, titles=[xs(r) for r in trees], momenta=False, ncols=8,
                  size=1.5)
print(len(trees), "1-trees:  U =", kite.U())
""", {"fig": "fig"}),

("kite2trees", """
carry = [r for r, c in kite.two_forests() if momentum_in(kite, c)]
fig = draw_panels(kite, carry, titles=[xs(r) + "$\\\\,p^2$" for r in carry], momenta=False,
                  ncols=8, size=1.5)
print(len(carry), "2-trees with p^2:")
terms = str(kite.F()).split(" + ")
print("F =", " + ".join(terms[:4]))
print("  +", " + ".join(terms[4:]))
""", {"fig": "fig"}),

("fpformula", """
# the massless bubble from U and F with the general formula:
# N = 2 lines, L = 1 loop, powers 1
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
J = S.subs(kk=0, kp=0)                                     # what has no loop momentum
print("sum x_i D_i =", S)
print("M =", M, "    Q = (%s) p" % Qc, "    J =", J)
U_sq = M
# det M (J - Q.M^-1.Q), with Q.Q = Qc^2 p.p
F_sq = (M*(J - Qc^2*pp/M)).simplify_rational().expand()
print("U = det M =", U_sq)
print("F = det M (J - Q M^-1 Q) =", F_sq.collect(pp))
""", {}),

("squarecheck", """
fam2 = family([("k", "m1"), ("k + p", "m2")], kin={"p^2": "pp"}, euclidean=True)
U_f, F_f = fam2.UF()
print("fam.UF() gives the same:      ", bool(SR(str(U_f)) == U_sq),
      bool((SR(str(F_f)) - F_sq).expand() == 0))
print("and the 2-trees gave the same:", bool((SR(str(g2.F())) - F_sq).expand() == 0))
""", {}),

("sirnb", """
# sir's notebook completes the square, sum x_i D_i = l.M.l - 2 Q.l + J:
# U = det M, F = det M (J - Q.M^-1.Q)
tri = family(["k", "k + p1", "k + p1 + p2"],
             kin={"p1^2": "P1Sq", "p2^2": "P2Sq", "p1.p2": "(QSq - P1Sq - P2Sq)/2"})
bb = family([("k", "m1"), ("k + p1", "m2")], kin={"p1^2": "p1sq"})
bx = family(["k", "k + p1", "k + p1 + p2", "k - p3"],
            kin={"p1^2": 0, "p2^2": 0, "p3^2": 0,
                 "p1.p2": "s/2", "p1.p3": "t/2", "p2.p3": "-(s+t)/2"})
for name, fam in [("triangle", tri), ("bubble", bb), ("box", bx)]:
    U_, F_ = fam.UF()
    print("%-9s U = %-18s F = %s" % (name, U_, F_))
""", {}),

]

# ------------------------------------------------------------------ a real calculation: pi0 -> gamma gamma
CELLS += [

("piongraph", """
# the quark triangle: P = pion vertex, A and B = photon vertices;
# all three lines are quarks of mass mq
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
                  titles=["1-tree " + lab(r) for r in trees]
                         + ["2-tree " + term(r, c) for r, c in two],
                  momenta=False, ncols=3, size=2.0)
""", {"fig": "fig"}),

("piondelta", """
x, y = var('x y')
m, m_pi = var('m m_pi')
Delta = SR(str(tri.F())).subs(x1=x, x2=y, x3=1 - x - y).expand()   # U = 1 on the simplex
Delta = Delta.subs({SR.var('mq'): m, SR.var('mpi2'): m_pi^2})
Delta
""", {}),

("pionUF", """
# the general formula with N = 3 lines (powers 1), L = 1 loop, D = 4, Minkowski:
# Int d^4l/(i pi^2) 1/(D1 D2 D3)
#   = (-1)^N Gamma(N - L D/2) Int dx delta(1 - sum x) U^(N-(L+1)D/2) / F^(N-L D/2)
N, L, Dv = 3, 1, 4
print("prefactor (-1)^N Gamma(N - L D/2) =", (-1)^N*gamma(N - L*Dv/2),
      "   power of U:", N - (L + 1)*Dv/2, "   power of F:", -(N - L*Dv/2))
Ux, Fx = SR(str(tri.U())), SR(str(tri.F()))
integrand = (-1)^N*gamma(N - L*Dv/2) * Ux^(N - (L + 1)*Dv/2) / Fx^(N - L*Dv/2)
integrand = integrand.subs(x3=1 - x1 - x2)
integrand = integrand.subs({SR.var('mq'): m, SR.var('mpi2'): m_pi^2}).simplify_full()
print("integrand on the simplex x1 + x2 + x3 = 1:", integrand)
""", {}),

("pionnum", """
from scipy.integrate import dblquad
mv, rv = 1.0, 0.5                                 # m = 1 and m_pi^2 = r m^2 with r = 0.5
f_num = fast_callable(integrand.subs(m=mv, m_pi=sqrt(rv)),
                      vars=[SR.var('x1'), SR.var('x2')])
val = dblquad(lambda b, a: f_num(a, b), 0, 1, 0, lambda a: 1 - a)[0]
print("Int d^4l/(i pi^2) 1/(D1 D2 D3) from U and F:  ", val)
print("-I(r)/m^2 with I(r) = (2/r) arcsin^2(sqrt(r)/2):",
      -(2/rv*arcsin(sqrt(rv)/2)^2)/mv^2)
""", {}),

("pionnorm", """
Ir = var('I_r')
# the loop measure of the amplitude is d^4l/(2 pi)^4 = (i pi^2/(2 pi)^4) d^4l/(i pi^2)
(I*pi^2/(2*pi)^4 * (-Ir/m^2)).simplify_full()
""", {}),

("pionIr", """
r = var('r')
term_k = lambda k: integrate(integrate((x*y)^k, y, 0, 1 - x), x, 0, 1)
series_I = sum(r^k * term_k(k) for k in range(10))
closed = 2/r * arcsin(sqrt(r)/2)^2
print("I(r) term by term:  ", series_I.series(r, 4).truncate())
print("closed form, series:", closed.taylor(r, 0, 3))
print("heavy quark, r -> 0:", limit(closed, r=0))
print("m = 330 MeV:  2 I(r) =", (2*closed.subs(r=(134.9768/330)^2)).n(digits=5))
""", {}),

("piontrace", """
# the trace above, written as on paper; FORM does the work
q, k1, k2 = momenta("q k1 k2")
mu, nu = lorentz_indices("mu nu")
m = var('m')
dirac_trace(gamma5() * (slash(q - x*k1 - (1-y)*k2) + m) * gamma(nu)
            * (slash(q - x*k1 + y*k2) + m) * gamma(mu)
            * (slash(q + (1-x)*k1 + y*k2) + m))
""", {}),

("pionamp", """
Nc, e, Q, g, f, mq_ = var('N_c e Q g f_pi m_q')
# one quark, both orderings of the photons, heavy quark
A_one = Nc * e^2 * Q^2 * g / (4*pi^2*mq_) * 2*limit(closed, r=0)
print("one quark:             ", A_one)
print("with g = m_q/f_pi:     ", A_one.subs(g=mq_/f))
alpha = var('alpha')
A = A_one.subs(g=mq_/f, Q=2/3) - A_one.subs(g=mq_/f, Q=-1/3)          # u minus d
A = A.subs(e=sqrt(4*pi*alpha))
print("u and d, e^2 = 4 pi alpha:", A.simplify_full())
""", {}),

("pioneps", """
# the contraction above, written as on paper
contract(epsilon(mu, nu, k1, k2) * epsilon(mu, nu, k1, k2),
         rules={dot(k1, k1): 0, dot(k2, k2): 0}, debug=True)
""", {}),

("pionrate", """
mpi = var('m_pi')
Gamma = (A^2 * mpi^3 / (64*pi)).simplify_full()
print("Gamma =", Gamma)
vals = {alpha: 1/137.035999, mpi: 134.9768, f: 130.2/sqrt(2)}          # MeV
G3 = (Gamma.subs(vals).subs({Nc: 3}) * 1e6).n()                        # eV
G1 = (Gamma.subs(vals).subs({Nc: 1}) * 1e6).n()
print("width     feynsage %.2f eV" % G3)
print("with N_c = 1: %.2f eV, nine times too small" % G1)
vals_f = dict(vals); vals_f[f] = 1.01*130.2/sqrt(2)
Gf = (Gamma.subs(vals_f).subs({Nc: 3})*1e6).n()
print("f_pi 1%% larger: Gamma = %.2f eV, 2%% lower (Gamma ~ 1/f_pi^2)" % Gf)
""", {}),

]

CELLS += [
("fa_tops", """
tops = topologies(0, 2, 2)                 # tree level, 2 -> 2
len(tops), [len(topologies(L, a, b)) for (L, a, b) in ((0, 2, 3), (1, 2, 2), (2, 1, 1))]
""", {}),

("fa_insert", """
sm = SM()                                  # FeynArts' {"SM", "SMQCD"}, Feynman gauge
diags = insert_fields(tops, ["e-", "e+"], ["mu-", "mu+"], sm)
diags.counts(), [d.propagator_fields() for d in diags]
""", {}),

("fa_draw", """
fig = diags.draw(size=1.9)
""", {"fig": "fig"}),

("fa_loop", """
one_loop = insert_fields(topologies(1, 2, 2, exclude=("tadpoles", "wf")), ["e-", "e+"], ["mu-", "mu+"], sm)
len(one_loop)
""", {}),

("fa_m2", """
from feynsage.models import EL, SW, CW, MW, MZ, MH, ME, MM
p1, p2, p3, p4 = momenta("p1 p2 p3 p4")
s, t, u = var('s t u')
mandelstam = {dot(p1, p1): ME^2, dot(p2, p2): ME^2, dot(p3, p3): MM^2, dot(p4, p4): MM^2,
              dot(p1, p2): (s - 2*ME^2)/2, dot(p3, p4): (s - 2*MM^2)/2,
              dot(p1, p3): (ME^2 + MM^2 - t)/2, dot(p2, p4): (ME^2 + MM^2 - t)/2,
              dot(p1, p4): (ME^2 + MM^2 - u)/2, dot(p2, p3): (ME^2 + MM^2 - u)/2}
M2 = (diags.squared([p1, p2], [p3, p4]) / 4).subs(mandelstam)       # spins summed, 1/4 for the average
len(M2.expand().operands())
""", {"text": True}),

("fa_numbers", """
values = {EL: sqrt(4*pi/137.036), ME: 0.000511, MM: 0.105658, MZ: 91.1876, MW: 80.379, MH: 125.25,
          SW: sqrt(0.23122), CW: sqrt(1 - 0.23122)}
def point(rs, c, me=0.000511, mmu=0.105658):
    sv = rs^2; pin = sqrt(sv/4 - me^2); pout = sqrt(sv/4 - mmu^2)
    tv = me^2 + mmu^2 - sv/2 + 2*pin*pout*c
    return {s: sv, t: tv, u: 2*me^2 + 2*mmu^2 - sv - tv}
[M2.subs(values).subs(point(70, c)).n(digits=9) for c in (0.5, -0.5)]
""", {}),

("fa_gauge", """
M2u = (diags.squared([p1, p2], [p3, p4], gauge="unitary") / 4).subs(mandelstam)
(M2u - M2).subs({MW: MZ*CW}).subs({SW: sqrt(1 - CW^2)}).subs(u=2*ME^2 + 2*MM^2 - s - t).simplify_rational()
""", {}),

("fa_sigma", """
GZ = var('GZ')
M2w = (diags.squared([p1, p2], [p3, p4], widths={"Z": GZ, "G0": GZ}) / 4).subs(mandelstam)
f = M2w.subs(values).subs({GZ: 2.4952})
c = var('c')
def sigma_pb(rs, me=0.000511, mmu=0.105658):
    sv = rs^2; pin = sqrt(sv/4 - me^2); pout = sqrt(sv/4 - mmu^2)
    g = fast_callable(f.subs(point(rs, c)), vars=[c], domain=CDF)
    return 3.894e8/(32*pi.n()*sv) * (pout/pin) * numerical_integral(lambda x: g(x).real(), -1, 1)[0]
[round(sigma_pb(rs), 2) for rs in (20, 60, 91.1876, 160)]
""", {}),

("fa_plot", """
import numpy as np, matplotlib.pyplot as plt
roots = np.linspace(20, 160, 141)
fig, ax = plt.subplots(figsize=(5.2, 3.0))
ax.semilogy(roots, [sigma_pb(x) for x in roots], color="C3", lw=1.6, label=r"$\\gamma + Z + H + G^0$")
ax.semilogy(roots, [4*np.pi/(3*x**2)/137.036**2*3.894e8 for x in roots], color="C0", ls="--", lw=1, label=r"$4\\pi\\alpha^2/(3s)$")
ax.set_xlabel(r"$\\sqrt{s}$ (GeV)"); ax.set_ylabel(r"$\\sigma$ (pb)"); ax.legend(fontsize=8)
fig.tight_layout()
""", {"fig": "fig"}),
]

CELLS += [
("h_setup", """
from feynsage.models import EL, SW, CW, MW, MZ, MH, MT, MB, GS
from feynsage.pv import finite_part, uv_part
from feynsage import metric as g_, comp as c_, dot as dot_
sm = SM()
k1, k2 = momenta("k1 k2")                      # the two decay products; the Higgs has k1 + k2
peskin = {EL: sqrt(4*pi/128), SW: sqrt(0.23), CW: sqrt(0.77), MW: 80, MZ: 91, MT: 175, MB: 5, GS: sqrt(4*pi*0.12)}
def width(final, masses, identical=False):
    # Gamma(h -> final) at tree level: diagrams, |M|^2 summed over spins, polarizations and colours
    dg = insert_fields(topologies(0, 1, 2), ["H"], final, sm)
    r = dg.squared([k1 + k2], [k1, k2])
    for leg in (2, 3):
        e = momenta("eps%d" % leg)
        if dg[0].model.field(final[leg - 2]).kind == "V":
            r = polarization_sum(r, e)
    r = color_factor(r, N=3)
    m1, m2 = masses
    r = r.subs({dot_(k1, k1): m1^2, dot_(k2, k2): m2^2, dot_(k1, k2): (MH^2 - m1^2 - m2^2)/2})
    p = sqrt((MH^2 - (m1 + m2)^2)*(MH^2 - (m1 - m2)^2))/(2*MH)
    return (p/(8*pi*MH^2) * r * (1/2 if identical else 1)).simplify_full()
G_bb = width(["b", "b~"], (MB, MB))
G_bb
""", {}),

("h_ff_check", """
al = var('al', latex_name=r'\\alpha')
assume(MB > 0, MW > 0, MZ > 0, MH > 2*MB, MH > 2*MW, MH > 2*MZ)
peskin_a = al*MH/(8*SW^2) * MB^2/MW^2 * (1 - 4*MB^2/MH^2)^(3/2) * 3
(G_bb.subs(EL=sqrt(4*pi*al)) - peskin_a).canonicalize_radical().simplify_full()
""", {}),

("h_vv", """
G_WW = width(["W+", "W-"], (MW, MW))
G_ZZ = width(["Z", "Z"], (MZ, MZ), identical=True)
tW, tZ = MH^2/MW^2, MH^2/MZ^2
x1559 = al*MH^3/(16*MW^2*SW^2) * (1 - 4/tW + 12/tW^2) * sqrt(1 - 4/tW)
x1560 = al*MH^3/(32*MW^2*SW^2) * (1 - 4/tZ + 12/tZ^2) * sqrt(1 - 4/tZ)
[(G.subs(EL=sqrt(4*pi*al)).subs(CW=MW/MZ) - x).canonicalize_radical().simplify_full() for G, x in ((G_WW, x1559), (G_ZZ, x1560))]
""", {}),

("h_goldstone", """
G_GG = width(["G+", "G-"], (0, 0))          # the Goldstone bosons of the ungauged theory (massless)
[(G_WW/G_GG).subs(peskin).subs(MH=m).n(digits=6) for m in (300, 1000, 5000)]
""", {}),

("h_gg", """
light = ['e', 'mu', 'ta', 'ne', 'nm', 'nt', 'u', 'c', 'd', 's']
dgg = insert_fields(topologies(1, 1, 2, exclude=("tadpoles", "wf")), ["H"], ["g", "g"], sm, exclude_fields=light)
def transverse(num, ext):                     # T_{mu nu} M^{mu nu} = (d - 2) (k1.k2) F
    mu, nu = ext[1], ext[2]
    return num*(g_(mu, nu) - (c_(k2, mu)*c_(k1, nu) + c_(k1, mu)*c_(k2, nu))/dot_(k1, k2))
kin = {"k1^2": 0, "k2^2": 0, "k1.k2": "MH^2/2"}
Mgg = sum(dgg.loop_amplitude([k1 + k2], [k1, k2], transverse, kin))
a2, a3 = color_indices("a2 a3")
F_gg = (finite_part(Mgg, 1) / MH^2).subs({color_delta(a2, a3): 1})      # F = T.M/((d - 2) k1.k2); colour delta^{ab} stripped
len(dgg), uv_part(Mgg).simplify_full()
""", {}),

("h_I", """
import mpmath
def I_x(tau):
    # Xianyu (21.67), I = 3 int dx dy (1 - 4xy)/(1 - xy tau - i0), with the y integral done by hand:
    # 3 int_0^1 dx [ (4/tau)(1 - x) - (1 - 4/tau) log(1 - x(1 - x) tau - i0)/(x tau) ]
    tau = mpmath.mpf(float(tau))
    lg = lambda x: mpmath.log(mpmath.mpc(1 - x*(1 - x)*tau, -mpmath.mpf(10)**-40))
    pts = [0, 0.5, 1]
    if tau > 4:                                   # the zeros of 1 - x(1 - x) tau, where the log is singular
        r = mpmath.sqrt(1 - 4/tau)/2
        pts = [0, 0.5 - r, 0.5, 0.5 + r, 1]
    return complex(3*mpmath.quad(lambda x: 4/tau*(1 - x) - (1 - 4/tau)*lg(x)/(x*tau), pts))
G_gg = abs(F_gg)^2 * MH^3/(8*pi)        # sum |M|^2 = 8 colours x 2 (k1.k2)^2 |F|^2, Gamma = sum |M|^2/(32 pi m_h)
def peskin_c(m):
    I = sum(I_x(float(m^2/mq^2)) for mq in (5, 175))
    return float((1/128)*m/(8*0.23)*m^2/80^2*0.12^2/(9*pi.n()^2)*abs(complex(I))^2)
[(float(G_gg.subs(peskin).subs(MH=m).n()), peskin_c(m)) for m in (100, 300, 400)]
""", {}),

("h_aa", """
fermions_light = ['e', 'mu', 'ne', 'nm', 'nt', 'u', 'c', 'd', 's']
daa = insert_fields(topologies(1, 1, 2, exclude=("tadpoles", "wf")), ["H"], ["A", "A"], sm, exclude_fields=fermions_light + ['ta'])
bosonic = DiagramList([x for x in daa if not any(sm.field(f).kind == "F" for f in x.propagator_fields())], daa.tops, [])
fermionic = DiagramList([x for x in daa if any(sm.field(f).kind == "F" for f in x.propagator_fields())], daa.tops, [])
kin = {"k1^2": 0, "k2^2": 0, "k1.k2": "MH^2/2"}
M_W = sum(bosonic.loop_amplitude([k1 + k2], [k1, k2], transverse, kin))
M_f = sum(fermionic.loop_amplitude([k1 + k2], [k1, k2], transverse, kin))
len(daa), len(bosonic), uv_part(M_W).simplify_full(), uv_part(M_f).simplify_full()
""", {}),

("h_aa_draw", """
fig = bosonic.draw(ncols=7, size=1.25)
""", {"fig": "fig"}),

("h_IW", """
def I_W(tau):                                  # Xianyu (21.97)-(21.100)
    tau = mpmath.mpf(float(tau))
    I1 = mpmath.quad(lambda x: mpmath.log(1 - x*(1 - x)*tau), [0, 1])
    I2 = 2*mpmath.quad(lambda x: mpmath.quad(lambda y: mpmath.log(1 - x*y*tau), [0, 1 - x]), [0, 1])
    I3 = mpmath.quad(lambda x: mpmath.quad(lambda y: (8 - 3*x + y + 4*x*y)*tau/(1 - x*y*tau), [0, 1 - x]), [0, 1])
    return complex((6*I1 - 8*I2 + tau*(I1 - I2) + I3)/tau)
v_ev, alpha = 2*MW*SW/EL, EL^2/(4*pi)
R_W = (finite_part(M_W, 1)/MH^2) * pi*v_ev/alpha        # = I_W if M_W = alpha m_h^2/(2 pi v) I_W (below the W W threshold)
[(float(R_W.subs(peskin).subs(MH=m).n().real()), I_W(m^2/80^2).real) for m in (50, 100, 150)]
""", {}),

("h_light", """
R_f = (finite_part(M_f, 1)/MH^2) * pi*v_ev/alpha
small = {MH: 1/10}
R_W.subs(peskin).subs(small).n(digits=40).real().n(digits=8), R_f.subs(peskin).subs(small).n(digits=40).real().n(digits=8)
""", {}),

("h_aa_width", """
F_aa = (finite_part(M_W + M_f, 1)/MH^2).subs(peskin)
def compare(m):
    G = float(abs(complex(F_aa.subs(MH=m).n()))^2 * m^3/(64*pi.n()))
    Sf = (4/9)*3*complex(I_x(m^2/175^2)) + (1/9)*3*complex(I_x(m^2/5^2))
    pref = (1/128)*m/(8*0.23)*m^2/80^2*(1/128)^2/(18*pi.n()^2)
    IW = I_W(m^2/80^2)
    return G, float(pref*abs(Sf - 1.5*IW)^2), float(pref*abs(Sf - IW)^2)
[compare(m) for m in (60, 100, 150)]
""", {}),

("h_pp", """
def sigma_pp_nb(rs, mh=30, I=1):
    # sigma(pp -> h + X) = int dx1 dx2 f_g(x1) f_g(x2) sigma(gg -> h), f_g(x) = 8 (1 - x)^7/x,
    # sigma(gg -> h) = pi^2/(8 m_h) Gamma(h -> gg) delta(s_hat - m_h^2)   (Xianyu (21.71))
    G = float((1/128)*mh/(8*0.23)*mh^2/80^2*0.12^2/(9*pi.n()^2)*I^2)
    tau0 = mh^2/rs^2
    lum = mpmath.quad(lambda x1: (8*(1 - x1)**7/x1)*(8*(1 - tau0/x1)**7/(tau0/x1))/x1, [tau0, 1])
    return float(pi.n()^2/(8*mh)*G*lum/rs^2 * 0.3894e6)      # GeV^-2 -> nb
roots = [1000*x for x in (1, 2, 5, 10, 20, 40)]
[(r/1000, round(sigma_pp_nb(r), 3)) for r in roots]
""", {}),

("h_plots", """
import numpy as np, matplotlib.pyplot as plt
mhs = np.linspace(50, 500, 46)
I_top = [complex(I_x(m**2/175**2)) for m in mhs]
G_gg_keV = [1e6*float(G_gg.subs(peskin).subs(MH=m).n()) for m in mhs]
fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.4, 2.9))
a1.plot(mhs, [z.real for z in I_top], label=r"Re $I(m_h^2/m_t^2)$"); a1.plot(mhs, [z.imag for z in I_top], ls="--", label=r"Im")
a1.set_xlabel(r"$m_h$ (GeV)"); a1.legend(fontsize=8)
a2.semilogy(mhs, G_gg_keV, color="C3"); a2.set_xlabel(r"$m_h$ (GeV)"); a2.set_ylabel(r"$\\Gamma(h \\to gg)$ (keV)")
fig.tight_layout()
""", {"fig": "fig"}),

("h_br", """
Gs = {"bb": G_bb.subs(EL=sqrt(4*pi*al)).subs(al=1/128), "WW": G_WW, "ZZ": G_ZZ, "tt": width(["t", "t~"], (MT, MT)), "gg": G_gg}
def partial(name, m):
    thr = {"WW": 2*80, "ZZ": 2*91, "tt": 2*175, "bb": 10, "gg": 0}[name]
    return 0.0 if m <= thr else float(abs(complex(Gs[name].subs(peskin).subs(MH=m).n())))
table = {n: [partial(n, m) for m in mhs] for n in Gs}
total = [sum(table[n][i] for n in Gs) for i in range(len(mhs))]
fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.4, 3.0))
a1.semilogy(mhs, total, color="k"); a1.set_xlabel(r"$m_h$ (GeV)"); a1.set_ylabel(r"$\\Gamma_h$ (GeV)")
labels = {"bb": r"$b\\bar b$", "tt": r"$t\\bar t$", "gg": r"$gg$", "WW": r"$W^+W^-$", "ZZ": r"$Z^0Z^0$"}
colours = {"bb": "#7b2d8e", "tt": "#d62728", "gg": "#ff7f0e", "WW": "#1f5fbf", "ZZ": "#2ca02c"}
for n in Gs:
    a2.semilogy(mhs, [max(table[n][i]/total[i], 1e-6) for i in range(len(mhs))], label=labels[n], color=colours[n])
print("BR(gg) at m_h = 100, 200, 300, 400 GeV:", [round(table["gg"][i]/total[i], 5) for i in (5, 15, 25, 35)])
a2.set_ylim(1e-3, 1.2); a2.set_xlabel(r"$m_h$ (GeV)"); a2.set_ylabel("branching fraction"); a2.legend(fontsize=7)
fig.tight_layout()
""", {"fig": "fig"}),

]

CELLS += [
("c21_qed", """
from feynsage.models import MM
from feynsage import u as u_spinor
p1, p2 = momenta("p1 p2")                       # muon in (p1), photon in (p2 - p1), muon out (p2)
q2 = var("q2")                                  # q^2 = (p2 - p1)^2, sent to 0 at the end
def F2_projector(num, ext):                     # Tr[(p1/ + m) P (p2/ + m) X] picks F2 out of X^mu
    P = vertex_projector(p1, p2, MM, ext[1], "F2")
    return dirac_trace((slash(p1) + MM)*P*(slash(p2) + MM)*num, dim="d", gamma5_scheme="NDR-even")
vertex_kin = {"p1^2": "MM^2", "p2^2": "MM^2", "p1.p2": "MM^2 - q2/2"}
def at_q2_zero(F):
    return ((F.subs(q2=10^-10) + F.subs(q2=-10^-10))/2)
tops_v = topologies(1, 2, 1, exclude=("tadpoles", "wf"))
ALL = ["A", "Z", "W+", "H", "G0", "G+", "u+", "u-", "uA", "uZ", "g", "ug", "e", "mu", "ta", "ne", "nm", "nt",
       "u", "c", "t", "d", "s", "b"]
only = lambda keep: [f for f in ALL if f not in keep]   # fields that may not run inside the loop
dq = insert_fields(tops_v, ["mu", "A"], ["mu"], sm, exclude_fields=only(["A", "mu"]))
F2 = finite_part(sum(dq.loop_amplitude([p1, p2 - p1], [p2], F2_projector, kin=vertex_kin)), 1)/EL
F2_qed = at_q2_zero(F2.subs({EL: 1, MM: 1})).n(digits=40)
F2_qed.n(digits=20), (1/(8*pi^2)).n(digits=20)     # F2(0)/e^2 and alpha/(2 pi) with e = 1
""", {}),
("c21_w", """
GF = EL^2/(4*sqrt(2)*SW^2*MW^2)                 # G_F/sqrt(2) = g^2/(8 m_W^2), g = e/sin(theta_w)
unit = GF*MM^2/(8*pi^2*sqrt(2))
heavy = {EL: 1, MM: 1, MW: 1000, SW: sqrt(23/100), CW: sqrt(77/100), MZ: 1000/sqrt(77/100)}   # m_W = 1000 m_mu
dW = insert_fields(tops_v, ["mu", "A"], ["mu"], sm, exclude_fields=only(["W+", "G+", "nm"]))
def a_mu(diags, xi=1):
    res = diags.loop_amplitude([p1, p2 - p1], [p2], F2_projector, kin=vertex_kin, xi=xi)
    return [(at_q2_zero(finite_part(r, 1)/EL/unit).subs(heavy)).n(digits=60).n(digits=8) for r in res]
[x.propagator_fields() for x in dW], a_mu(dW)
""", {}),
("c21_xi", """
[(xi, sum(a_mu(dW, xi))) for xi in (1, 3, 1/2, 1/100, 100)]
""", {}),
("c21_z", """
dZ = insert_fields(tops_v, ["mu", "A"], ["mu"], sm, exclude_fields=only(["Z", "G0", "mu"]))
sw2 = 23/100
[x.propagator_fields() for x in dZ], sum(a_mu(dZ)), -(4/3 + 8/3*sw2 - 16/3*sw2^2).n(digits=8)
""", {}),
("c21_ww", """
from feynsage.models import ME
from feynsage.qft import tensors as T
from feynsage.diagrams import conjugate_amplitude
dww = insert_fields(topologies(0, 2, 2), ["e-", "e+"], ["W+", "W-"], sm)
q1, q2_, q3, q4 = momenta("q1 q2 q3 q4")       # e-(q1) e+(q2) -> W+(q3) W-(q4)
Mww = dww.amplitude([q1, q2_], [q3, q4])
s_, t_, u_ = var("s t u")
kin_ww = {dot(q1,q1): 0, dot(q2_,q2_): 0, dot(q3,q3): MW^2, dot(q4,q4): MW^2, dot(q1,q2_): s_/2,
          dot(q3,q4): (s_ - 2*MW^2)/2, dot(q1,q3): (MW^2 - t_)/2, dot(q2_,q4): (MW^2 - t_)/2,
          dot(q1,q4): (MW^2 - u_)/2, dot(q2_,q3): (MW^2 - u_)/2}
# longitudinal polarization vectors, real:  eps_L(q3) = [q3 (q3.q4) - q4 m_W^2]/(m_W sqrt((q3.q4)^2 - m_W^4))
P34 = (s_ - 2*MW^2)/2; NL = MW*sqrt(P34^2 - MW^4)
eL3, eL4 = (P34*q3 - MW^2*q4)/NL, (P34*q4 - MW^2*q3)/NL
longitudinal = {T._declare(n, "vector"): v for n, v in (("eps3", eL3), ("eps3_c", eL3), ("eps4", eL4), ("eps4_c", eL4))}
def squared_00(hand):                           # |M(e-_hand e+ -> W+_0 W-_0)|^2
    Mh = chiral(Mww, q1, hand)
    return SR(spin_sum(Mh*conjugate_amplitude(Mh), rules=longitudinal)).subs(kin_ww).subs({ME: 0})
LR, RL = squared_00("L"), squared_00("R")
mw_ = 804/10; sw2 = 23/100; mz_ = mw_/sqrt(1 - sw2)
num_ww = {EL: 1, SW: sqrt(sw2), CW: sqrt(1 - sw2), MW: mw_, MZ: mz_}
def pt_ww(rs, c):                               # theta = angle between e- and W+
    E = rs/2; p = sqrt(E^2 - mw_^2); tv = mw_^2 - rs^2/2 + 2*E*p*c
    return {s_: rs^2, t_: tv, u_: 2*mw_^2 - rs^2 - tv}
LR.subs(num_ww).subs(pt_ww(500, -1/2)).n(digits=20)
""", {}),
("c21_xianyu", """
def X24_first(rs, c):                           # first line of Xianyu (21.24), without i e^2
    sv = rs^2; E = rs/2; p = sqrt(E^2 - mw_^2); sn = sqrt(1 - c^2); uv = mw_^2 - sv/2 - 2*E*p*c
    return (mz_^2/(sv*(sv - mz_^2)) - 1/(2*sw2)/(sv - mz_^2))*(-4*E*p*(E^2 + p^2)/mw_^2 + 16*E^3*p/mw_^2)*sn + 1/(2*sw2)/uv*2*E*(-3*E^2*p + p^3 - 2*E^3*c)*sn/mw_^2
def X24_second(rs, c, sign):                    # second line of (21.24); printed with sign = -1
    sv = rs^2; b = sqrt(1 - 4*mw_^2/sv); D = 1 + b^2 + 2*b*c
    return -sv/(4*mw_^2)*(mz_^2/(sv - mz_^2)*b*(3 - b^2) + sign/(2*sw2)*((2/D - sv/(sv - mz_^2))*b*(3 - b^2) + 4*c/D))*sqrt(1 - c^2)
def X29(rs, c):                                 # Xianyu (21.29)
    sv = rs^2; b = sqrt(1 - 4*mw_^2/sv)
    return sv/(sv - mz_^2)*mz_^2/(4*mw_^2)*b*(b^2 - 3)*sqrt(1 - c^2)
rows = []
for rs, c in [(200, 3/10), (500, -1/2), (1000, 4/5)]:
    a, r = LR.subs(num_ww).subs(pt_ww(rs, c)), RL.subs(num_ww).subs(pt_ww(rs, c))
    rows.append([rs, c] + [(a/x^2).n(digits=60).n(digits=15) for x in (X24_first(rs, c), X24_second(rs, c, -1), X24_second(rs, c, +1))]
                + [(r/X29(rs, c)^2).n(digits=60).n(digits=15)])
print("sqrt(s)  cos   first line     printed 2nd    2nd with +     (21.29)")
for row in rows:
    print("%-8s %-5s " % (row[0], row[1]) + "  ".join("%-14.10f" % float(x) for x in row[2:]))
""", {}),
("c21_peskin", """
J = chiral(vbar(q2_)*slash(q3 - q4)*u_spinor(q1), q1, "L")      # vbar_L gamma_l u_L (k+ - k-)^l
J2 = SR(spin_sum(J*conjugate_amplitude(J))).subs(kin_ww)
def P108(rs, cP):                               # the bracket of Peskin (21.108); theta_P = angle of W-
    sv = rs^2; b = sqrt(1 - 4*mw_^2/sv)
    return (1/(2*sw2)*(-sv/(sv - mz_^2)*(mz_^2/(2*mw_^2) + 1) + 2/b^2 - 8*mw_^2/(sv*b^2*(1 + b^2 - 2*b*cP)))
            + mz_^2/mw_^2*((sv/2 + mw_^2)/(sv - mz_^2)))/sv
[(rs, c, (LR.subs(num_ww).subs(pt_ww(rs, c))/(J2.subs(num_ww).subs(pt_ww(rs, c))*P108(rs, -c)^2)).n(digits=60).n(digits=15))
 for rs, c in [(200, 3/10), (500, -1/2), (1000, 4/5)]]
""", {}),
("c21_equiv", """
dgg = insert_fields(topologies(0, 2, 2), ["e-", "e+"], ["G+", "G-"], sm)
Mg = chiral(dgg.amplitude([q1, q2_], [q3, q4]), q1, "L")
GG = SR(spin_sum(Mg*conjugate_amplitude(Mg))).subs(kin_ww).subs({ME: 0})
c = 3/10
def r(x, y, rs):
    return (x.subs(num_ww).subs(pt_ww(rs, c)).n(digits=80)/y.subs(num_ww).subs(pt_ww(rs, c)).n(digits=80)).n(digits=8)
lim = (1 - c^2)/(4*sw2*(1 - sw2))^2                         # (e^2 sin(theta)/(4 s_w^2 c_w^2))^2
X28 = (1 + 2*c)^2*(1 - c^2)/(2*sw2*(1 + c))^2               # Xianyu (21.28) squared
print("sqrt(s)/GeV   W0W0/phi phi   W0W0/limit   W0W0/(21.28)")
for rs in (10^3, 10^4, 10^5, 10^6):
    x = LR.subs(num_ww).subs(pt_ww(rs, c)).n(digits=80)
    print("%-13d %-14.8f %-12.8f %.8f" % (rs, float(r(LR, GG, rs)), float(x/lim), float(x/X28)))
""", {}),
]

CELLS += [
("c21_fp", """
# Xianyu's Feynman-parameter forms at q^2 = 0 with the full Delta, at the real m_W/m_mu
import mpmath
mpmath.mp.dps = 30
r_real = (80379/1000)/(1056584/10^7)
real = {EL: 1, MM: 1, MW: r_real, SW: sqrt(sw2), CW: sqrt(1 - sw2), MZ: r_real/sqrt(1 - sw2)}
def a_real(diags):
    res = diags.loop_amplitude([p1, p2 - p1], [p2], F2_projector, kin=vertex_kin)
    return [(at_q2_zero(finite_part(r, 1)/EL/unit).subs(real)).n(digits=40) for r in res]
def dbl(f):
    return mpmath.quad(lambda x: mpmath.quad(f, [0, 1 - x]), [0, 1])
R2, z2, c2 = mpmath.mpf(float(r_real))**2, mpmath.mpf(float(r_real))**2/(1 - mpmath.mpf(0.23)), 1 - mpmath.mpf(0.23)
x21_2 = 2*R2*dbl(lambda y: (1 - y)*(3 - 2*y)/((1 - y)*R2 - y*(1 - y)))                # (21.2)
x21_9 = -(R2/(2*c2))*dbl(lambda y: (2*y*(3 + y) - (4*mpmath.mpf(0.23) - 1)**2*2*y*(1 - y))/((1 - y)**2 + y*z2))   # (21.9)
[(a_real(dW)[3].n(digits=12), mpmath.nstr(x21_2, 12)), (a_real(dZ)[1].n(digits=12), mpmath.nstr(x21_9, 12))]
""", {}),
("c21_trans_setup", """
# every helicity state with process(): electrons massless, the W helicities +1, -1, 0 (Jacob-Wick vectors)
W = process("e- e+ -> W+ W-", massless=["e"])
num_W = {EL: 1, SW: sqrt(sw2), CW: sqrt(1 - sw2), MW: mw_, MZ: mz_}       # e = 1, as above
def hel(he, h3, h4, rs, c, diagrams=None):           # he = helicity of the e- (the e+ has the opposite)
    r = W.squared(helicities={"e-": he, "e+": -he, "W+": h3, "W-": h4}, diagrams=diagrams)
    return r.subs(num_W).subs({W.sqrt_s: rs, W.cos_theta: c}).n(digits=60)
hel(-1, 0, 0, 500, -1/2)                             # the same number as Out above for e-_L e+_R -> W0 W0
""", {}),
("c21_trans", """
def X(name, rs, c, pm):              # Xianyu's printed amplitudes without i e^2 (first and second lines)
    sv = rs^2; E0 = rs/2; p = sqrt(E0^2 - mw_^2); b = p/E0; sth = sqrt(1 - c^2); D = 1 + b^2 + 2*b*c
    uv = mw_^2 - 2*E0^2 - 2*E0*p*c
    br = mz_^2/(sv*(sv - mz_^2)) - 1/(2*sw2)/(sv - mz_^2)
    return {"25a": br*(8*E0^2*p/mw_)*(-pm + c)/sqrt(2) - 1/(2*sw2)/uv*(2*E0/mw_)*(E0^2*(2*c - pm) + 2*E0*p + pm*p^2)*(pm + c)/sqrt(2),
            "25b": (mz_^2/(sv - mz_^2)*b - 1/(2*sw2)*(sv/(sv - mz_^2)*b + (pm - 2*c - 2*b - pm*b^2)/D))*(rs/mw_)*(pm + c)/sqrt(2),
            "26a": br*(-4*E0*p*sth) + 1/(2*sw2)/uv*2*E0*(p + E0*c)*sth,
            "26b": (-mz_^2/(sv - mz_^2)*b + 1/(2*sw2)*(sv/(sv - mz_^2)*b - 2*(b + c)/D))*sth,
            "27a": -1/(2*sw2)/uv*2*E0^2*(-pm + c)*sth, "27b": 1/(2*sw2)*2*(pm - c)*sth/D,
            "30": mz_^2/(sv - mz_^2)*(rs/mw_)*b*(pm - c)/sqrt(2), "31": mz_^2/(sv - mz_^2)*b*sth}[name]
rs, c = 500, -1/2
print("e-  W+ W-   eq.      feynsage |M|^2 / Xianyu^2: first line   second line")
for pm in (1, -1):
    for hand, h3, h4, eq in (("L", 0, pm, "25"), ("L", -pm, 0, "25"), ("L", pm, pm, "26"), ("L", pm, -pm, "27"),
                             ("R", 0, pm, "30"), ("R", -pm, 0, "30"), ("R", pm, pm, "31"), ("R", pm, -pm, "32")):
        v = hel({"L": -1, "R": +1}[hand], h3, h4, rs, c)
        if eq == "32":
            print("%s   %+d  %+d   (21.32)  |M|^2 = %.1e" % (hand, h3, h4, float(v)))
            continue
        cols = [X(eq + "a", rs, c, pm), X(eq + "b", rs, c, pm)] if eq in ("25", "26", "27") else [X(eq, rs, c, pm)]
        print("%s   %+d  %+d   (21.%s)  " % (hand, h3, h4, eq) + "   ".join("%.12f" % float(v/x^2) for x in cols))
""", {}),
("c21_algebra", """
# only the printed formulas: does each line equal the next one?  (no feynsage, no diagrams)
Ev, pv, cv, wv = var("E p c w")
mwv = sqrt(Ev^2 - pv^2); mzv = mwv^2/(1 - wv); sv = 4*Ev^2; bv = pv/Ev; snv = sqrt(1 - cv^2); Dv = 1 + bv^2 + 2*bv*cv
uv = mwv^2 - 2*Ev^2 - 2*Ev*pv*cv
brv = mzv/(sv*(sv - mzv)) - 1/(2*wv)/(sv - mzv)
first24 = brv*(-4*Ev*pv*(Ev^2 + pv^2)/mwv^2 + 16*Ev^3*pv/mwv^2)*snv + 1/(2*wv)/uv*2*Ev*(-3*Ev^2*pv + pv^3 - 2*Ev^3*cv)*snv/mwv^2
A = mzv/(sv - mzv)*bv*(3 - bv^2); B = 1/(2*wv)*((2/Dv - sv/(sv - mzv))*bv*(3 - bv^2) + 4*cv/Dv)
second24 = lambda sA, sB: -sv/(4*mwv^2)*(sA*A + sB*B)*snv
zero = lambda x: bool(x.simplify_full() == 0)
first27 = lambda pm: -1/(2*wv)/uv*2*Ev^2*(-pm + cv)*snv
second27 = lambda pm: 1/(2*wv)*2*(pm - cv)*snv/Dv
tt = var("tt")                       # m_W/E -> 0
print("(21.24) first line = second line as printed, {A - B}:", zero(first24 - second24(1, -1)))
print("(21.24) first line = second line with {-A - B}:      ", zero(first24 - second24(-1, -1)))
print("(21.27) first = second:", zero(first27(1) - second27(1)), "   first = -(second):", zero(first27(1) + second27(1)))
lim = (first24/snv).subs({pv: Ev*sqrt(1 - tt^2)}).simplify_full().limit(tt=0).factor()
print("high-energy limit of the (21.24) first line, divided by sin(theta):", lim)
print("Peskin (21.107) bracket 1/(2c^2) - 1/(4c^2 s^2) + 1/(2s^2):", (1/(2*(1 - wv)) - 1/(4*(1 - wv)*wv) + 1/(2*wv)).factor())
""", {}),
("c21_fig", """
import numpy as np, matplotlib.pyplot as plt
# Peskin Fig. 21.10: e-_L e+_R at E_cm = 1000 GeV, theta = angle of the W- (cos = -cos_theta of the W+),
# curves (h-, h+); dsigma/dcos in units of R = 4 pi alpha^2/(3 E_cm^2) is (3/8) beta |M|^2 with e = 1
rs = 1000; beta = sqrt(1 - 4*mw_^2/rs^2).n()
pairs = [(hm, hp) for hm in (1, -1, 0) for hp in (1, -1, 0)]
cfg = [{"e-": he, "e+": -he, "W+": hp, "W-": hm} for he in (-1, 1) for hm, hp in pairs]
tab = W.helicity_table(cfg)                          # all 18 states, the traces in parallel
f = [fast_callable(SR(x.subs(num_W).subs({W.sqrt_s: rs})), vars=[W.cos_theta], domain=CDF) for x in tab]
cPs = [k/100 for k in range(-99, 100, 2)]
ds = lambda k, cP: 3/8*float(beta)*f[k](-cP).real()
cur = {pq: [ds(k, x) for x in cPs] for k, pq in enumerate(pairs)}
tot = [sum(cur[pq][i] for pq in pairs) for i in range(len(cPs))]
totR = [sum(ds(9 + k, x) for k in range(9)) for x in cPs]
fig, ax = plt.subplots(figsize=(5.2, 3.8))
ax.semilogy(cPs, tot, "k", lw=1.5, label="total")
for pq, lab in (((0, 0), "(0,0)"), ((1, -1), "(+,-)"), ((-1, 1), "(-,+)")):
    ax.semilogy(cPs, cur[pq], lw=1.1, label=lab)
ax.semilogy(cPs, [a + b for a, b in zip(cur[(-1, 0)], cur[(0, 1)])], lw=1.1, label="(-,0)+(0,+)")
ax.semilogy(cPs, [a + b for a, b in zip(cur[(1, 0)], cur[(0, -1)])], lw=1.1, label="(+,0)+(0,-)")
ax.semilogy(cPs, totR, "k--", lw=1, label=r"$e^-_R e^+_L$ total")
ax.set_ylim(0.01, 200); ax.set_xlim(-1, 1); ax.set_xlabel(r"$\\cos\\theta$"); ax.set_ylabel(r"$d\\sigma/d\\cos\\theta$ (units of R)")
ax.legend(fontsize=7, ncol=2); fig.tight_layout()
i0 = cPs.index(1/100)
print("e-_L / e-_R at cos = 0.01:", round(tot[i0]/totR[i0], 1), "   (+,+) and (-,-) largest for |cos| < 0.9:",
      "%.1e" % max(max(cur[(1, 1)][i], cur[(-1, -1)][i]) for i in range(len(cPs)) if abs(cPs[i]) < 0.9))
""", {"fig": "fig"}),
("c21_growth", """
# one W longitudinal, one transverse (21.2a): single diagrams grow like sqrt(s), the sum falls like 1/sqrt(s)
names = [str(x.propagator_fields()[0]) for x in W.diagrams]
print("|M|^2 for e-_L e+_R -> W+_0 W-_(+), cos = 0.3;  diagrams:", names)
for rs in (10^3, 10^4, 10^5):
    print("sqrt(s) = %6d GeV: " % rs + "  ".join("%9.3e" % float(hel(-1, 0, 1, rs, 3/10, diagrams=[k])) for k in range(1, 5))
          + "   sum: %9.3e" % float(hel(-1, 0, 1, rs, 3/10)))
""", {}),
]

CELLS += [
("p_qed", """
from feynsage.models import QED
P = process("e- e+ -> mu- mu+", model=QED(), massless=["e", "mu"])
P, P.squared().factor()
""", {}),
("p_hel", """
[(he, hm, P.squared(helicities={"e-": he, "e+": -he, "mu-": hm, "mu+": -hm}).factor()) for he in (1, -1) for hm in (1, -1)]
""", {}),
("p_sigma", """
P.dsigma_dcos().factor(), P.sigma(sqrt_s=10, unit="pb"), float(4*pi*(1/137.035999084)^2/(3*10^2)*0.3893793721e9)
""", {}),
("p_sm", """
Pz = process("e- e+ -> mu- mu+", widths={"Z": 2.4952})       # the full Standard Model, Z with its width
[(rs, round(Pz.sigma(sqrt_s=rs, unit="pb"), 3)) for rs in (20, 60, 91.1876, 120, 200)]
""", {}),
("p_decay", """
Hbb = process("H -> b b~")
G = Hbb.width()
print("Gamma(H -> b b~) divided by alpha m_h m_b^2 N_c (1 - 4 m_b^2/m_h^2)^(3/2)/(8 sin^2 m_W^2):",
      (G/(3*EL^2/(4*pi)*MH*MB^2*(1 - 4*MB^2/MH^2)^(3/2)/(8*SW^2*MW^2))).subs({MH: 125, MB: 5, MW: 80, SW: sqrt(0.23)}).n(digits=15))
for text in ("H -> b b~", "Z -> mu- mu+", "Z -> nu_e nu_e~", "W+ -> e+ nu_e", "t -> W+ b"):
    D = process(text)
    print("%-16s %.6f GeV" % (text, float(D.numeric(D.width()).n())))
""", {}),
("p_ww", """
W = process("e- e+ -> W+ W-", massless=["e"])
LL = W.squared(helicities={"e-": -1, "e+": +1, "W+": 0, "W-": 0})
W.numeric(LL, sqrt_s=500, cos_theta=-1/2, alpha=1/(4*pi)).n(digits=20)       # e = 1
""", {}),
("p_qcd", """
Q = process("u u~ -> g g", exclude_fields=["A", "Z", "H", "G0", "G+", "W+"], massless=["u"])
s, t, u = Q.s, Q.t, Q.u
from feynsage.models import GS
(Q.squared() - 32*GS^4/27*(t/u + u/t - 9*(t^2 + u^2)/(4*s^2))).subs(s=-t-u).simplify_rational()
""", {}),
]
