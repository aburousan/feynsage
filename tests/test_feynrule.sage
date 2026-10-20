# Feynman rules from Lagrangians (feynsage.feynrule) against the textbook rules (Peskin and Schroeder) and against
# FeynCalc 10.2's FeynRule (tests/feyncalc_ref/feynrule.wl, its output in feynrule_fc.txt), compared number by
# number at random values of the Lorentz and colour indices with the SU(3) structure constants.
#     sage tests/test_feynrule.sage
import itertools, random, re, time, os
from feynsage import *
from feynsage.feynrule import d
from feynsage.qft.color import COLF, KRON, f_color
t_start = time.time()
FAIL = []
def check(name, ok):
    print(("ok   " if ok else "FAIL ") + name)
    if not ok:
        FAIL.append(name)
def same(x, y):
    return bool((SR(x) - SR(y)).expand() == 0)

p1, p2, p3, p4 = momenta("p1 p2 p3 p4")
mu, nu, rho, sg, al, be, mu1, nu1 = lorentz_indices("mu nu rho sg al be mu1 nu1")
a, b, c, dd, aa, a1, e9, x, y = color_indices("a b c dd aa a1 e9 x y")
i1, i2 = quark_colors("i1 i2")
lam, e, g, gs, m, CH = var('lam e g gs m CH')

phi = field("phi")
check("phi^4 with L = -lam/4! phi^4: -i lam (Peskin p. 95)", same(feynman_rule(-lam/24*phi()^4, [(phi, p1), (phi, p2), (phi, p3), (phi, p4)]), -I*lam))
check("the two-point function of L = (d phi)^2/2 - m^2 phi^2/2 is i (p^2 - m^2), p2 = -p1",
      same(feynman_rule(d(mu, phi())*d(mu, phi())/2 - m^2/2*phi()^2, [(phi, p1), (phi, p2)]), -I*m^2 - I*dot(p1, p2)))
psi = field("psi", "dirac"); A = field("A", "vector")
r = feynman_rule(-e*psi.bar()*gamma(mu)*psi()*A(mu), [(A, p1, nu), (psi, p2), (psi.bar, p3)])
check("QED, L = -e psibar gamma^mu psi A_mu: -i e gamma^nu (Peskin p. 124: -i Q|e| gamma^mu with Q|e| = e)", str(r) == str(-I*e*gamma(nu)))
q = field("q", "dirac", colour="fundamental"); G = field("G", "vector", colour="adjoint")
r = feynman_rule(I*q.bar(i1)*gamma(mu)*covariant_derivative(mu, q, G, gs, i=i1), [(G, p1, nu, a), (q, p2, i2), (q.bar, p3, i1)])
check("quark-gluon from i qbar D/ q: i g_s gamma^nu T^a (FeynCalc's FeynRule documentation)", str(r) == str(I*gs*T_color(a, i1, i2)*gamma(nu)))
YM = -field_strength(G, al, be, aa, gs)*field_strength(G, al, be, aa, gs)/4
r3 = feynman_rule(YM, [(G, p1, mu, a), (G, p2, nu, b), (G, p3, rho, c)])
ref3 = gs*f_color(a, b, c)*(metric(mu, nu)*(comp(p1, rho) - comp(p2, rho)) + metric(nu, rho)*(comp(p2, mu) - comp(p3, mu))
                            + metric(mu, rho)*(comp(p3, nu) - comp(p1, nu)))
check("three-gluon vertex g f^abc [g^{mu nu}(k - p)^rho + ...], all momenta incoming, from (16.6) (Peskin Fig. 16.1, p. 507)", same(r3, ref3))
gh = field("c", "ghost", colour="adjoint")
Lgh = -gh.bar(x)*d(mu, d(mu, gh(x)) + gs*f_color(x, y, c)*G(mu, y)*gh(c))
check("ghost-gluon vertex from (16.32): -g f^abc p^mu with p the momentum going out along the ghost line (Peskin Fig. 16.5, p. 515)",
      same(feynman_rule(Lgh, [(gh.bar, p1, a), (G, p2, mu, b), (gh, p3, c)]).subs({comp(p2, mu): -comp(p1, mu) - comp(p3, mu)}),
           gs*f_color(a, b, c)*comp(p1, mu)))
cphi = field("phi", "complex")
Dphi = d(mu, cphi()) + I*e*A(mu)*cphi(); Dphid = d(mu, cphi.bar()) - I*e*A(mu)*cphi.bar()
Ls = Dphid*Dphi - m^2*cphi.bar()*cphi()
check("scalar QED, D = d + i e A: -i e (p + p')^mu and the seagull 2 i e^2 g^{mu nu} (Peskin Problem 9.1, p. 312)",
      same(feynman_rule(Ls, [(cphi, p1), (cphi.bar, p2), (A, p3, mu)]), -I*e*(comp(p1, mu) - comp(p2, mu)))
      and same(feynman_rule(Ls, [(cphi, p1), (cphi.bar, p2), (A, p3, mu), (A, p4, nu)]), 2*I*e^2*metric(mu, nu)))

# ---------------------------------------------------------------------------- against FeynCalc, number by number
H = field("H")
heft = -CH/4*field_strength(G, mu1, nu1, a1, gs)*field_strength(G, mu1, nu1, a1, gs)*H()
ours = {"G4": feynman_rule(YM, [(G, p1, mu, a), (G, p2, nu, b), (G, p3, rho, c), (G, p4, sg, dd)]),
        "HGG": feynman_rule(heft, [(G, p1, mu, a), (G, p2, nu, b), (H, p3)]),
        "HGGG": feynman_rule(heft, [(G, p1, mu, a), (G, p2, nu, b), (G, p3, rho, c), (H, p4)])}
Kron = color_delta
fc = {}
for line in open(os.path.join("tests", "feyncalc_ref", "feynrule_fc.txt")):
    k, v = line.split(": ", 1)
    v = v.strip().replace('FCGV["sun441"]', 'e9').replace("[", "(").replace("]", ")").replace("fcolor", "f_color").replace("SD(", "Kron(")
    v = re.sub(r"(\w+) \. (\w+)", r"Dot(\1, \2)", v).replace(" d,", " dd,").replace("(d,", "(dd,").replace(", d)", ", dd)")
    fc[k] = SR(eval(preparse(v)))
f = {}
for (u, v, w), val in {(1, 2, 3): 1, (1, 4, 7): 1/2, (1, 5, 6): -1/2, (2, 4, 6): 1/2, (2, 5, 7): 1/2, (3, 4, 5): 1/2, (3, 6, 7): -1/2,
                       (4, 5, 8): sqrt(3)/2, (6, 7, 8): sqrt(3)/2}.items():
    for perm, sgn in [((u, v, w), 1), ((v, w, u), 1), ((w, u, v), 1), ((v, u, w), -1), ((u, w, v), -1), ((w, v, u), -1)]:
        f[perm] = sgn*val
gmn = diagonal_matrix([1, -1, -1, -1])
rnd = random.Random(int(5))
P = {str(k): vector(QQ, [rnd.randint(-9, 9) for _ in range(4)]) for k in (p1, p2, p3, p4)}
def value(ex, L, C):
    ex = SR(ex).expand()
    terms = ex.operands() if ex.operator() is not None and "add" in str(ex.operator()) else [ex]
    tot = 0
    for t in terms:
        fs = t.operands() if t.operator() is not None and "mul" in str(t.operator()) else [t]
        dummies = sorted({str(z) for F in fs if F.operator() in (COLF, KRON) for z in F.operands()} - set(C))
        for dv in itertools.product(range(1, 9), repeat=len(dummies)):
            cv = dict(C); cv.update(dict(zip(dummies, dv)))
            v = 1
            for F in fs:
                op = F.operator(); args = [str(z) for z in F.operands()] if op is not None else []
                if op is COLF: v *= f.get(tuple(cv[z] for z in args), 0)
                elif op is KRON: v *= 1 if cv[args[0]] == cv[args[1]] else 0
                elif str(op) == "Metric": v *= gmn[L[args[0]], L[args[1]]]
                elif str(op) == "Comp": v *= P[args[0]][L[args[1]]]
                elif str(op) == "Dot": v *= P[args[0]] * gmn * P[args[1]]
                else: v *= F.subs(gs=3/7, CH=5/11)
            tot += v
    return tot
# Peskin Fig. 16.1 (p. 507): legs (a, mu), (b, nu), (c, rho), (d, sigma)
pk4 = -I*gs^2*(f_color(a, b, e9)*f_color(c, dd, e9)*(metric(mu, rho)*metric(nu, sg) - metric(mu, sg)*metric(nu, rho))
              + f_color(a, c, e9)*f_color(b, dd, e9)*(metric(mu, nu)*metric(rho, sg) - metric(mu, sg)*metric(nu, rho))
              + f_color(a, dd, e9)*f_color(b, c, e9)*(metric(mu, nu)*metric(rho, sg) - metric(mu, rho)*metric(nu, sg)))
fc["G4peskin"] = pk4; ours["G4peskin"] = ours["G4"]
colour_sets = {2: [(1, 1), (4, 4), (2, 3)], 3: [(1, 2, 3), (4, 5, 8), (1, 4, 7), (3, 6, 7), (1, 1, 2)],
               4: [(1, 2, 1, 2), (1, 2, 2, 1), (4, 5, 4, 5), (1, 4, 2, 6), (1, 1, 2, 2), (3, 8, 3, 8)]}
for k, label in [("G4", "four-gluon vertex from -F^2/4"), ("HGG", "H g g from -C_H F^2 H/4 (HEFT)"), ("HGGG", "H g g g (HEFT)"),
                 ("G4peskin", "four-gluon vertex from -F^2/4")]:
    n = {"G4": 4, "HGG": 2, "HGGG": 3, "G4peskin": 4}[k]
    ok, nonzero, total = True, 0, 0
    for lv in itertools.product(range(4), repeat=n):              # every choice of the Lorentz indices
        L = {str([mu, nu, rho, sg][i]): lv[i] for i in range(n)}
        for cvals in colour_sets[n]:
            C = {str([a, b, c, dd][i]): cvals[i] for i in range(n)}
            u_, v_ = value(ours[k], L, C), value(fc[k], L, C)
            ok = ok and bool(SR(u_ - v_) == 0)
            nonzero += int(bool(SR(u_) != 0)); total += 1
    ref = "Peskin Fig. 16.1 (p. 507)" if k == "G4peskin" else "FeynCalc's FeynRule"
    check("%s = %s for every Lorentz index and %d colour choices (%d of %d values non-zero)"
          % (label, ref, len(colour_sets[n]), nonzero, total), ok and nonzero > 0)

print("time %.0f s" % (time.time() - t_start))
print("ALL FEYNMAN-RULE CHECKS PASS" if not FAIL else "FAILED: %s" % FAIL)
