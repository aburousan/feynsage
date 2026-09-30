# Two textbook one-loop results of QED (Peskin and Schroeder), computed with feynsage and
# checked against the book's formulas:
#   1. vacuum polarisation Pi(q^2) (section 7.5): Ward identity, the Feynman-parameter
#      formula (7.90), and the leptonic running of alpha up to M_Z;
#   2. the anomalous magnetic moment F2(0) = alpha/(2 pi) (section 6.3), with U and F read
#      off the vertex graph, the Dirac algebra in FORM and the integrals done exactly.
import sys, time; sys.path.insert(0, '.')
from feynsage import *
import mpmath

ok = True
t0 = time.time()
# ------------------------------------------------------------------ 1. vacuum polarisation
T = form.dirac_trace(['mu', 'l + m', 'nu', 'l + q + m'], vectors=['l', 'q'], dim='D')
num = form.to_loop(T).replace('D', 'd')
res = loop(num, ["l", "m"], ["l + q", "m"], kin={"q^2": "s"})
s = var('s')
Cg, Cqq = res["g^{mu nu}"], res["q^mu q^nu"]
ward = (Cg + s * Cqq).expand()
print("Ward identity  C_g + q^2 C_qq =", ward)
ok &= bool(ward == 0)
# Pi^{mu nu} = (q^2 g - q q) Pi(q^2) with Pi = (alpha/4pi) C_qq  (see the notebook for the factors)
res0 = loop(num, ["l", "m"], ["l + q", "m"], kin={"q^2": 0})
Cqq0 = res0["q^mu q^nu"]
m_ = var('m')

def Pi_hat_ours(sv, mv=1):
    """(Pi(s) - Pi(0)) / (alpha/4pi), from feynsage, finite part at mu = 1."""
    a = finite_part(Cqq, 1).subs({s: sv, m_: mv}).n()
    b = finite_part(Cqq0, 1).subs({m_: mv}).n()
    return complex(a - b)

sys.path.insert(0, 'tests')
from peskin_ref import pi_hat as Pi_hat_peskin          # Peskin (7.91) with mpmath, plain Python


worst = 0
for sv in (-7, -1, 1, 3, 3.9, 5, 12, 40):
    a, b = Pi_hat_ours(sv), Pi_hat_peskin(sv)
    worst = max(worst, abs(a - b) / max(1, abs(b)))
    print("   s = %5s   feynsage %s   Peskin %s" % (sv, a, b))
print("Pi_hat: feynsage vs Peskin (7.91), worst relative difference over 8 points: %.1e" % worst)
ok &= bool(worst < 1e-10)

alpha = 1 / 137.035999084
MZ = 91187.6                                            # MeV
leptons = {"e": 0.51099895, "mu": 105.6583755, "tau": 1776.86}
dalpha = sum(alpha / (4 * pi.n()) * Pi_hat_ours(MZ ** 2, ml).real for ml in leptons.values())
approx = sum(alpha / (3 * pi.n()) * (log(MZ ** 2 / ml ** 2) - 5 / 3) for ml in leptons.values())
print("leptonic Delta alpha(M_Z^2) at one loop: %.6f   (large-s formula (alpha/3pi)(ln s/m^2 - 5/3): %.6f)" % (dalpha, approx))
print("   so 1/alpha runs from 137.036 to %.2f at M_Z from leptons alone" % (1 / (alpha * 1 / (1 - dalpha))))
ok &= bool(abs(dalpha - approx) < 2e-4)

# ------------------------------------------------------------------ 2. F2(0) = alpha/(2 pi)
Q2, m = var('Q2 m')
vert = graph("A-B, C-A:m, B-C:m", {"A": "p", "B": "-pp", "C": "pp - p"},
             kin={"p^2": "m^2", "pp^2": "m^2", "p.pp": "m^2 - Q2/2"})
U, F = vert.U(), vert.F()
print("vertex graph: U =", U, "  F =", F)
x2s, x3s = SR.var('x2'), SR.var('x3')
Delta = SR(str(F)).subs({SR.var('x1'): 1 - x2s - x3s}).expand()      # U = 1 on the simplex
vec = ['p', 'pp', 'k']
Dd, K2, Kp, Kpp = SR.var('D'), SR.var('K2'), SR.var('Kp'), SR.var('Kpp')
dot_ = function('dot')
table = {('p', 'p'): m ** 2, ('pp', 'pp'): m ** 2, ('p', 'pp'): m ** 2 - Q2 / 2, ('k', 'k'): K2, ('k', 'p'): Kp, ('k', 'pp'): Kpp}

def tr(factors):
    e = form.dirac_trace(factors, vectors=vec, dim='D')
    return e.substitute_function(dot_, lambda a, b: table.get((str(a), str(b)), table.get((str(b), str(a))))).expand()

# Gordon: Gamma^mu = A gamma^mu + B (p + pp)^mu/(2m) with F2 = -B; project with X = gamma_mu and X = (p+pp)_mu
M = matrix(SR, [[tr(['pp + m', 'mu', 'p + m', 'mu']), tr(['pp + m', 'p + pp', 'p + m']) / (2 * m)],
                [tr(['pp + m', 'p + pp', 'p + m']), tr(['pp + m', 'p + m']) * (4 * m ** 2 - Q2) / (2 * m)]])
lp = 'k - x2*p - x3*pp'                     # l = k - x2 p - x3 pp: the shift that makes F/U appear
NG = tr(['pp + m', 'nu', 'pp + %s + m' % lp, 'mu', 'p + %s + m' % lp, 'nu', 'p + m', 'mu'])
NV = tr(['pp + m', 'nu', 'pp + %s + m' % lp, 'p + pp', 'p + %s + m' % lp, 'nu', 'p + m'])

def symint(e):
    """Odd powers of k vanish; (k.a)(k.b) -> (a.b) k^2/D."""
    out = 0
    for t in (e.operands() if 'add' in str(e.operator()) else [e]):
        a, b = t.degree(Kp), t.degree(Kpp)
        c = t.subs({Kp: 1, Kpp: 1})
        if a + b == 0:
            out += c
        elif (a, b) in ((2, 0), (0, 2)):
            out += c * m ** 2 * K2 / Dd
        elif (a, b) == (1, 1):
            out += c * (m ** 2 - Q2 / 2) * K2 / Dd
    return out

AB = M.solve_right(vector(SR, [symint(NG), symint(NV)]))
F2num = (-AB[1]).expand()
eps_, Dl = SR.var('eps_'), SR.var('Dl')

def kint(j):                                 # Int d^dk/(i pi^(d/2)) (k^2)^j/(k^2 - Delta)^3
    d = 4 - 2 * eps_
    return (-1) ** (3 + j) * gamma(j + d / 2) * gamma(3 - j - d / 2) / (gamma(d / 2) * gamma(3)) * Dl ** (d / 2 + j - 3)

integrand = (F2num.coefficient(K2, 0) * kint(0) + F2num.coefficient(K2, 1) * kint(1)).subs({Dd: 4 - 2 * eps_})
ser = integrand.series(eps_, 1)
pole = ser.coefficient(eps_, -1).simplify_full()
fin = ser.coefficient(eps_, 0).subs({Dl: Delta})
fin0 = fin.series(Q2, 1).truncate().coefficient(Q2, 0).simplify_full()
print("F2: UV pole =", pole, "   integrand at q^2 = 0:", fin0)
assume(x2s > 0, x2s < 1, m > 0)
J = 2 * integrate(integrate(fin0, x3s, 0, 1 - x2s), x2s, 0, 1)
F2 = (J / (16 * pi ** 2) * 4 * pi).simplify_full()       # e^2/(16 pi^2) J with e^2 = 4 pi alpha
print("F2(0) =", F2, "* alpha")
ok &= bool(pole == 0) and bool(F2 == 1 / (2 * pi))
a_exp = 0.00115965218059
print("  alpha/(2 pi) = %.10f   measured a_e = %.14f (Fan et al. 2023): the one-loop term is %.2f%% above it"
      % (alpha / (2 * pi.n()), a_exp, 100 * (alpha / (2 * pi.n()) / a_exp - 1)))
print("time %.1fs" % (time.time() - t0))
print("ALL PESKIN CHECKS PASS" if ok else "SOME CHECK FAILED")
