# Sector symmetries: symmetries that hold inside one sector only, as LiteRed and Kira find them.
# Two-loop sunset, Euclidean, p^2 = pp, with the numerator lines k + p and l.
#  - three different masses: 7 masters (4 in the sunset sector, 3 tadpole products)
#  - three equal masses: the three tadpole products are one integral and the sunset sector has
#    two masters, so 3 in all (without sector symmetries the reducer kept 6)
# The reductions of the equal-mass case are then checked against the Feynman-parameter integral
# at D = 2.6, where every integral converges.
import sys, math; sys.path.insert(0, '.')
from feynsage import *
from scipy.integrate import dblquad, quad

def sunset(m1, m2, m3):
    return family([("k", m1), ("k - l", m2), ("l + p", m3), "k + p", "l"], kin={"p^2": "pp"}, euclidean=True)

targets = ["F(2,1,1,0,0)", "F(1,2,1,0,0)", "F(1,1,2,0,0)", "F(2,2,1,0,0)", "F(1,1,1,-1,0)"]
n3 = len(ibp_reduce(sunset("1", "2", "3"), targets).masters)
eq = sunset("1", "1", "1")
red = ibp_reduce(eq, targets[:4] + ["F(2,2,2,0,0)"])
print("three masses:", n3, "masters   equal masses:", len(red.masters), red.masters)
assert n3 == 7 and len(red.masters) == 3

U, F = eq.UF()
Dv, ppv, L = 2.6, 7/3, 2
X = [SR.var('x%d' % i) for i in range(1, 6)]
Uf = fast_callable(SR(str(U)), vars=X)
Ff = fast_callable(SR(str(F)).subs(pp=ppv), vars=X)

def value(a):
    lines = [i for i in range(5) if a[i] > 0]; A = sum(a)
    pre = math.gamma(A - L*Dv/2) / math.prod(math.gamma(a[i]) for i in lines)
    def f(*y):
        x = [0.0]*5
        for i, v in zip(lines, y): x[i] = v
        return math.prod(x[i]**(a[i] - 1) for i in lines) * Uf(*x)**(A - (L + 1)*Dv/2) / Ff(*x)**(A - L*Dv/2)
    if len(lines) == 2:
        return pre * quad(lambda u: f(u, 1 - u), 0, 1, epsabs=1e-12, epsrel=1e-10, limit=200)[0]
    return pre * dblquad(lambda v, u: f(u, v, 1 - u - v), 0, 1, 0, lambda u: 1 - u, epsabs=1e-11, epsrel=1e-9)[0]

M = {m: value(m) for m in red.masters}
worst = 0
for a, row in red.table.items():
    rhs = sum(float(SR(str(c)).subs({SR.var('d'): Dv, SR.var('pp'): ppv})) * M[m] for m, c in row.items())
    lhs = value(a)
    worst = max(worst, abs(lhs - rhs) / abs(lhs))
    print(a, "direct %.10f  from masters %.10f" % (lhs, rhs))
print("largest relative difference: %.1e" % worst)
assert worst < 1e-6
print("SECTOR SYMMETRIES OK")
