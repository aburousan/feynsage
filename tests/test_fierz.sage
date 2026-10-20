# Products of two spinor chains in four dimensions (feynsage.fierz): Fierz rearrangement, Clifford decomposition and
# Sirlin's relations, checked with explicit 4x4 matrices (Peskin's chiral basis, g5 = i g0 g1 g2 g3, eps^0123 = +1)
# and random complex spinors, summing every repeated index with the metric; and against FeynCalc's SirlinSimplify.
#     sage tests/test_fierz.sage
import itertools, random, time
from feynsage import *
from feynsage.fierz import fierz, sirlin, clifford_decompose
from feynsage.qft.objects import DiracExpr
t_start = time.time()
FAIL = []
def check(name, ok):
    print(("ok   " if ok else "FAIL ") + name)
    if not ok:
        FAIL.append(name)

# explicit matrices
s0 = identity_matrix(QQbar, 2); sx = matrix(QQbar, [[0, 1], [1, 0]]); sy = matrix(QQbar, [[0, -I], [I, 0]]); sz = matrix(QQbar, [[1, 0], [0, -1]])
Z = zero_matrix(QQbar, 2)
G = [block_matrix([[Z, s0], [s0, Z]])] + [block_matrix([[Z, s], [-s, Z]]) for s in (sx, sy, sz)]
G5 = I * G[0] * G[1] * G[2] * G[3]
g = diagonal_matrix([1, -1, -1, -1])
def eps(*ix):
    if len(set(ix)) < 4: return 0
    return Permutation([i + 1 for i in ix]).sign()          # eps^{0123} = +1, upper indices
rnd = random.Random(int(11))
cplx = lambda: QQbar(rnd.randint(-5, 5) + I*rnd.randint(-5, 5))
spin = {}
def spinor(s):
    key = (s.kind, str(s.p))
    if key not in spin:
        spin[key] = vector(QQbar, [cplx() for _ in range(4)])
    return spin[key]

def factor_matrix(F, iv):
    M = zero_matrix(QQbar, 4)
    for atom, c in F.items:
        cc = QQbar(num_coef(c, iv))
        if atom[0] == 'i': M += cc * G[iv[str(atom[1])]]
        elif atom[0] == '5': M += cc * G5
        elif atom[0] == '1': M += cc * identity_matrix(QQbar, 4)
        else: raise ValueError(atom)
    return M

def num_coef(c, iv):
    """a coefficient with Metric / Epsilon of the indices in iv, upper indices"""
    c = SR(c)
    w = [SR.wild(k) for k in range(4)]
    for F in c.find(SR(metric(*lorentz_indices("zza zzb"))).operator()(w[0], w[1])):
        a_, b_ = [str(x) for x in F.operands()]
        c = c.subs({F: g[iv[a_], iv[b_]]})
    for F in c.find(SR(epsilon(*lorentz_indices("zzc zzd zze zzf"))).operator()(w[0], w[1], w[2], w[3])):
        c = c.subs({F: eps(*[iv[str(x)] for x in F.operands()])})
    return c

def indices_of(term):
    c, o, cl = term
    names = [str(v) for v in SR(c).variables() if T.is_index(v)]
    for ch in cl:
        for F in ch.factors:
            for atom, cc in F.items:
                if atom[0] == 'i': names.append(str(atom[1]))
                names += [str(v) for v in SR(cc).variables() if T.is_index(v)]
    return names

def value(M, free):
    """numerical value of a product of two bilinears; free = {index name: value}.  Every monomial of a
    coefficient is summed over its own repeated indices (those it shares with the chains included)."""
    tot = 0
    for c, o, cl in M.terms:
        ce = SR(c).expand()
        monos = ce.operands() if ce.operator() is not None and "add" in str(ce.operator()) else [ce]
        for mono in monos:
            term = (mono, o, cl)
            names = indices_of(term)
            summed = sorted({n for n in names if n not in free})
            for vals in itertools.product(range(4), repeat=len(summed)):
                iv = dict(free); iv.update(zip(summed, vals))
                w = prod(g[v, v] for v in vals)        # one index up, one down: g_{aa} per summed index
                x = QQbar(num_coef(mono, iv)) * w
                for ch in cl:
                    Mx = identity_matrix(QQbar, 4)
                    for F in ch.factors:
                        Mx = Mx * factor_matrix(F, iv)
                    x *= spinor(ch.left) * Mx * spinor(ch.right)
                tot += x
    return tot

from feynsage.qft import tensors as T
p1, p2, p3, p4 = momenta("p1 p2 p3 p4"); mu, nu, rho, tau = lorentz_indices("mu nu rho tau")

def M2(GA, GB):
    return ubar(p1) * GA * u(p2) * ubar(p3) * GB * u(p4)

cases = [("S x S", one(), one(), {}), ("V x V", gamma(mu), gamma(mu), {}), ("(V-A) x (V-A)", gamma(mu)*PL(), gamma(mu)*PL(), {}),
         ("T x T", sigma(mu, nu), sigma(mu, nu), {}), ("g^mu P_L x g^nu P_R (free mu nu)", gamma(mu)*PL(), gamma(nu)*PR(), {"mu": 2, "nu": 3})]
for name, GA, GB, free in cases:
    lhs = M2(GA, GB); rhs = fierz(lhs)
    check("Fierz %s: both sides equal for random spinors (explicit matrices)" % name, value(lhs, free) == value(rhs, free))
# the classic: (V-A) x (V-A) is Fierz invariant for anticommuting fields
lhs = M2(gamma(mu)*PL(), gamma(mu)*PL())
swapped = ubar(p1) * gamma(mu)*PL() * u(p4) * ubar(p3) * gamma(mu)*PL() * u(p2)
check("(V-A) x (V-A): [1 g^mu P_L 2][3 g_mu P_L 4] = -[1 g^mu P_L 4][3 g_mu P_L 2] for commuting spinors", value(lhs, {}) == -value(swapped, {}))
# Sirlin, FeynCalc's example: [3 g^mu g^rho g^nu P_L 1][4 g_mu g_tau g_nu P_L 2] = 4 g^{rho tau} [3 g^a P_L 1][4 g_a P_L 2]
al = lorentz_indices("al")
lhs = ubar(p3)*gamma(mu)*gamma(rho)*gamma(nu)*PL()*u(p1) * ubar(p4)*gamma(mu)*gamma(tau)*gamma(nu)*PL()*u(p2)
ref = 4*metric(rho, tau) * (ubar(p3)*gamma(al)*PL()*u(p1) * ubar(p4)*gamma(al)*PL()*u(p2))
ok = all(value(sirlin(lhs), {"rho": r, "tau": t_}) == value(ref, {"rho": r, "tau": t_}) == value(lhs, {"rho": r, "tau": t_})
         for r, t_ in [(0, 0), (1, 1), (2, 3), (3, 3)])
check("SirlinSimplify documentation example: = 4 g^{rho tau} [g^a P_L][g_a P_L], and sirlin() keeps the value", ok)
s1 = lorentz_indices("s1")
compact = metric(rho, tau)*(ubar(p3)*gamma(s1)*u(p1)*ubar(p4)*gamma(s1)*u(p2) - ubar(p3)*gamma(s1)*u(p1)*ubar(p4)*gamma(s1)*gamma5()*u(p2)
                            - ubar(p3)*gamma(s1)*gamma5()*u(p1)*ubar(p4)*gamma(s1)*u(p2) + ubar(p3)*gamma(s1)*gamma5()*u(p1)*ubar(p4)*gamma(s1)*gamma5()*u(p2))
r = sirlin(lhs)
check("sirlin() gives that compact form: g^{rho tau} [VV - VA - AV + AA] = 4 g^{rho tau} [g^a P_L][g_a P_L] (4 terms)",
      len(r.terms) == 4 and all(value(r - compact, {"rho": a_, "tau": b_}) == 0 for a_, b_ in [(0, 0), (1, 1), (1, 2)]))

print("time %.0f s" % (time.time() - t_start))
print("ALL FIERZ CHECKS PASS" if not FAIL else "FAILED: %s" % FAIL)
