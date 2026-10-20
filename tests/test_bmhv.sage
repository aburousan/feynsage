# BMHV traces (feynsage.bmhv) against explicit gamma matrices in d = 6, where the BMHV algebra exists exactly:
# G^mu = gamma^mu x 1 (mu = 0..3), G^4 = i g5 x s1, G^5 = i g5 x s2, gamma5 = g5 x 1 anticommutes with G^0..G^3 and
# commutes with G^4, G^5; traces divided by 2 (Tr 1 = 4).  Every result is evaluated at random index values 0..5 and
# random six-momenta, and compared; and against FeynCalc's DiracTrace with $BreitMaison = True.
#     sage tests/test_bmhv.sage
import itertools, random, time
from feynsage import *
from feynsage.bmhv import bmhv_trace, to_d, MetricBar, MetricHat, CompBar, CompHat, DotBar, DotHat
from feynsage.qft import tensors as T
t_start = time.time()
FAIL = []
def check(name, ok):
    print(("ok   " if ok else "FAIL ") + name)
    if not ok:
        FAIL.append(name)
s0 = identity_matrix(QQbar, 2); Z = zero_matrix(QQbar, 2)
sx, sy, sz = matrix(QQbar, [[0, 1], [1, 0]]), matrix(QQbar, [[0, -I], [I, 0]]), matrix(QQbar, [[1, 0], [0, -1]])
g4 = [block_matrix([[Z, s0], [s0, Z]])] + [block_matrix([[Z, s], [-s, Z]]) for s in (sx, sy, sz)]
g5_4 = I * g4[0] * g4[1] * g4[2] * g4[3]
G = [g.tensor_product(s0) for g in g4] + [(I * g5_4).tensor_product(sx), (I * g5_4).tensor_product(sy)]
G5 = g5_4.tensor_product(s0)
met = diagonal_matrix([1, -1, -1, -1, -1, -1])
assert all(G[a]*G[b] + G[b]*G[a] == 2*met[a, b]*identity_matrix(QQbar, 8) for a in range(6) for b in range(6))
assert all(G5*G[a] == (-1 if a < 4 else 1)*G[a]*G5 for a in range(6))
rnd = random.Random(int(17))
mu, nu, rho, sg, tau, ka = lorentz_indices("mu nu rho sg tau ka")
p, q, l, k = momenta("p q l k")
P = {"p": vector(QQ, [rnd.randint(-5, 5) for _ in range(4)] + [0, 0]),        # p, q four-dimensional
     "q": vector(QQ, [rnd.randint(-5, 5) for _ in range(4)] + [0, 0]),
     "l": vector(QQ, [rnd.randint(-5, 5) for _ in range(6)]), "k": vector(QQ, [rnd.randint(-5, 5) for _ in range(6)])}
low = lambda v: vector(QQ, [met[a, a]*v[a] for a in range(6)])
def eps4(*ix):
    if any(x > 3 for x in ix) or len(set(ix)) < 4: return 0
    return Permutation([x + 1 for x in ix]).sign()
def explicit(atoms, iv):
    M = identity_matrix(QQbar, 8)
    for a in atoms:
        if a[0] == 'i': M = M * G[iv[a[1]]]
        elif a[0] == 'p': v = low(P[a[1]]); M = M * sum(v[x]*G[x] for x in range(6))
        else: M = M * G5
    return M.trace() / 2
def num(e, iv):
    e = SR(e).subs(d=6)
    w = [SR.wild(n) for n in range(4)]
    sub = {}
    opM, opC, opD, opE = [SR(x).operator() for x in (metric(mu, nu), comp(p, mu), dot(p, q), epsilon(mu, nu, rho, sg))]
    def val(x):
        s = str(x)
        return ('i', iv[s]) if s in iv else ('v', P[s])
    for F in e.find(MetricBar(w[0], w[1])): a, b = [iv[str(x)] for x in F.operands()]; sub[F] = met[a, b] if a < 4 and b < 4 else 0
    for F in e.find(MetricHat(w[0], w[1])): a, b = [iv[str(x)] for x in F.operands()]; sub[F] = met[a, b] if a > 3 and b > 3 else 0
    for F in e.find(opM(w[0], w[1])): a, b = [iv[str(x)] for x in F.operands()]; sub[F] = met[a, b]
    for F in e.find(CompBar(w[0], w[1])): v, a = F.operands(); sub[F] = P[str(v)][iv[str(a)]] if iv[str(a)] < 4 else 0
    for F in e.find(CompHat(w[0], w[1])): v, a = F.operands(); sub[F] = P[str(v)][iv[str(a)]] if iv[str(a)] > 3 else 0
    for F in e.find(opC(w[0], w[1])): v, a = F.operands(); sub[F] = P[str(v)][iv[str(a)]]
    for F in e.find(DotBar(w[0], w[1])): u, v = F.operands(); sub[F] = sum(met[x, x]*P[str(u)][x]*P[str(v)][x] for x in range(4))
    for F in e.find(DotHat(w[0], w[1])): u, v = F.operands(); sub[F] = sum(met[x, x]*P[str(u)][x]*P[str(v)][x] for x in range(4, 6))
    for F in e.find(opD(w[0], w[1])): u, v = F.operands(); sub[F] = P[str(u)] * met * P[str(v)]
    for F in e.find(opE(w[0], w[1], w[2], w[3])):
        args = [val(x) for x in F.operands()]
        tot = 0
        for combo in itertools.product(range(4), repeat=4):
            ok_, fac = True, 1
            for (kind, v), x in zip(args, combo):
                if kind == 'i':
                    if v != x: ok_ = False; break
                else:
                    fac *= low(v)[x]
            if ok_: tot += eps4(*combo) * fac
        sub[F] = tot
    return QQbar(e.subs(sub))
def agrees(G_, atoms, free_names, four_dim=(p, q), trials=40):
    r = bmhv_trace(G_, four_dim=four_dim)
    for _ in range(trials):
        iv = {n: rnd.randint(0, 5) for n in free_names}
        rep = {}
        for a in atoms:                                   # an index written twice: summed with the metric
            pass
        names = [a[1] for a in atoms if a[0] == 'i']
        summed = sorted({n for n in names if names.count(n) == 2})
        tot = 0
        for vals in itertools.product(range(6), repeat=len(summed)):
            iv2 = dict(iv); iv2.update(zip(summed, vals))
            tot += prod(met[v, v] for v in vals) * explicit(atoms, iv2)
        if num(r, iv) != tot:
            return False
    return True
cases = [("g^mu g^nu g^rho g^sg g5", gamma(mu)*gamma(nu)*gamma(rho)*gamma(sg)*gamma5(), [('i', 'mu'), ('i', 'nu'), ('i', 'rho'), ('i', 'sg'), ('5',)], ["mu", "nu", "rho", "sg"]),
         ("six gammas and g5", gamma(mu)*gamma(nu)*gamma(rho)*gamma(sg)*gamma(tau)*gamma(ka)*gamma5(),
          [('i', x) for x in ("mu", "nu", "rho", "sg", "tau", "ka")] + [('5',)], ["mu", "nu", "rho", "sg", "tau", "ka"]),
         ("g5 in the middle: g^mu g^nu g5 g^rho g^sg", gamma(mu)*gamma(nu)*gamma5()*gamma(rho)*gamma(sg), [('i', 'mu'), ('i', 'nu'), ('5',), ('i', 'rho'), ('i', 'sg')], ["mu", "nu", "rho", "sg"]),
         ("two g5: g5 g^mu g^nu g5 g^rho g^sg", gamma5()*gamma(mu)*gamma(nu)*gamma5()*gamma(rho)*gamma(sg), [('5',), ('i', 'mu'), ('i', 'nu'), ('5',), ('i', 'rho'), ('i', 'sg')], ["mu", "nu", "rho", "sg"]),
         ("d-dimensional loop momentum: l/ g^mu k/ g^nu l/ g^rho g5", slash(l)*gamma(mu)*slash(k)*gamma(nu)*slash(l)*gamma(rho)*gamma5(),
          [('p', 'l'), ('i', 'mu'), ('p', 'k'), ('i', 'nu'), ('p', 'l'), ('i', 'rho'), ('5',)], ["mu", "nu", "rho"]),
         ("a summed index: g^mu g^nu g^rho g^sg g_mu g5", gamma(mu)*gamma(nu)*gamma(rho)*gamma(sg)*gamma(mu)*gamma5(),
          [('i', 'mu'), ('i', 'nu'), ('i', 'rho'), ('i', 'sg'), ('i', 'mu'), ('5',)], ["nu", "rho", "sg"]),
         ("no g5, hats from a d-dim momentum: l/ g^mu l/ g^nu", slash(l)*gamma(mu)*slash(l)*gamma(nu), [('p', 'l'), ('i', 'mu'), ('p', 'l'), ('i', 'nu')], ["mu", "nu"])]
for name, G_, atoms, free in cases:
    check("BMHV trace %s = explicit 8x8 gammas in d = 6" % name, agrees(G_, atoms, free))
# FeynCalc 10.2, $BreitMaison = True: DiracTrace[GAD[mu, nu, rho, sig] . GA5] = -4 I LC[mu, nu, rho, sig]
check("= FeynCalc: Tr[g^mu g^nu g^rho g^sg g5] = -4 i eps^{mu nu rho sg} (4-dimensional epsilon)",
      bool((bmhv_trace(gamma(mu)*gamma(nu)*gamma(rho)*gamma(sg)*gamma5()) + 4*I*epsilon(mu, nu, rho, sg)).expand() == 0))
# FeynCalc's six-gamma and momentum traces ($BreitMaison = True: LC is the 4-dimensional epsilon, MTD, FVD, SPD
# d-dimensional), against the same explicit matrices
fc_T6 = (-4*I)*(epsilon(nu, rho, sg, tau)*metric(ka, mu) - epsilon(mu, rho, sg, tau)*metric(ka, nu) + epsilon(mu, nu, sg, tau)*metric(ka, rho)
  - epsilon(mu, nu, rho, tau)*metric(ka, sg) + epsilon(mu, nu, rho, sg)*metric(ka, tau) - epsilon(ka, rho, sg, tau)*metric(mu, nu)
  + epsilon(ka, nu, sg, tau)*metric(mu, rho) - epsilon(ka, nu, rho, tau)*metric(mu, sg) + epsilon(ka, nu, rho, sg)*metric(mu, tau)
  - epsilon(ka, mu, sg, tau)*metric(nu, rho) + epsilon(ka, mu, rho, tau)*metric(nu, sg) - epsilon(ka, mu, rho, sg)*metric(nu, tau)
  - epsilon(ka, mu, nu, tau)*metric(rho, sg) + epsilon(ka, mu, nu, sg)*metric(rho, tau) - epsilon(ka, mu, nu, rho)*metric(sg, tau))
def fc_agrees(fc, atoms, free_names, trials=40):
    for _ in range(trials):
        iv = {n: rnd.randint(0, 5) for n in free_names}
        if num(fc, iv) != explicit(atoms, iv):
            return False
    return True
check("FeynCalc's BMHV six-gamma trace (LC x MTD) = explicit 8x8 gammas in d = 6, and = feynsage after to_d()",
      fc_agrees(fc_T6, [('i', x) for x in ("mu", "nu", "rho", "sg", "tau", "ka")] + [('5',)], ["mu", "nu", "rho", "sg", "tau", "ka"])
      and all(num(to_d(bmhv_trace(gamma(mu)*gamma(nu)*gamma(rho)*gamma(sg)*gamma(tau)*gamma(ka)*gamma5())) - fc_T6,
                  {n: rnd.randint(0, 5) for n in ["mu", "nu", "rho", "sg", "tau", "ka"]}) == 0 for _ in range(30)))

# d = 8 (16x16 matrices), for the d-dependence that summed indices bring: hats i g5 x c_k with four anticommuting c_k
c8 = [sx.tensor_product(s0), sy.tensor_product(s0), sz.tensor_product(sx), sz.tensor_product(sy)]
G8 = [g.tensor_product(identity_matrix(QQbar, 4)) for g in g4] + [(I*g5_4).tensor_product(c) for c in c8]
G58 = g5_4.tensor_product(identity_matrix(QQbar, 4)); met8 = diagonal_matrix([1, -1, -1, -1, -1, -1, -1, -1])
def explicit8(atoms, iv):
    M = identity_matrix(QQbar, 16)
    for a in atoms:
        M = M * (G8[iv[a[1]]] if a[0] == 'i' else G58)
    return M.trace() / 4
def num8(e, iv):
    e = SR(e).subs(d=8)
    w = [SR.wild(n) for n in range(4)]
    opM, opE = [SR(x).operator() for x in (metric(mu, nu), epsilon(mu, nu, rho, sg))]
    sub = {}
    for F in e.find(MetricBar(w[0], w[1])): a, b = [iv[str(x)] for x in F.operands()]; sub[F] = met8[a, b] if a < 4 and b < 4 else 0
    for F in e.find(MetricHat(w[0], w[1])): a, b = [iv[str(x)] for x in F.operands()]; sub[F] = met8[a, b] if a > 3 and b > 3 else 0
    for F in e.find(opM(w[0], w[1])): a, b = [iv[str(x)] for x in F.operands()]; sub[F] = met8[a, b]
    for F in e.find(opE(w[0], w[1], w[2], w[3])): sub[F] = eps4(*[iv[str(x)] for x in F.operands()])
    return QQbar(e.subs(sub))
ok8 = True
for atoms, free in [([('i', 'mu'), ('i', 'nu'), ('i', 'rho'), ('i', 'sg'), ('i', 'mu'), ('5',)], ["nu", "rho", "sg"]),
                    ([('i', 'mu'), ('i', 'nu'), ('i', 'rho'), ('i', 'sg'), ('i', 'tau'), ('i', 'mu'), ('i', 'nu'), ('5',)], ["rho", "sg", "tau"]),
                    ([('i', 'mu'), ('i', 'nu'), ('i', 'mu'), ('i', 'nu')], [])]:
    Gx = prod([gamma(SR(T.lorentz_indices(a[1]))) if a[0] == 'i' else gamma5() for a in atoms])
    r = bmhv_trace(Gx)
    names = [a[1] for a in atoms if a[0] == 'i']
    summed = sorted({n for n in names if names.count(n) == 2})
    for _ in range(5):
        iv = {n: rnd.randint(0, 7) for n in free}
        tot = sum(prod(met8[v, v] for v in vals) * explicit8(atoms, dict(iv, **dict(zip(summed, vals))))
                  for vals in itertools.product(range(8), repeat=len(summed)))
        ok8 = ok8 and num8(r, iv) == tot
check("summed indices (one, two, and no g5) = explicit 16x16 gammas in d = 8, so the d-dependence is right too", ok8)

print("time %.0f s" % (time.time() - t_start))
print("ALL BMHV CHECKS PASS" if not FAIL else "FAILED: %s" % FAIL)
