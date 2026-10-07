# The Pythonic Dirac/Lorentz layer (feynsage.qft) against explicit 4x4 Dirac matrices and direct FORM.
# Run from the repository folder:  sage tests/test_qft_dirac.sage
import sys, time, itertools, random
sys.path.insert(0, '.')
from feynsage import *
from feynsage import form
from feynsage.qft.tensors import DOT, COMP, G, DELTA, EPS, TENSOR_NAMES, conjugate_vector
from sage.symbolic.function_factory import function as sfunction    # function() would rebind the names dot, g, comp here
t0 = time.time()
FAIL = []


def check(name, ok, extra=""):
    print("%-62s %s %s" % (name, "ok" if ok else "FAILED", extra))
    if not ok:
        FAIL.append(name)


def raises(f, exc=Exception):
    try:
        f()
    except exc:
        return True
    return False


# ---------------------------------------------------------------- explicit matrices
CC = QQbar
s1 = matrix(CC, [[0, 1], [1, 0]]); s2 = matrix(CC, [[0, -I], [I, 0]]); s3 = matrix(CC, [[1, 0], [0, -1]])
one2, zero2 = identity_matrix(CC, 2), zero_matrix(CC, 2)
gM = [block_matrix([[one2, zero2], [zero2, -one2]])] + [block_matrix([[zero2, s], [-s, zero2]]) for s in (s1, s2, s3)]
gM = [matrix(CC, x) for x in gM]
etaM = diagonal_matrix(QQ, [1, -1, -1, -1])
g5M = I * gM[0] * gM[1] * gM[2] * gM[3]                      # Peskin (3.68)
gE = [-I * gM[1], -I * gM[2], -I * gM[3], gM[0]]               # Euclidean: hermitian, square to 1
etaE = identity_matrix(QQ, 4)
g5E = gE[0] * gE[1] * gE[2] * gE[3]
Id4 = identity_matrix(CC, 4)
check("Euclidean gammas: {g_a, g_b} = 2 delta_ab",
      all(gE[a] * gE[b] + gE[b] * gE[a] == 2 * etaE[a, b] * Id4 for a in range(4) for b in range(4)))


def eps4(*ix):
    return sign(Permutation([i + 1 for i in ix])) if len(set(ix)) == 4 else 0


class World:
    """Numbers for a check: explicit gammas, metric (upper indices), vectors (upper components)."""

    def __init__(self, euclid, vecs):
        self.eu = euclid
        self.g = gE if euclid else gM
        self.g5 = g5E if euclid else g5M
        self.eta = etaE if euclid else etaM
        self.vecs = vecs                                        # name -> list of 4 upper components

    def low(self, p):
        return [self.eta[a, a] * self.vecs[str(p)][a] for a in range(4)]

    def slash(self, p):
        return sum(self.g[a] * self.low(p)[a] for a in range(4))

    def gam(self, a):
        return self.g[a]

    def ev(self, e, ix):
        """A Sage tensor expression at fixed values ix of its free indices."""
        e = SR(e)
        op = e.operator()
        if op is None:
            if e.is_numeric():
                return CC(e.pyobject()) if e.pyobject() in CC else CC(e)
            return CC(ix[str(e)]) if str(e) in ix else CC(e)
        args = e.operands()
        name = getattr(op, '__name__', str(op))
        name = TENSOR_NAMES.get(name, name)
        if name in ('add_vararg', 'add'):
            return sum(self.ev(a, ix) for a in args)
        if name in ('mul_vararg', 'mul'):
            return prod(self.ev(a, ix) for a in args)
        if name == 'pow':
            return self.ev(args[0], ix) ** int(args[1])
        if name == 'dot':
            a, b = args
            return sum(self.vecs[str(a)][k] * self.low(b)[k] for k in range(4))
        if name == 'comp':
            return CC(self.vecs[str(args[0])][ix[str(args[1])]])
        if name in ('g', 'delta'):
            return CC(self.eta[ix[str(args[0])], ix[str(args[1])]])
        if name == 'Eps':
            slots = []
            for a in args:
                if str(a) in ix:
                    slots.append([(ix[str(a)], 1)])
                else:                                           # a momentum: contract with p_a (lower)
                    slots.append([(k, self.low(a)[k]) for k in range(4)])
            tot = 0
            for combo in itertools.product(*slots):
                w = prod(c for _, c in combo)
                if w:
                    tot += eps4(*[k for k, _ in combo]) * w
            return CC(tot)
        if e.is_numeric():
            return CC(e)
        raise ValueError(name)


def rvec():
    return [QQ(random.randint(-9, 9)) / random.randint(1, 4) for _ in range(4)]


random.seed(int(7))
p, q, k, l, p1, p2, k1, k2 = momenta("p q k l p1 p2 k1 k2")
mu, nu, rho, sig, al = lorentz_indices("mu nu rho sig al")
m, M_, s, t = var('m M_ s t')
VEC = {str(x): rvec() for x in (p, q, k, l, p1, p2, k1, k2)}
WM, WE = World(False, VEC), World(True, VEC)
IDX = [str(x) for x in (mu, nu, rho, sig)]


def all_index_values(res, explicit, W, free):
    """Compare res (Sage) and explicit(ix) (a number) for every value of the free indices."""
    worst = 0
    for vals in itertools.product(range(4), repeat=len(free)):
        ix = dict(zip(free, vals))
        d = abs(W.ev(res, ix) - explicit(ix))
        worst = max(worst, d)
    return worst == 0


# ---------------------------------------------------------------- Minkowski traces
r = dirac_trace(gamma(mu) * gamma(nu))
check("Tr(g^mu g^nu) = 4 g^{mu nu}", bool((r - 4 * metric(mu, nu)).is_trivial_zero()), str(r))
check("  ... all index values vs matrices", all_index_values(r, lambda ix: (WM.gam(ix['mu']) * WM.gam(ix['nu'])).trace(), WM, ['mu', 'nu']))
r = dirac_trace(gamma(mu) * gamma(nu) * gamma(rho) * gamma(sig))
ref = 4 * (metric(mu, nu) * metric(rho, sig) - metric(mu, rho) * metric(nu, sig) + metric(mu, sig) * metric(nu, rho))
check("Tr(4 gammas) = 4(g g - g g + g g)  (Peskin 5.5)", bool((r - ref).expand().is_trivial_zero()))
r5 = dirac_trace(gamma(mu) * gamma(nu) * gamma(rho) * gamma(sig) * gamma5())
check("Tr(g^mu g^nu g^rho g^sig g5) = -4 i eps^{mu nu rho sig}  (Peskin 5.5)",
      bool((r5 + 4 * I * epsilon(mu, nu, rho, sig)).expand().is_trivial_zero()), str(r5))
check("  ... all 256 index values vs matrices with g5 = i g0 g1 g2 g3",
      all_index_values(r5, lambda ix: (prod(WM.gam(ix[n]) for n in IDX) * g5M).trace(), WM, IDX))
check("contract eps^{abcd} eps_{abcd} = -24  (Peskin 5.6)",
      bool(contract(epsilon(mu, nu, rho, sig) * epsilon(mu, nu, rho, sig)) == -24))
r = contract(epsilon(al, nu, rho, mu) * epsilon(al, nu, rho, sig))
check("eps^{a b c mu} eps_{a b c sig} = -6 g^mu_sig  (Peskin 5.6)", bool((r + 6 * metric(mu, sig)).is_trivial_zero()), str(r))
r = dirac_trace(slash(p) * slash(q))
check("Tr(p/ q/) = 4 p.q", bool((r - 4 * dot(p, q)).is_trivial_zero()))
r = dirac_trace(slash(p) * gamma(mu) * slash(q) * gamma(nu))
ref = form.compute(lines=[["p", "mu", "q", "nu"]], vectors=["p", "q"])
check("Tr(p/ g^mu q/ g^nu) vs direct FORM (form.compute)",
      all_index_values(r - ref.substitute_function(sfunction('g'), G).substitute_function(sfunction('dot'), DOT)
                       .substitute_function(sfunction('comp'), COMP), lambda ix: 0, WM, ['mu', 'nu']))
check("  ... vs matrices",
      all_index_values(r, lambda ix: (WM.slash(p) * WM.gam(ix['mu']) * WM.slash(q) * WM.gam(ix['nu'])).trace(), WM, ['mu', 'nu']))
x = [p, q, k, l, p1, p2]
r = dirac_trace(gamma5() * prod(slash(y) for y in x))
check("Tr(g5 a/ b/ c/ d/ e/ f/) vs matrices",
      WM.ev(r, {}) == (g5M * prod(WM.slash(y) for y in x)).trace())
expr = (slash(p) + m) * gamma(mu) * PL() * (slash(q) - M_) * gamma(nu) * PR() * slash(k)
r = dirac_trace(expr, rules={m: 3, M_: 2})
explicit = lambda ix: ((WM.slash(p) + 3) * WM.gam(ix['mu']) * (Id4 - g5M) / 2 * (WM.slash(q) - 2) * WM.gam(ix['nu'])
                       * (Id4 + g5M) / 2 * WM.slash(k)).trace()
check("Tr((p/+m) g^mu PL (q/-M) g^nu PR k/) vs matrices (rules m=3, M=2)", all_index_values(r, explicit, WM, ['mu', 'nu']))
r = dirac_trace(gamma(mu) * slash(p) * gamma(nu) * slash(q) * gamma(mu) * slash(k) * gamma(nu) * slash(l))
check("Tr with repeated indices summed",
      WM.ev(r, {}) == sum(etaM[a, a] * etaM[b, b] * (WM.gam(a) * WM.slash(p) * WM.gam(b) * WM.slash(q) * WM.gam(a)
                                                     * WM.slash(k) * WM.gam(b) * WM.slash(l)).trace()
                          for a in range(4) for b in range(4)))
r = dirac_trace(sigma(mu, nu) * sigma(rho, sig))
check("Tr(sigma^{mu nu} sigma^{rho sig}) vs matrices",
      all_index_values(r, lambda ix: ((I / 2) ** 2 * (WM.gam(ix['mu']) * WM.gam(ix['nu']) - WM.gam(ix['nu']) * WM.gam(ix['mu']))
                                      * (WM.gam(ix['rho']) * WM.gam(ix['sig']) - WM.gam(ix['sig']) * WM.gam(ix['rho']))).trace(), WM, IDX))
r = dirac_trace(gamma5() * gamma(mu) * gamma(nu) * slash(p) * slash(q))
check("Tr(g5 g^mu g^nu p/ q/) = -4 i eps(mu,nu,p,q)",
      bool((r + 4 * I * epsilon(mu, nu, p, q)).expand().is_trivial_zero()), str(r))
rs = dirac_trace([gamma(mu) * gamma(mu), slash(p) * slash(p), 7])
print("   batch:", rs)
check("batch: [g^mu g_mu, p/ p/, 7] -> [16, 4 p^2, 28]",
      bool(rs[0] == 16) and bool((rs[1] - 4 * dot(p, p)).is_trivial_zero()) and bool(rs[2] == 28))

# ---------------------------------------------------------------- kinematics
r = dirac_trace(slash(p) * slash(k), rules={dot(p, k): s / 2})
check("rules: dot(p, k) -> s/2", bool((r - 2 * s).is_trivial_zero()))
r = dirac_trace(slash(p) * slash(k) * slash(p) * slash(k), rules={dot(p, p): m ^ 2, dot(k, k): 0, dot(p, k): s / 2})
check("rules on a longer trace", bool((r - (8 * (s / 2) ^ 2 - 4 * m ^ 2 * 0)).expand().is_trivial_zero()), str(r))
r1 = dirac_trace(slash(k2) * slash(p1), rules={k2: p1 + p2 - k1})
r2 = dirac_trace(slash(p1 + p2 - k1) * slash(p1))
check("momentum rule k2 -> p1 + p2 - k1 inside a slash", bool((r1 - r2).expand().is_trivial_zero()))
r = dirac_trace(slash(k2) * slash(p1) / dot(p1, k2), rules={k2: p1 + p2 - k1})
check("momentum rule also inside a denominator 1/(p1.k2)", bool((r - 4).simplify_rational() == 0), str(r))
r = dirac_trace(slash(p) * slash(q)).subs({dot(p, q): s})
check("result.subs({dot(p, q): s})", bool(r == 4 * s), str(r))

# ---------------------------------------------------------------- d dimensions
d = var('d')
r = dirac_trace(gamma(mu) * gamma(nu) * gamma(mu) * gamma(nu), dim=d)
check("Tr(g^mu g^nu g_mu g_nu) = 4 d (2 - d) in d dims", bool((r - 4 * d * (2 - d)).expand().is_trivial_zero()), str(r))
r = dirac_trace(gamma5() * gamma(mu) * gamma5() * gamma(nu), dim=d)
check("NDR: Tr(g5 g^mu g5 g^nu) = -4 g^{mu nu}", bool((r + 4 * metric(mu, nu)).is_trivial_zero()), str(r))
r = dirac_trace(PL() * gamma(mu) * PR() * gamma(nu), dim=d)
check("NDR: Tr(PL g^mu PR g^nu) = 2 g^{mu nu}", bool((r - 2 * metric(mu, nu)).is_trivial_zero()), str(r))
check("NDR refuses Tr(g^mu g^nu g^rho g^sig PR) in d dims",
      raises(lambda: dirac_trace(gamma(mu) * gamma(nu) * gamma(rho) * gamma(sig) * PR(), dim=d), DiracError))
check("NDR refuses Tr(g^mu g^nu g^rho g^sig g5) in d dims",
      raises(lambda: dirac_trace(gamma(mu) * gamma(nu) * gamma(rho) * gamma(sig) * gamma5(), dim=d), DiracError))
check("eps contracted in d dims needs eps_in_d=True",
      raises(lambda: contract(epsilon(mu, nu, rho, sig) * epsilon(mu, nu, rho, sig), dim=d), DiracError))

# ---------------------------------------------------------------- Euclidean
set_convention(euclidean=True)
r = dirac_trace(gamma(mu) * gamma(nu) * gamma(rho) * gamma(sig) * gamma5())
check("Euclidean: Tr(g_mu g_nu g_rho g_sig g5) = 4 eps_{mu nu rho sig}",
      bool((r - 4 * epsilon(mu, nu, rho, sig)).expand().is_trivial_zero()), str(r))
check("  ... all 256 index values vs matrices with g5 = g1 g2 g3 g4",
      all_index_values(r, lambda ix: (prod(WE.gam(ix[n]) for n in IDX) * g5E).trace(), WE, IDX))
r = dirac_trace(gamma5() * prod(slash(y) for y in x))
check("Euclidean: Tr(g5 a/ ... f/) vs matrices", WE.ev(r, {}) == (g5E * prod(WE.slash(y) for y in x)).trace())
r = contract(epsilon(mu, nu, rho, sig) * epsilon(mu, nu, rho, sig))
check("Euclidean: eps_{abcd} eps_{abcd} = +24", bool(r == 24), str(r))
r = dirac_trace(slash(p) * gamma(mu) * slash(q) * gamma(nu))
check("Euclidean: Tr(p/ g_mu q/ g_nu) vs matrices",
      all_index_values(r, lambda ix: (WE.slash(p) * WE.gam(ix['mu']) * WE.slash(q) * WE.gam(ix['nu'])).trace(), WE, ['mu', 'nu']))
check("Euclidean refuses a Minkowski metric", raises(lambda: contract(G(mu, nu) * comp(p, mu))))
set_convention(euclidean=False)

# ---------------------------------------------------------------- spinors, conjugation, spin sums
J = ubar(p2, m) * gamma(mu) * u(p1, m)
check("ubar(p2) g^mu u(p1) is a number (scalar kind)", J.kind() == 'scalar')
check("u(p) * ubar(p) is refused (use spin_sum)", raises(lambda: u(p, m) * ubar(p, m), DiracError))
check("gamma(mu) * ubar(p) is refused", raises(lambda: gamma(mu) * ubar(p, m), DiracError))

# conjugation: [abar G b]^* = bbar Gbar a for any spinors a, b (abar = a^dagger g0): explicit numbers
A = vector(CC, [1 + 2 * I, -1, 3 * I, 2]); B = vector(CC, [2, I, -1 - I, 5])
bar = lambda x: x.conjugate() * gM[0]


def num_bilinear(expr, ix, sp):
    tot = 0
    for c, o, cl in expr.terms:
        val = WM.ev(c, ix)
        for ch in cl:
            mat = Id4
            for F in ch.factors:
                fm = 0
                for a, cf in F.items:
                    am = {'1': Id4, '5': g5M}.get(a[0])
                    if am is None:
                        am = WM.gam(ix[str(a[1])]) if a[0] == 'i' else WM.slash(a[1])
                    fm = fm + WM.ev(cf, ix) * am
                mat = mat * fm
            left = sp[ch.left.key()]
            right = sp[ch.right.key()]
            val *= left * mat * right
        tot += val
    return tot


sp = {('ubar', 'p2', 'm'): bar(A), ('u', 'p1', 'm'): B, ('ubar', 'p1', 'm'): bar(B), ('u', 'p2', 'm'): A}
for nm, G_ in [("1", one()), ("g5", gamma5()), ("g^mu", gamma(mu)), ("g^mu g5", gamma(mu) * gamma5()),
               ("sigma^{mu nu}", sigma(mu, nu)), ("PL g^mu PR p/", PL() * gamma(mu) * PR() * slash(p)),
               ("(1+2i) g^mu q/ g^nu", (1 + 2 * I) * gamma(mu) * slash(q) * gamma(nu))]:
    Mx = ubar(p2, m) * G_ * u(p1, m)
    ok = all(num_bilinear(conjugate(Mx), {'mu': a, 'nu': b}, sp) == num_bilinear(Mx, {'mu': a, 'nu': b}, sp).conjugate()
             for a in range(4) for b in range(4))
    check("conjugate([ubar %s u]) vs explicit spinors" % nm, ok)
Mx = ubar(p2, m) * gamma(mu) * u(p1, m) * ubar(k2) * gamma(mu) * u(k1)
check("conjugate renames summed indices (mu -> mu_c)", [str(i) for i in dummy_indices(conjugate(Mx))] == ['mu_c'])

# e+ e- -> mu+ mu- : M = vbar(p2) g^mu u(p1) * ubar(k1) g_mu v(k2), electrons massless
Mee = vbar(p2) * gamma(mu) * u(p1) * ubar(k1, M_) * gamma(mu) * v(k2, M_)
r = spin_sum(Mee * conjugate(Mee))
ref = form.compute(lines=[["p2", "mu", "p1", "nu"], ["k1 + Mmu", "mu", "k2 - Mmu", "nu"]], vectors=["p1", "p2", "k1", "k2"])
ref = contract(ref.subs({var('Mmu'): M_}))      # contract() reads form.compute output and puts it in canonical order
check("spin_sum |M|^2 for e+e- -> mu+mu- vs direct FORM", bool((r - ref).expand().is_trivial_zero()), str(r))
ref2 = 32 * (dot(p1, k1) * dot(p2, k2) + dot(p1, k2) * dot(p2, k1) + M_ ^ 2 * dot(p1, p2))
check("  ... = 32[(p1.k1)(p2.k2) + (p1.k2)(p2.k1) + M^2 p1.p2]", bool((r - ref2).expand().is_trivial_zero()))
r4 = spin_sum(Mee * conjugate(Mee), dim=d)
check("  ... in d dims, at d = 4 the same", bool((r4.subs(d=4) - r).expand().is_trivial_zero()))
check("spin_sum needs partners", raises(lambda: spin_sum(Mee * conjugate(ubar(p2) * u(p1))), DiracError))


# ---------------------------------------------------------------- open lines (contract / simplify_dirac)
check("contract(g^{mu nu} g_mu g_nu) = 4", bool(contract(metric(mu, nu) * gamma(mu) * gamma(nu)) == 4))
check("  ... = d in d dims", bool(contract(metric(mu, nu) * gamma(mu) * gamma(nu), dim=d) == d))
r1 = contract(gamma(mu) * gamma(nu) * gamma(mu))
r2 = contract(gamma(mu) * gamma(nu) * gamma(rho) * gamma(mu))
r3 = contract(gamma(mu) * gamma(nu) * gamma(rho) * gamma(sig) * gamma(mu))
check("Peskin (5.9): g^mu g^nu g_mu = -2 g^nu", repr(r1) == "-2*gamma(nu)", repr(r1))
check("Peskin (5.9): g^mu g^nu g^rho g_mu = 4 g^{nu rho}", bool((r2 - 4 * metric(nu, rho)).is_trivial_zero()), repr(r2))
X = r3 + 2 * gamma(sig) * gamma(rho) * gamma(nu)
check("Peskin (5.9): g^mu g^nu g^rho g^sig g_mu = -2 g^sig g^rho g^nu (reversed order)",
      bool(SR(contract(X)) == 0), repr(contract(X)))
r = contract(gamma(mu) * gamma(nu) * gamma(mu), dim=d)
check("g^mu g^nu g_mu = (2 - d) g^nu", repr(r) == "(-d + 2)*gamma(nu)", repr(r))
r = contract(gamma(mu) * slash(p) * slash(q) * gamma(mu))
check("g^mu p/ q/ g_mu = 4 p.q (4 dims)", bool((r - 4 * dot(p, q)).is_trivial_zero()), repr(r))
r = contract(ubar(p2, m) * slash(p2) * gamma(mu) * slash(p1) * u(p1, m))
check("Dirac equation at both ends: ubar(p2) p2/ g^mu p1/ u(p1) = m^2 ubar g^mu u",
      repr(r) == "m^2*ubar(p2, m)*gamma(mu)*u(p1, m)", repr(r))
r = contract(vbar(p2, m) * gamma(mu) * gamma5() * slash(p1) * v(p1, m))
check("v spinor through gamma5: vbar g^mu g5 p1/ v(p1) = -m vbar g^mu g5 v",
      repr(r) == "-m*vbar(p2, m)*gamma(mu)*gamma5()*v(p1, m)", repr(r))
# open-line result against the trace: Tr[contract(G) X] = Tr[G X]
Gx = gamma(mu) * slash(p) * gamma(nu) * slash(q) * gamma(mu) * gamma5()
r1 = dirac_trace(contract(Gx) * slash(k) * gamma(nu))
r2 = dirac_trace(Gx * slash(k) * gamma(nu))
check("Tr[contract(G) X] = Tr[G X] (with gamma5)", bool((r1 - r2).expand().is_trivial_zero()))

# ---------------------------------------------------------------- Compton scattering, Peskin (5.87)
kp, pp = momenta("kp pp")
e_in = polarization("e_in", k)                     # incoming photon eps(k)
e_out = polarization("e_out", kp)                  # outgoing photon: the amplitude has eps*(k') = e_out_c
eo = conjugate_vector(e_out)
Mc = (ubar(pp, m) * slash(eo) * (slash(p) + slash(k) + m) * slash(e_in) * u(p, m) / (2 * dot(p, k))
      + ubar(pp, m) * slash(e_in) * (slash(p) - slash(kp) + m) * slash(eo) * u(p, m) / (-2 * dot(p, kp)))
kin = {dot(p, p): m ^ 2, dot(pp, pp): m ^ 2, dot(k, k): 0, dot(kp, kp): 0}
r = spin_sum(Mc * conjugate(Mc), rules={pp: p + k - kp})
r = polarization_sum(r, [e_in, e_out])
r = (r / 4).subs({dot(p, p): m ^ 2, dot(k, k): 0, dot(kp, kp): 0, dot(k, kp): dot(p, k) - dot(p, kp)})
pk, pkp = dot(p, k), dot(p, kp)
ref = 2 * (pkp / pk + pk / pkp + 2 * m ^ 2 * (1 / pk - 1 / pkp) + m ^ 4 * (1 / pk - 1 / pkp) ^ 2)
check("Compton (1/4) sum |M|^2 = Peskin (5.87), Feynman-gauge photon sums",
      bool((r - ref).simplify_rational().is_trivial_zero()))
r2 = polarization_sum(spin_sum(Mc * conjugate(Mc), rules={pp: p + k - kp}), [e_in, e_out], gauge="physical", n=p)
r2 = (r2 / 4).subs({dot(p, p): m ^ 2, dot(k, k): 0, dot(kp, kp): 0, dot(k, kp): dot(p, k) - dot(p, kp)})
check("  ... the same with physical polarization sums (Ward identity)", bool((r2 - ref).simplify_rational().is_trivial_zero()))
check("polarization_sum refuses a term quadratic in eps",
      raises(lambda: polarization_sum(dot(e_in, p) ^ 2, e_in), DiracError))


# ---------------------------------------------------------------- conventions and notation
set_convention("note")
Dn = spacetime_dimension()
check('convention "note": the dimension is D', str(Dn) == 'D')
r = dirac_trace(gamma(mu) * gamma(nu) * gamma(mu) * gamma(nu), dim="D")
check('dim="D" works like dim=D', bool((r - 4 * Dn * (2 - Dn)).expand().is_trivial_zero()), str(r))
set_convention("standard")
check('convention "standard": the dimension is d', str(spacetime_dimension()) == 'd')
check("to_euclidean: s = p^2 -> -p_E^2", bool(to_euclidean(dot(p, p) - m ^ 2) == -dot(p, p) - m ^ 2))
check("to_minkowski undoes to_euclidean", bool((to_minkowski(to_euclidean(dot(p, q) * dot(k, l))) - dot(p, q) * dot(k, l)).is_trivial_zero()))
check("to_euclidean refuses free indices", raises(lambda: to_euclidean(comp(p, mu)), TypeError))
check("LaTeX: gamma^mu p/ (p/+m)", latex(gamma(mu) * slash(p) * (slash(q) + m)) == r"\gamma^{\mu} {\not{p}} \left({\not{q}} + m\right)",
      latex(gamma(mu) * slash(p) * (slash(q) + m)))
check("LaTeX: [ubar(p2) g^mu u(p1)]", latex(ubar(p2, m) * gamma(mu) * u(p1, m)) == r"\left[\bar{u}(p_{2}) \gamma^{\mu} u(p_{1})\right]")
check("LaTeX: p.q, p^2, p^mu, g^{mu nu}, eps",
      latex(dot(p, q) + dot(p, p) + comp(p, mu) * comp(q, nu) + metric(mu, nu) + epsilon(mu, nu, p, q)).count(r"\cdot") == 1)
check("LaTeX: (p.q)^2 and (p^2)^2 are valid",
      r"\left({p} \cdot {q}\right)^{2}" in latex(dot(p, q) ^ 2) and latex(dot(p, p) ^ 2) == "{{p}^{2}}^{2}",
      latex(dot(p, q) ^ 2) + "   " + latex(dot(p, p) ^ 2))
set_convention(euclidean=True)
check("LaTeX in the Euclidean convention: p_E^2", latex(dot(p, p)) == "{{p}_{E}^{2}}", latex(dot(p, p)))
set_convention(euclidean=False)


# ---------------------------------------------------------------- SU(N) colour
ca, cb, cc, cd = color_indices("ca cb cc cd")
ci, cj, ck, cl = quark_colors("ci cj ck cl")
Nc = var('N')
from feynsage.qft.color import KRON
r = color_factor(T_color(ca, ci, cj) * T_color(cb, cj, ci))
check("Tr(T^a T^b) = delta^{ab}/2", bool((r - KRON(ca, cb) / 2).is_trivial_zero()), str(r))
r = color_factor(color_chain(ca, ca, i=ci, j=cj))
check("(T^a T^a)_{ij} = C_F delta_ij", bool((r - (Nc ^ 2 - 1) / (2 * Nc) * KRON(ci, cj)).expand().is_trivial_zero()), str(r))
r = color_factor(f_color(ca, cc, cd) * f_color(cb, cc, cd))
check("f^{acd} f^{bcd} = N delta^{ab}", bool((r - Nc * KRON(ca, cb)).expand().is_trivial_zero()), str(r))
r = color_factor(color_trace(ca, cb, ca, cb))
check("Tr(T^a T^b T^a T^b) = -(N^2-1)/(4N)", bool((r + (Nc ^ 2 - 1) / (4 * Nc)).expand().is_trivial_zero()), str(r))
r = color_factor(f_color(ca, cb, cc) ^ 2, N=3)
check("f^{abc} f^{abc} = N (N^2 - 1) = 24 for SU(3)", bool(r == 24), str(r))
# explicit Gell-Mann matrices: T^a = lambda^a / 2
lam = [matrix(CC, 3, 3, x) for x in (
    [0, 1, 0, 1, 0, 0, 0, 0, 0], [0, -I, 0, I, 0, 0, 0, 0, 0], [1, 0, 0, 0, -1, 0, 0, 0, 0],
    [0, 0, 1, 0, 0, 0, 1, 0, 0], [0, 0, -I, 0, 0, 0, I, 0, 0], [0, 0, 0, 0, 0, 1, 0, 1, 0],
    [0, 0, 0, 0, 0, -I, 0, I, 0])] + [matrix(CC, 3, 3, [1, 0, 0, 0, 1, 0, 0, 0, -2]) / sqrt(QQbar(3))]
Tg = [x / 2 for x in lam]
TT = [[Tg[x] * Tg[y] for y in range(8)] for x in range(8)]
CFm = sum(Tg[x] * Tg[x] for x in range(8))
expl = sum((TT[a][b] * Tg[c]).trace() * (Tg[c] * TT[b][a] * CFm).trace() for a in range(8) for b in range(8) for c in range(8))
r = color_factor(color_trace(ca, cb, cc) * color_trace(cc, cb, ca, cd, cd), N=3)
check("Tr(T^a T^b T^c) Tr(T^c T^b T^a T^d T^d) vs Gell-Mann matrices", bool(QQbar(r) == expl), "%s vs %s" % (r, expl))
check("Dirac: an index three times is refused (trace and contract)",
      raises(lambda: dirac_trace(gamma(mu) * gamma(mu) * slash(p) * gamma(mu)), DiracError)
      and raises(lambda: contract(gamma(mu) * gamma(mu) * gamma(mu)), DiracError))
check("contract(..., eps_in_d=True) reaches the open-line path",
      not raises(lambda: contract(epsilon(mu, nu, rho, sig) * gamma(mu) * gamma(nu), dim=d, eps_in_d=True)))
check("mixed list keeps eps_in_d", not raises(lambda: contract([epsilon(mu, nu, rho, sig) * epsilon(mu, nu, rho, sig), gamma(mu) * gamma(mu)],
                                                              dim=d, eps_in_d=True)))
check("an index three times is refused", raises(lambda: color_factor(color_trace(ca, cb) * color_trace(ca, cb, cb)), Exception))
fabc = lambda a, b, c: -2 * I * ((Tg[a] * Tg[b] * Tg[c]).trace() - (Tg[a] * Tg[c] * Tg[b]).trace())
expl = sum(fabc(a, b, c) * (Tg[a] * Tg[b] * Tg[c]).trace() for a in range(8) for b in range(8) for c in range(8))
r = color_factor(f_color(ca, cb, cc) * color_trace(ca, cb, cc), N=3)
check("f^{abc} Tr(T^a T^b T^c) = i N (N^2-1)/4 = 6i, vs Gell-Mann matrices", bool(QQbar(r) == expl), "%s vs %s" % (r, expl))
# colour inside an amplitude: q qbar -> gluon* -> mu mu like current with T^a; sum over colours
Mq = vbar(p2) * gamma(mu) * u(p1) * T_color(ca, cj, ci) * ubar(k1) * gamma(mu) * v(k2)
r = spin_sum(Mq * conjugate(Mq))
M0 = vbar(p2) * gamma(mu) * u(p1) * ubar(k1) * gamma(mu) * v(k2)
r0 = spin_sum(M0 * conjugate(M0))
check("spin_sum with colour: sum |T^a_ji|^2 = (N^2-1)/2 times the Dirac part",
      bool((r - (Nc ^ 2 - 1) / 2 * r0).expand().is_trivial_zero()), str(r))

check("latex_name: momenta('pq', latex_name=\"p'\") prints p'", latex(slash(momenta("pq", latex_name="p'"))) == r"{\not{{p'}}}",
      latex(slash(momenta("pq", latex_name="p'"))))
ex1 = polarization("ex1", k, latex_name=r"\epsilon")
check("polarization latex: conjugate prints with a star", latex(conjugate_vector(ex1)) == r"{{\epsilon}^{*}}", latex(conjugate_vector(ex1)))
check("color_delta: delta_ij delta_ji = N", bool(color_factor(color_delta(ci, cj) * color_delta(cj, ci)) == Nc))

# ---------------------------------------------------------------- massive vector boson (P&S problem 5.5)
# e- (q1) e+ (q2) -> gamma (q3) B (q4), massless electrons; Xianyu's solution eq. (5.87); full massive sum for B
q1, q2, q3, q4 = momenta("q1 q2 q3 q4")
gB, MBm, sB, tB, uB = SR.var('gB MBm sB tB uB')
eph = polarization("eph", q3); eBm = polarization("eBm", q4, mass=MBm)
phc, Bc = conjugate_vector(eph), conjugate_vector(eBm)
MBa = (vbar(q2) * slash(Bc) * (slash(q1) - slash(q3)) * slash(phc) * u(q1) / (-2 * dot(q1, q3))
       + vbar(q2) * slash(phc) * (slash(q1) - slash(q4)) * slash(Bc) * u(q1) / (MBm ^ 2 - 2 * dot(q1, q4)))
kinB = {dot(q1, q1): 0, dot(q2, q2): 0, dot(q3, q3): 0, dot(q4, q4): MBm ^ 2, dot(q1, q2): sB / 2,
        dot(q1, q3): -tB / 2, dot(q2, q4): (MBm ^ 2 - tB) / 2, dot(q1, q4): (MBm ^ 2 - uB) / 2,
        dot(q2, q3): -uB / 2, dot(q3, q4): (sB - MBm ^ 2) / 2}
rB = (polarization_sum(spin_sum(MBa * conjugate(MBa)), [eph, eBm]) / 4).subs(kinB).subs(sB=MBm ^ 2 - tB - uB)
x587 = 2 * (uB / tB + tB / uB + 2 * (MBm ^ 2 - tB - uB) * MBm ^ 2 / (tB * uB))
check("e+e- -> gamma B with -g + kk/M^2: 2[u/t + t/u + 2 s M^2/(t u)]", bool((rB - x587).simplify_rational() == 0), str(rB))
MBk = (vbar(q2) * slash(q4) * (slash(q1) - slash(q3)) * slash(phc) * u(q1) / (-2 * dot(q1, q3))
       + vbar(q2) * slash(phc) * (slash(q1) - slash(q4)) * slash(q4) * u(q1) / (MBm ^ 2 - 2 * dot(q1, q4)))
rk = polarization_sum(spin_sum(MBk * conjugate(MBk)), eph).subs(kinB).subs(sB=MBm ^ 2 - tB - uB)
check("Ward identity for B: |k_B . M|^2 = 0", bool(rk.simplify_rational() == 0), str(rk))

# ---------------------------------------------------------------- errors a student may hit
check("slash(p + m) is refused (write slash(p) + m)", raises(lambda: slash(p + m), TypeError))
check("gamma(p) for a momentum is refused (write slash(p))", raises(lambda: gamma(p), TypeError))
check("gamma(x) of an undeclared symbol is caught in a trace",
      raises(lambda: dirac_trace(gamma(var('muu')) * gamma(nu)), Exception))
check("gamma(5) is still Euler's Gamma", bool(gamma(5) == 24))
check("names clashing with FORM internals (symbol FSE1, E1)",
      bool(dirac_trace(var('FSE1') * var('E1') * slash(p) * slash(p)) == 4 * var('FSE1') * var('E1') * dot(p, p)))

print("time %.1f s" % (time.time() - t0))
print("ALL QFT-LAYER CHECKS PASS" if not FAIL else "FAILED: %s" % FAIL)
