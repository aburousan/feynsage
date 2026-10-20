r"""
Tensor decomposition of loop integrals with any number of loops (FeynCalc's Tdec).

    tdec([("l1", "mu"), ("l2", "nu"), ...], ["p1", "p2"], rules=None)

Int ... l1^mu l2^nu ... f(l, p) = sum_T c_T T^{mu nu ...}: the T are all tensors made of metrics and the external
momenta, the c_T scalar integrals (numerators l_i.p_j, l_i.l_j).  The c_T solve G c = R with G_ab = T_a . T_b
and R_a = T_a . (l1^mu l2^nu ...), all indices contracted in d dimensions.  Two tensors are contracted as a
graph: their index positions, joined by the metrics of both, form closed cycles (a factor d each) and open
paths between two momenta (a scalar product).  Zero Gram determinants are not supported (as in FeynCalc).
"""
from sage.all import SR, matrix, vector

from .qft import tensors as T


def _matchings(items):
    """All ways to pick a set of disjoint pairs from items (the metric part of a tensor)."""
    if not items:
        yield []
        return
    first, rest = items[0], items[1:]
    for m in _matchings(rest):                       # first is not in a pair
        yield m
    for k, other in enumerate(rest):                 # first is paired with other
        for m in _matchings(rest[:k] + rest[k + 1:]):
            yield [(first, other)] + m


def _basis(n, externals):
    """Tensors of rank n: (pairs of positions joined by a metric, {position: external momentum})."""
    from itertools import product
    out = []
    for pairs in _matchings(list(range(n))):
        free = [i for i in range(n) if all(i not in p for p in pairs)]
        for moms in product(externals, repeat=len(free)):
            out.append((pairs, dict(zip(free, moms))))
    return out


def _attach(Tn, n):
    """slot -> ("e", partner) or ("m", momentum) for one tensor."""
    at = {}
    for a, b in Tn[0]:
        at[a], at[b] = ("e", b), ("e", a)
    for i, p in Tn[1].items():
        at[i] = ("m", p)
    return at


def _contract(A, B, n, sp, d):
    """A . B with all n index positions contracted.  Every slot has one attachment from A and one from B (a
    metric partner or a momentum), so the slots fall into alternating paths from a momentum to a momentum
    (their scalar product) and closed alternating cycles of metrics (a factor d each)."""
    at = (_attach(A, n), _attach(B, n))
    done, val = set(), 1
    for i in range(n):
        for side in (0, 1):
            if (i, side) in done or at[side][i][0] != "m":
                continue
            start = at[side][i][1]                     # a path starts at this momentum
            slot, s = i, 1 - side                      # look at the other tensor's attachment of the slot
            done.add((i, side))
            while True:
                kind, x = at[s][slot]
                done.add((slot, s))
                if kind == "m":
                    val *= sp(start, x)
                    break
                slot = x                               # follow the metric to its partner slot,
                done.add((slot, s))                    # then the partner's attachment in the other tensor
                s = 1 - s
    for i in range(n):                                 # what is left are closed cycles
        if (i, 0) in done:
            continue
        val *= d
        slot, s = i, 0
        while (slot, s) not in done:
            done.add((slot, s))
            slot = at[s][slot][1]
            done.add((slot, s))
            s = 1 - s
    return val


def _key(A):
    return (tuple(sorted(tuple(sorted(p)) for p in A[0])), tuple(sorted((i, str(m)) for i, m in A[1].items())))


def _permute(A, pm):
    return ([(pm[a], pm[b]) for a, b in A[0]], {pm[i]: m for i, m in A[1].items()})


def _group_perms(groups, n):
    """All permutations of the positions that only exchange positions within a group."""
    from itertools import permutations, product
    out = []
    for choice in product(*[list(permutations(g)) for g in groups]):
        pm = list(range(n))
        for g, img in zip(groups, choice):
            for a, b in zip(g, img):
                pm[a] = b
        out.append(pm)
    return out


def tdec(indices, externals, rules=None, dim=None, check=True):
    """
    The tensor decomposition  l_{i1}^{mu1} ... l_{in}^{mun} -> sum_T c_T T^{mu1 ... mun}  under a loop integral.
    indices: [(loop momentum, Lorentz index), ...] as strings or feynsage symbols; externals: the external
    momenta the result may contain; rules: scalar products of the externals, e.g. {dot(p1, p1): 0} or
    {"p1.p2": "s/2"} (FeynCalc's FinalSubstitutions); dim: the dimension (default the symbol d of the qft layer).
    Returns a Sage expression in metric(mu, nu), comp(p, mu) and dot(l, p), dot(l1, l2).  Positions with the same
    loop momentum are symmetrised first (Pak), so l^mu l^nu l^rho l^sigma needs few tensors.  check=True
    contracts the result with basis tensors again, independently (feynsage's contract, FORM), and compares.
    """
    from sage.all import PolynomialRing, QQ
    from .qft.tensors import dot, comp, metric
    from ._parse import sr
    if dim is None:
        dim = T.dimension()
    vec = lambda x: SR(T.momenta(x)) if isinstance(x, str) else SR(x)
    loops = [vec(l) for l, _ in indices]
    idx = [SR(T.lorentz_indices(m)) if isinstance(m, str) else SR(m) for _, m in indices]
    ps = [vec(p) for p in externals]
    n = len(indices)
    rl = {}
    for k, v in (rules or {}).items():
        if isinstance(k, str):
            a_, b_ = [vec(x.strip()) for x in k.replace("^2", "." + k.split("^")[0]).split(".")]
            k = dot(a_, b_)
        rl[SR(k)] = sr(v) if isinstance(v, str) else SR(v)
    # the ring: d, every scalar product without a rule, the symbols of the rules
    names = sorted({str(x) for x in loops} | {str(x) for x in ps})
    pairs = [(a_, b_) for i, a_ in enumerate(names) for b_ in names[i:]]
    free = [pr for pr in pairs if SR(dot(SR(pr[0]), SR(pr[1]))) not in rl]
    extra = sorted({str(v) for val in rl.values() for v in SR(val).variables()} - {str(dim)})
    gens = [str(dim)] + ["dot_%s_%s" % pr for pr in free] + extra
    P = PolynomialRing(QQ, gens)
    Fq = P.fraction_field()
    g = dict(zip(gens, P.gens()))

    def sp(a_, b_):
        a_, b_ = sorted([str(a_), str(b_)])
        key = SR(dot(SR(a_), SR(b_)))
        if key in rl:
            return P(str(rl[key]))
        return g["dot_%s_%s" % (a_, b_)]
    d_ = g[str(dim)]
    basis = _basis(n, ps)
    # Pak's symmetrisation: positions of the same loop momentum can be exchanged, so the coefficients of the
    # tensors of one orbit are equal; decompose onto the orbit sums (much smaller systems)
    groups = {}
    for i in range(n):
        groups.setdefault(str(loops[i]), []).append(i)
    perms = _group_perms(list(groups.values()), n)
    orbits, where = [], {}
    for A in basis:
        if _key(A) in where:
            continue
        orb = {}
        for pm in perms:
            B = _permute(A, pm)
            orb[_key(B)] = B
        for kk in orb:
            where[kk] = len(orbits)
        orbits.append(list(orb.values()))
    qten = ([], {i: loops[i] for i in range(n)})
    G = matrix(Fq, [[sum(_contract(A, B, n, sp, d_) for A in OA for B in OB) for OB in orbits] for OA in orbits])
    R = vector(Fq, [sum(_contract(A, qten, n, sp, d_) for A in OA) for OA in orbits])
    if G.change_ring(P).det() == 0:                    # (det over the fraction field is very slow)
        raise ZeroDivisionError("tdec: the Gram determinant of the tensor basis vanishes for these kinematics")
    c = G.solve_right(R)
    back = {str(dim): SR(dim)}
    back.update({"dot_%s_%s" % pr: SR(dot(SR(pr[0]), SR(pr[1]))) for pr in free})
    back.update({v: SR(v) for v in extra})
    to_sr = lambda x: SR(str(x.numerator())).subs({SR.var(k): v for k, v in back.items()}) / \
        SR(str(x.denominator())).subs({SR.var(k): v for k, v in back.items()})

    def to_frac(e):
        """A scalar Sage expression in Dot(a, b), d and the rule symbols as an element of the fraction field."""
        e = SR(e)
        w0, w1 = SR.wild(0), SR.wild(1)
        op = SR(dot(loops[0], loops[0])).operator()       # the symbolic function Dot itself (dot() checks its input)
        sub = {}
        for D_ in e.find(op(w0, w1)):
            a_, b_ = sorted(str(x) for x in D_.operands())
            key = SR(dot(SR(a_), SR(b_)))
            sub[D_] = SR(rl[key]) if key in rl else SR.var("dot_%s_%s" % (a_, b_))
        e = e.subs(sub)
        return Fq(P(str(e.numerator()))) / Fq(P(str(e.denominator())))

    def tensor(A):
        t = SR(1)
        for a_, b_ in A[0]:
            t *= metric(idx[a_], idx[b_])
        for i, p in A[1].items():
            t *= comp(p, idx[i])
        return t
    out = sum(to_sr(ci).factor() * sum(tensor(A) for A in O) for ci, O in zip(c, orbits) if ci != 0)
    if check:
        from .qft.operations import contract
        for t in [tensor(A) for A in basis[:6]]:      # FORM: (result . T_a) = (l ... l . T_a)
            rhs = SR(1)
            for i in range(n):
                rhs *= comp(loops[i], idx[i])
            diff = SR(contract(out * t, rules=rl, dim=dim)) - SR(contract(rhs * t, rules=rl, dim=dim))
            if to_frac(diff) != 0:                     # exact, in the rational functions (not Maxima)
                raise RuntimeError("tdec: the decomposition does not contract back; please report")
    return out
