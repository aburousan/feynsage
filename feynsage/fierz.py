r"""
Products of two spinor chains in four dimensions: the Clifford-basis expansion, Fierz rearrangement and
Sirlin's relations (FeynCalc's SirlinSimplify).

    clifford_decompose(G)        G = sum_C c_C Gamma_C over 1, gamma^a, sigma^{ab}, gamma^a g5, g5 (4 dimensions)
    fierz(M, anticommuting=False)    [ubar1 G_A u2][ubar3 G_B u4] -> sum_C w_C [ubar1 Gamma_C u4][ubar3 G_B Gamma~_C G_A u2]
    sirlin(M)                    both chains expanded in the Clifford basis and their coefficients contracted, so
                                 contracted index pairs between the chains disappear (Sirlin's identities)

The rearrangement is the completeness relation of the 16 basis matrices,
    delta_ii' delta_ll' = sum_C w_C (Gamma_C)_il (Gamma~_C)_l'i'
(w = 1/4 for 1, gamma^a, gamma^a g5, g5 and 1/8 for sigma^{ab}, the dual matrices Gamma~ = 1, gamma_a, sigma_ab,
g5 gamma_a, g5), so no basis has to be chosen for the second chain.  Field operators anticommute: then the
rearranged product has an extra minus sign (anticommuting=True).  Only valid in four dimensions.
"""
from sage.all import SR, Integer, Permutation

from .qft import tensors as T
from .qft.objects import DiracExpr, Chain, one, gamma, gamma5, sigma


_n = [0]


def _idx():
    _n[0] += 1
    return SR(T.lorentz_indices("fzi%d" % _n[0]))


def _basis():
    """[(Gamma_C, dual Gamma~_C, weight)] with fresh summed indices."""
    a, b = _idx(), _idx()
    c, d = _idx(), _idx()
    return [(one(), one(), Integer(1) / 4),
            (gamma(a), gamma(a), Integer(1) / 4),
            (sigma(b, c), sigma(b, c), Integer(1) / 8),
            (gamma(d) * gamma5(), gamma5() * gamma(d), Integer(1) / 4),
            (gamma5(), gamma5(), Integer(1) / 4)]


def _matrix(ch):
    return DiracExpr([(SR(1), Chain(None, ch.factors, None), ())])


def _attach(G, left, right):
    """The bilinear left G right as a list of (coefficient, Chain) terms."""
    return [(c, Chain(left, o.factors if o is not None else (), right)) for c, o, cl in G.terms]


def clifford_decompose(G):
    """The 4 x 4 matrix G (a DiracExpr without spinors) as a list of (coefficient, basis matrix): the coefficient
    is w_C Tr[G Gamma~_C] with the indices of the basis matrix free (summed with it)."""
    from .qft.operations import dirac_trace
    out = []
    for B, Bd, w in _basis():
        c = dirac_trace(G * Bd)
        if not SR(c).is_trivial_zero():
            out.append((w * c, B))
    return out


def _bilinears(M):
    for c, o, cl in M.terms:
        if o is not None or len(cl) != 2:
            raise ValueError("a product of exactly two spinor chains is needed, as ubar(p1)*G*u(p2) * ubar(p3)*G'*u(p4)")
        yield c, cl[0], cl[1]


def _tidy(M):
    """Summed indices renamed s1, s2, ... in the order they first appear (chains first, then the coefficient),
    epsilon tensors with their arguments in a fixed order (with the sign of the permutation), equal terms
    combined: so the same structure written with other dummy names adds up, and cancelling terms vanish."""
    from .qft.objects import _combine, Factor
    from .qft.operations import _rename_factor
    from .qft.tensors import epsilon
    eps_op = SR(epsilon(*T.lorentz_indices("fze1 fze2 fze3 fze4"))).operator()
    out = []
    for c, o, cl in M.terms:
        order = []
        for ch in cl:
            for F in ch.factors:
                for atom, cc in F.items:
                    if atom[0] == 'i' and str(atom[1]) not in order:
                        order.append(str(atom[1]))
        for v in SR(c).variables():
            if T.is_index(v) and str(v) not in order:
                order.append(str(v))
        counts = {}
        for ch in cl:
            for F in ch.factors:
                for atom, cc in F.items:
                    if atom[0] == 'i':
                        counts[str(atom[1])] = counts.get(str(atom[1]), 0) + 1
        for v in SR(c).variables():
            if T.is_index(v):
                counts[str(v)] = counts.get(str(v), 0) + _count(SR(c), v)
        summed = [n for n in order if counts.get(n, 0) >= 2]
        sub = {SR(n): SR(T.lorentz_indices("s%d" % (k + 1))) for k, n in enumerate(summed)}
        c2 = SR(c).subs(sub)
        es = {}
        for E in c2.find(eps_op(*[SR.wild(k) for k in range(4)])):
            args = list(E.operands())
            idx = sorted(range(4), key=lambda k: str(args[k]))
            sign = Permutation([k + 1 for k in idx]).sign()
            es[E] = sign * eps_op(*[args[k] for k in idx])
        c2 = c2.subs(es) if es else c2
        chains = tuple(Chain(ch.left, tuple(_rename_factor(F, sub) for F in ch.factors), ch.right) for ch in cl)
        out.append((c2.expand(), o, chains))
    return DiracExpr(_combine(out))


def _count(e, v):
    """How often the index v occurs in the monomials of e (the largest number in one monomial)."""
    e = SR(e).expand()
    terms = e.operands() if e.operator() is not None and "add" in str(e.operator()) else [e]
    best = 0
    for t in terms:
        best = max(best, str(t).count(str(v)))
    return best


def fierz(M, anticommuting=False):
    """The Fierz rearrangement of a product of two spinor chains (the right spinors exchanged), exact in four
    dimensions.  anticommuting=True for field operators (an extra minus sign)."""
    from .qft.operations import contract
    terms = []
    s = -1 if anticommuting else 1
    for c, ch1, ch2 in _bilinears(M):
        GA, GB = _matrix(ch1), _matrix(ch2)
        for B, Bd, w in _basis():
            for c1, n1 in _attach(B, ch1.left, ch2.right):
                for c2, n2 in _attach(GB * Bd * GA, ch2.left, ch1.right):
                    terms.append((s * w * c * c1 * c2, None, (n1, n2)))
    return _tidy(contract(DiracExpr(terms), dim=4))


def sirlin(M):
    """Both chains of a product expanded in the Clifford basis, the coefficients contracted in four dimensions:
    index pairs summed between the two chains (as in [g^mu g^rho g^nu P_L] x [g_mu g_tau g_nu P_L]) are
    removed, which is what Sirlin's relations do.  FeynCalc's SirlinSimplify."""
    from .qft.operations import contract
    terms = []
    for c, ch1, ch2 in _bilinears(M):
        for k1, B1 in clifford_decompose(_matrix(ch1)):
            for k2, B2 in clifford_decompose(_matrix(ch2)):
                for c1, n1 in _attach(B1, ch1.left, ch1.right):
                    for c2, n2 in _attach(B2, ch2.left, ch2.right):
                        terms.append((c * k1 * k2 * c1 * c2, None, (n1, n2)))
    return _tidy(contract(DiracExpr(terms), dim=4))
