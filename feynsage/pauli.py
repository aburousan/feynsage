r"""
Pauli matrices in three dimensions (FeynCalc's PauliTrace, PauliSimplify, PauliSigmaCombine), for
nonrelativistic expansions (bound states, NRQCD).

    pauli_simplify(f1, f2, ...)      the product sigma.f1 sigma.f2 ... as  S + V^c sigma^c
    pauli_trace(f1, f2, ...)         Tr[sigma.f1 sigma.f2 ...] = 2 S

A factor f is a Cartesian index (cartesian_indices("i j")) for sigma^i, or a 3-vector (vectors3("a b")) for
sigma.a.  The product is reduced one factor at a time with  sigma^a sigma^b = delta^{ab} + i eps^{abc} sigma^c,
and the coefficients are contracted in three dimensions: delta substitution, eps^{abc} eps^{ade} = delta^{bd}
delta^{ce} - delta^{be} delta^{cd}, a^c b^c = a.b.  Results use delta3(i, j), eps3(i, j, k), comp3(a, i),
dot3(a, b), cross products as eps3 contracted with comp3.
"""
from sage.all import SR, I, function, Permutation

_KIND = {}
DELTA3 = function("delta3", nargs=2, print_latex_func=lambda self, a, b: r"\delta^{%s %s}" % (a, b))
EPS3 = function("eps3", nargs=3, print_latex_func=lambda self, a, b, c: r"\epsilon^{%s %s %s}" % (a, b, c))
COMP3 = function("comp3", nargs=2, print_latex_func=lambda self, a, i: r"%s^{%s}" % (a, i))
DOT3 = function("dot3", nargs=2, print_latex_func=lambda self, a, b: r"\vec{%s}\cdot\vec{%s}" % (a, b))


def cartesian_indices(names):
    """cartesian_indices("i j k") -> Cartesian (three-dimensional) indices."""
    out = []
    for n in names.replace(",", " ").split():
        _KIND[n] = "index"
        out.append(SR.var(n))
    return out[0] if len(out) == 1 else tuple(out)


def vectors3(names):
    """vectors3("a b") -> three-vectors a, b for sigma.a, comp3(a, i), dot3(a, b)."""
    out = []
    for n in names.replace(",", " ").split():
        _KIND[n] = "vector"
        out.append(SR.var(n))
    return out[0] if len(out) == 1 else tuple(out)


_n = [0]


def _fresh():
    _n[0] += 1
    return cartesian_indices("pz%d" % _n[0])


def _is_index(x):
    return _KIND.get(str(x)) == "index"


def _monomials(e):
    e = SR(e).expand()
    if e.is_zero():
        return []
    op = e.operator()
    return e.operands() if op is not None and "add" in getattr(op, "__name__", str(op)) else [e]


def _factors(m):
    op = m.operator()
    if op is not None and "mul" in getattr(op, "__name__", str(op)):
        out = []
        for f in m.operands():
            if f.operator() is not None and "pow" in str(f.operator()) and f.operands()[1].is_integer() and f.operands()[1] > 0:
                out += [f.operands()[0]] * int(f.operands()[1])
            else:
                out.append(f)
        return out
    if op is not None and "pow" in str(op) and m.operands()[1].is_integer() and m.operands()[1] > 0:
        return [m.operands()[0]] * int(m.operands()[1])
    return [m]


def _contract_mono(fs):
    """Contract repeated Cartesian indices in a list of factors; returns a list of (coefficient, factors)."""
    for k, F in enumerate(fs):
        if F.operator() is DELTA3:
            a, b = F.operands()
            if a == b:
                return [(SR(3), fs[:k] + fs[k + 1:])]
            for x, y in ((a, b), (b, a)):
                if _is_index(x):
                    rest = fs[:k] + fs[k + 1:]
                    if any(x in G.operands() for G in rest if G.operator() in (DELTA3, EPS3, COMP3)):
                        return [(SR(1), [G.subs({x: y}) for G in rest])]
    # comp3(a, c) comp3(b, c) -> dot3(a, b)
    comps = [(k, F) for k, F in enumerate(fs) if F.operator() is COMP3]
    for (k1, F1) in comps:
        for (k2, F2) in comps:
            if k1 < k2 and F1.operands()[1] == F2.operands()[1]:
                a, b = sorted([F1.operands()[0], F2.operands()[0]], key=str)
                rest = [G for j, G in enumerate(fs) if j not in (k1, k2)]
                return [(SR(1), rest + [DOT3(a, b)])]
    # eps eps with a common index
    epses = [(k, F) for k, F in enumerate(fs) if F.operator() is EPS3]
    for (k1, E1) in epses:
        for (k2, E2) in epses:
            if k1 >= k2:
                continue
            common = [x for x in E1.operands() if x in E2.operands()]
            if not common:
                continue
            c = common[0]
            a1, s1 = _bring_first(E1.operands(), c)
            a2, s2 = _bring_first(E2.operands(), c)
            b, cc = a1[1], a1[2]
            d, e = a2[1], a2[2]
            rest = [G for j, G in enumerate(fs) if j not in (k1, k2)]
            return [(SR(s1 * s2), rest + [DELTA3(b, d), DELTA3(cc, e)]),
                    (SR(-s1 * s2), rest + [DELTA3(b, e), DELTA3(cc, d)])]
    # eps with a repeated index, or eps(a, b, c) with two equal entries
    for k, E in epses:
        ops = E.operands()
        if len(set(str(x) for x in ops)) < 3:
            return [(SR(0), [])]
    return None


def _bring_first(ops, c):
    ops = list(ops)
    j = ops.index(c)
    order = [j] + [i for i in range(3) if i != j]
    sign = Permutation([i + 1 for i in order]).sign()
    return [ops[i] for i in order], sign


def contract3(e):
    """Contract the repeated Cartesian indices of a Sage expression in delta3, eps3, comp3, dot3."""
    work = [(SR(1), _factors(m)) for m in _monomials(e)]
    out = SR(0)
    while work:
        c, fs = work.pop()
        r = _contract_mono(fs)
        if r is None:
            t = c
            for F in fs:
                t *= F
            out += t
        else:
            for cc, gs in r:
                if not cc.is_zero():
                    work.append((c * cc, gs))
    return _canon(out.expand())


def _canon(e):
    """delta3 and dot3 with sorted arguments, eps3 sorted with the sign of the permutation."""
    w = [SR.wild(k) for k in range(3)]
    sub = {}
    for F in e.find(DELTA3(w[0], w[1])):
        a, b = F.operands()
        if str(b) < str(a):
            sub[F] = DELTA3(b, a)
    for F in e.find(DOT3(w[0], w[1])):
        a, b = F.operands()
        if str(b) < str(a):
            sub[F] = DOT3(b, a)
    for F in e.find(EPS3(w[0], w[1], w[2])):
        ops = list(F.operands())
        order = sorted(range(3), key=lambda k: str(ops[k]))
        sub[F] = Permutation([k + 1 for k in order]).sign() * EPS3(*[ops[k] for k in order])
    return e.subs(sub) if sub else e


def _slot(f, c):
    """The coefficient of sigma^c in the factor f: delta3(i, c) for an index i, comp3(a, c) for a vector a."""
    f = SR(f)
    if _KIND.get(str(f)) == "vector":
        return COMP3(f, c)
    if _is_index(f):
        return DELTA3(f, c)
    raise TypeError("a Pauli factor is a Cartesian index or a three-vector, got %s" % f)


def pauli_simplify(*factors):
    """sigma.f1 sigma.f2 ... = S + V^c sigma^c: returns (S, V, c) with V a Sage expression with the free
    Cartesian index c (the product of no factors is the unit matrix: S = 1)."""
    S, c = SR(1), _fresh()
    V = SR(0)
    for f in factors:
        b = _fresh()
        Fb = _slot(f, b)                      # sigma.f = Fb sigma^b (b summed)
        d = _fresh()
        # (S + V^c sigma^c) sigma^b = V^b + (S delta^{bd} + i V^c eps^{cbd}) sigma^d, times Fb, b summed
        newS = contract3(V.subs({c: b}) * Fb)
        newV = contract3((S * DELTA3(b, d) + I * V * EPS3(c, b, d)) * Fb)
        S, V, c = newS, newV, d
    return S, V, c


def pauli_trace(*factors):
    """Tr[sigma.f1 sigma.f2 ...] = 2 S (Tr 1 = 2, Tr sigma = 0)."""
    S, V, c = pauli_simplify(*factors)
    return (2 * S).expand()
