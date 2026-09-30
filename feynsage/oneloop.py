r"""
Closed one-loop results (Euclidean, measure [dl] = d^D l / pi^(D/2)), all derived in
the lecture notes, plus an epsilon-expansion helper.

- tadpole(n, m2):                Int [dl] (l^2+m^2)^(-n)
- bubble_massless(a, b, p2):     Int [dl] (l^2)^(-a) ((l-p)^2)^(-b)
- triangle_onshell(a1, a2, a3, Q2): massless on-shell triangle with arbitrary powers
- bubble_equal_mass(p2, m2):     equal-mass bubble as a 2F1
- bubble_one_mass(a, b, q2, m2): massless line (power a) and massive line (power b)
"""
from sage.all import SR, ZZ, var, gamma, hypergeometric, exp, euler_gamma, pi

D, eps = var('D eps')


def tadpole(n, m2):
    return gamma(n - D/2)/gamma(n) * m2**(D/2 - n)


def G(a, b):
    return gamma(a + b - D/2)*gamma(D/2 - a)*gamma(D/2 - b)/(gamma(a)*gamma(b)*gamma(D - a - b))


def bubble_massless(a, b, p2):
    return G(a, b) * p2**(D/2 - a - b)


def triangle_onshell(a1, a2, a3, Q2):
    a = a1 + a2 + a3
    return Q2**(D/2 - a) * gamma(a - D/2)*gamma(D/2 - a1 - a3)*gamma(D/2 - a2 - a3) / (gamma(a1)*gamma(a2)*gamma(D - a))


def bubble_equal_mass(p2, m2):
    x = p2/(4*m2)
    return gamma(2 - D/2) * m2**(D/2 - 2) * hypergeometric([2 - D/2, 1], [SR(3)/2], -x)


def bubble_one_mass(a, b, q2, m2):
    z = q2/m2
    return (gamma(a + b - D/2)*gamma(D/2 - a)/(gamma(b)*gamma(D/2)) * m2**(D/2 - a - b)
            * hypergeometric([a + b - D/2, a], [D/2], -z))


def _regularise_gamma(e):
    """Rewrite Gamma(-k + c eps) = Gamma(1 + c eps) / [(-k + c eps)(-k + 1 + c eps)...(c eps)]
    (from Gamma(z+1) = z Gamma(z) used k+1 times), so every Gamma left is regular at eps = 0."""
    from sage.all import prod
    from sage.symbolic.operators import add_vararg, mul_vararg
    op = e.operator()
    if op is None:
        return e
    args = [_regularise_gamma(a) for a in e.operands()]
    if op == gamma:
        z = args[0]
        z0 = z.subs(eps=0)
        if z0 in ZZ and z0 <= 0:
            k = -ZZ(z0)
            return gamma(z + k + 1) / prod(z + j for j in range(k + 1))
        return gamma(z)
    return op(*args)


def expand_eps(expr, order=0, loops=1):
    """Laurent expansion at D = 4 - 2 eps, with the conventional factor e^(L eps gamma_E)."""
    e = (exp(loops*eps*euler_gamma) * expr).subs(D=4 - 2*eps)
    e = _regularise_gamma(e)
    return e.series(eps, order + 1).truncate().simplify_full()
