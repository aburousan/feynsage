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

# `from feynsage.oneloop import *` gives only these; re-exporting Sage's gamma would hide feynsage's gamma(mu)
__all__ = ['D', 'eps', 'tadpole', 'G', 'bubble_massless', 'triangle_onshell', 'bubble_equal_mass',
           'bubble_one_mass', 'expand_eps']

D, eps = var('D eps')


def tadpole(n, m2):
    """The Euclidean tadpole Int [dl] (l^2 + m^2)^(-n) = Gamma(n - D/2)/Gamma(n) (m^2)^(D/2 - n), with
    [dl] = d^D l/pi^(D/2).  m2 is the squared mass, n any power (a number or a symbol).
        tadpole(1, m2)     # Gamma(1 - D/2) (m^2)^(D/2 - 1)"""
    return gamma(n - D/2)/gamma(n) * m2**(D/2 - n)


def G(a, b):
    """The coefficient of the massless bubble: Int [dl] (l^2)^(-a) ((l - p)^2)^(-b) = G(a, b) (p^2)^(D/2 - a - b),
    G(a, b) = Gamma(a + b - D/2) Gamma(D/2 - a) Gamma(D/2 - b)/(Gamma(a) Gamma(b) Gamma(D - a - b)).
    Euclidean, [dl] = d^D l/pi^(D/2)."""
    return gamma(a + b - D/2)*gamma(D/2 - a)*gamma(D/2 - b)/(gamma(a)*gamma(b)*gamma(D - a - b))


def bubble_massless(a, b, p2):
    """The massless bubble Int [dl] (l^2)^(-a) ((l - p)^2)^(-b) = G(a, b) (p^2)^(D/2 - a - b), Euclidean p^2 = p2.
        bubble_massless(1, 1, p2)"""
    return G(a, b) * p2**(D/2 - a - b)


def triangle_onshell(a1, a2, a3, Q2):
    """The massless triangle with two on-shell legs (p1^2 = p2^2 = 0) and Q = p1 + p2, Q^2 = Q2 (Euclidean),
    lines with powers a1, a2, a3:
        (Q^2)^(D/2 - a) Gamma(a - D/2) Gamma(D/2 - a1 - a3) Gamma(D/2 - a2 - a3)/(Gamma(a1) Gamma(a2) Gamma(D - a)),
    a = a1 + a2 + a3.  Lines 1 and 2 meet at the vertex where Q enters, line 3 joins the two on-shell legs (at
    a3 = 0 it is the bubble G(a1, a2) of Q).  Measure [dl] = d^D l/pi^(D/2)."""
    a = a1 + a2 + a3
    return Q2**(D/2 - a) * gamma(a - D/2)*gamma(D/2 - a1 - a3)*gamma(D/2 - a2 - a3) / (gamma(a1)*gamma(a2)*gamma(D - a))


def bubble_equal_mass(p2, m2):
    """The equal-mass bubble Int [dl] 1/((l^2 + m^2)((l - p)^2 + m^2)) (Euclidean, powers 1 and 1) as a Gauss
    function: Gamma(2 - D/2) (m^2)^(D/2 - 2) 2F1(2 - D/2, 1; 3/2; -p^2/(4 m^2)).  p2 = p^2, m2 = m^2."""
    x = p2/(4*m2)
    return gamma(2 - D/2) * m2**(D/2 - 2) * hypergeometric([2 - D/2, 1], [SR(3)/2], -x)


def bubble_one_mass(a, b, q2, m2):
    """The bubble with one massless line (power a) and one massive line (power b, mass^2 m2), Euclidean q^2 = q2:
        Gamma(a + b - D/2) Gamma(D/2 - a)/(Gamma(b) Gamma(D/2)) (m^2)^(D/2 - a - b) 2F1(a + b - D/2, a; D/2; -q^2/m^2).
    Measure [dl] = d^D l/pi^(D/2)."""
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
