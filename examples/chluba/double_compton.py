"""
Double Compton scattering e(P) + gamma(K0) -> e(P') + gamma(K1) + gamma(K2) from its six Feynman
diagrams, with FORM for the traces (m = 1, exact rational kinematics).  The chains are written as in
the book, gamma(mu1) * (slash(p' + q) + m) * ..., and all 36 traces go to FORM in one call.

Used by the notebooks in this folder and by section 18 of examples/feynsage_walkthrough.ipynb:

    from double_compton import dc_table, dc_msq, mandl_skyrme_X, eikonal
    tab = dc_table(pk0, pk1, pk2, k0k1, k0k2)     # P.K0, P.K1, P.K2, K0.K1, K0.K2
    dc_msq(tab)            # sum over spins and polarisations of |M|^2 / e^6
    mandl_skyrme_X(tab)    # Mandl and Skyrme's X (Chluba's thesis eq. D.1)
    eikonal(tab)           # the soft-photon factor S of the photon K2
"""
import itertools
from sage.all import SR, var
from feynsage import momenta, lorentz_indices, gamma, slash, dirac_trace, dot

DC_VECTORS = ['p', 'pp', 'k0', 'k1', 'k2']
P, PP, K0, K1, K2 = momenta("p pp k0 k1 k2")
MU1, MU2, MU3 = lorentz_indices("mu1 mu2 mu3")
_m = var('m')


def dc_evaluate(expr, table):
    """Put the numbers of `table` into the scalar products of a trace (m = 1)."""
    vec = dict(zip(DC_VECTORS, (P, PP, K0, K1, K2)))
    return SR(expr).subs({dot(vec[a], vec[b]): val for (a, b), val in table.items()}).subs({_m: 1}).expand()


def dc_table(pk0, pk1, pk2, k0k1, k0k2):
    """All scalar products from five of them. P' = P + K0 - K1 - K2 on shell fixes K1.K2."""
    k1k2 = -(pk0 - pk1 - pk2 - k0k1 - k0k2)
    t = {('p', 'p'): 1, ('k0', 'k0'): 0, ('k1', 'k1'): 0, ('k2', 'k2'): 0, ('k0', 'p'): pk0, ('k1', 'p'): pk1,
         ('k2', 'p'): pk2, ('k0', 'k1'): k0k1, ('k0', 'k2'): k0k2, ('k1', 'k2'): k1k2}
    d = lambda a, b: t[tuple(sorted((a, b)))]
    pp = {'p': 1, 'k0': 1, 'k1': -1, 'k2': -1}
    for v in ['p', 'k0', 'k1', 'k2']:
        t[tuple(sorted(('pp', v)))] = sum(cc * d(u, v) for u, cc in pp.items())
    t[('pp', 'pp')] = sum(c1 * c2 * d(u1, u2) for u1, c1 in pp.items() for u2, c2 in pp.items())
    return {tuple(sorted(k)): v for k, v in t.items()}


def _square(expr, table):
    """(momentum)^2 - m^2 for a combination such as 'pp + -k0' or 'p - (k1)'."""
    e = SR(expr.replace('pp', 'PPx'))
    co = {v: e.coefficient(SR.var('PPx' if v == 'pp' else v)) for v in DC_VECTORS}
    return sum(co[a] * co[b] * table[tuple(sorted((a, b)))] for a in DC_VECTORS for b in DC_VECTORS) - 1


def dc_msq(table):
    """Sum over spins and polarisations of |M|^2 / e^6, all six orders of the three photons."""
    ph = {'a': (-K0, MU1, '-k0'), 'b': (K1, MU2, 'k1'), 'c': (K2, MU3, 'k2')}   # photons as outgoing momenta
    chains = []
    for o in itertools.permutations('abc'):
        (q1, i1, s1), (q2, i2, s2), (q3, i3, s3) = (ph[x] for x in o)
        line = gamma(i1) * (slash(PP + q1) + _m) * gamma(i2) * (slash(P - q3) + _m) * gamma(i3)
        back = gamma(i3) * (slash(P - q3) + _m) * gamma(i2) * (slash(PP + q1) + _m) * gamma(i1)
        chains.append((line, back, _square('pp + %s' % s1, table) * _square('p - (%s)' % s3, table)))
    # all 36 traces Tr[(p'/ + m) chain_1 (p/ + m) reversed chain_2] in one FORM run
    exprs = [(slash(PP) + _m) * c1 * (slash(P) + _m) * b2 for c1, _, _ in chains for _, b2, _ in chains]
    dens = [d1 * d2 for _, _, d1 in chains for _, _, d2 in chains]
    traces = dirac_trace(exprs)
    tot = sum(dc_evaluate(tr, table) / d for tr, d in zip(traces, dens))
    return -tot                                       # (-g) for each of the three photons


def mandl_skyrme_X(table):
    """Mandl and Skyrme's X, with |M|^2 = e^6 X in Chluba's thesis (eq. D.1)."""
    d = lambda a, b: table[tuple(sorted((a, b)))]
    k0, k1, k2 = -d('p', 'k0'), d('p', 'k1'), d('p', 'k2')
    k0p, k1p, k2p = d('pp', 'k0'), -d('pp', 'k1'), -d('pp', 'k2')
    a = 1 / k0 + 1 / k1 + 1 / k2
    b = 1 / k0p + 1 / k1p + 1 / k2p
    c = 1 / (k0 * k0p) + 1 / (k1 * k1p) + 1 / (k2 * k2p)
    x = k0 + k1 + k2
    z = k0 * k0p + k1 * k1p + k2 * k2p
    A = k0 * k1 * k2
    B = k0p * k1p * k2p
    rho = k0 / k0p + k0p / k0 + k1 / k1p + k1p / k1 + k2 / k2p + k2p / k2
    return (2 * (a * b - c) * ((a + b) * (2 + x) - (a * b - c) - 8) - 2 * x * (a ** 2 + b ** 2)
            - 2 * (a * b + c * (1 - x)) * rho - 8 * c
            + 4 * x / (A * B) * ((A + B) * (1 + x) + x ** 2 * (1 - z) + 2 * z - (a * A + b * B) * (2 + (1 - x) * z / x)))


def eikonal(table):
    """The soft factor S = 2 P.P'/((P.K2)(P'.K2)) - m^2/(P.K2)^2 - m^2/(P'.K2)^2 of the photon K2."""
    d = lambda a, b: table[tuple(sorted((a, b)))]
    return 2 * d('p', 'pp') / (d('p', 'k2') * d('pp', 'k2')) - 1 / d('p', 'k2') ** 2 - 1 / d('pp', 'k2') ** 2
