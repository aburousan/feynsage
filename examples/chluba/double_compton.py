"""
Double Compton scattering e(P) + gamma(K0) -> e(P') + gamma(K1) + gamma(K2) from its six Feynman
diagrams, with FORM for the traces (m = 1, exact rational kinematics).

Used by the notebooks in this folder and by section 18 of examples/feynsage_walkthrough.ipynb:

    from double_compton import dc_table, dc_msq, mandl_skyrme_X, eikonal
    tab = dc_table(pk0, pk1, pk2, k0k1, k0k2)     # P.K0, P.K1, P.K2, K0.K1, K0.K2
    dc_msq(tab)            # sum over spins and polarisations of |M|^2 / e^6
    mandl_skyrme_X(tab)    # Mandl and Skyrme's X (Chluba's thesis eq. D.1)
    eikonal(tab)           # the soft-photon factor S of the photon K2
"""
import itertools
from sage.all import SR, function, var
from feynsage import form

DC_VECTORS = ['p', 'pp', 'k0', 'k1', 'k2']
_dot = function('dot')
_m = var('m')


def dc_evaluate(expr, table):
    """Put the numbers of `table` into the scalar products of a FORM result (m = 1)."""
    return expr.substitute_function(_dot, lambda a, b: table[tuple(sorted((str(a), str(b))))]).subs({_m: 1}).expand()


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
    ph = {'a': ('-k0', 'mu1'), 'b': ('k1', 'mu2'), 'c': ('k2', 'mu3')}     # photons as outgoing momenta
    chains = []
    for o in itertools.permutations('abc'):
        (q1, i1), (q2, i2), (q3, i3) = (ph[x] for x in o)
        chains.append(([i1, 'pp + %s + m' % q1, i2, 'p - (%s) + m' % q3, i3],
                       _square('pp + %s' % q1, table) * _square('p - (%s)' % q3, table)))
    tot = 0
    for c1, d1 in chains:
        for c2, d2 in chains:
            tr = form.dirac_trace(['pp + m'] + c1 + ['p + m'] + c2[::-1], vectors=DC_VECTORS)
            tot += dc_evaluate(tr, table) / (d1 * d2)
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
