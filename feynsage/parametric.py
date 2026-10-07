r"""
Feynman parametrisation of a whole family with symbolic powers.

    feynman_parametrize(fam)                 powers a1, ..., at
    feynman_parametrize(fam, (1, 1, 1, 1))

returns (prefactor, integrand, variables): the integral of the family with powers a, with the measure
Int prod_r d^D l_r / pi^(D/2), equals

    prefactor * Int prod dx_i delta(1 - sum x_i) integrand
    integrand = prod x_i^(a_i - 1) U^(a - (L+1) D/2) F^(L D/2 - a),  a = sum a_i
    prefactor = Gamma(a - L D/2)/prod Gamma(a_i)          (Euclidean)
              = (-1)^a i^L Gamma(a - L D/2)/prod Gamma(a_i)  (Minkowski, propagators k^2 - m^2 + i0)

U and F are those of fam.UF() (in Minkowski space F already carries the sign of the parametric formula).
The delta function may keep any non-empty subset of the x_i (Cheng and Wu).
"""
from sage.all import SR, gamma, I, prod


def feynman_parametrize(fam, powers=None, dim=None):
    D = SR.var('D') if dim is None else SR(dim)
    U, F = fam.UF()
    xs = [SR(str(v)) for v in U.parent().gens()]
    U, F = SR(U), SR(F)
    a = [SR.var('a%d' % (i + 1)) for i in range(fam.t)] if powers is None else [SR(p) for p in powers]
    A, L = sum(a), len(fam.loops)
    integrand = prod(x**(ai - 1) for x, ai in zip(xs, a) if ai != 1) * U**(A - (L + 1)*D/2) * F**(L*D/2 - A)
    pref = gamma(A - L*D/2) / prod(gamma(ai) for ai in a)
    if not fam.kin.euclidean:
        pref = (-1)**A * I**L * pref
    return pref, integrand, xs
