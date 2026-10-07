r"""
Counting master integrals the way Lee and Pomeransky do: in a sector (the lines that are present) count
the points where every derivative of G = U + F vanishes, with every x_i of the sector non-zero and G
non-zero.  The count is the dimension of Q[x]/I, where I is the ideal of the derivatives saturated by
x_1 ... x_n G (so solutions on x_i = 0 or G = 0 are removed), i.e. the critical points with multiplicity.
Give numbers for the invariants: the count is for generic kinematics.

    critical_points(U + F, [x1, x2])                         # one sector
    master_count(fam, sectors, values={"s": -11/5, ...})     # several sectors of a family
"""
from sage.all import SR, QQ, PolynomialRing, prod


def critical_points(G, xs):
    """Number of critical points of G in the variables xs with all xs != 0 and G != 0 (with multiplicity)."""
    xs = [SR(x) for x in xs]
    if not xs:                                  # no variables left: one point, if G is not zero there
        return 0 if SR(G).is_zero() else 1
    R = PolynomialRing(QQ, len(xs), [str(x) for x in xs], order='degrevlex')    # multivariate even for one x
    g = R(SR(G))
    I = R.ideal([g.derivative(v) for v in R.gens()])
    J, _ = I.saturation(R.ideal([prod(R.gens()) * g]))
    if J.dimension() > 0:
        return float('inf')                     # a whole curve of critical points: this count does not apply
    return J.vector_space_dimension()


def master_count(fam, sectors, values=None):
    """{sector: number of critical points} for the sectors (tuples of 0/1) of a family, at the given values
    of the invariants ({name: number})."""
    U, F = fam.UF()
    xs = [SR(str(v)) for v in U.parent().gens()]
    G = SR(U) + SR(F)
    if values:
        G = G.subs({SR.var(k): QQ(v) for k, v in values.items()})
    out = {}
    for sec in sectors:
        sub = {x: 0 for x, on in zip(xs, sec) if not on}
        out[tuple(sec)] = critical_points(G.subs(sub), [x for x, on in zip(xs, sec) if on])
    return out
