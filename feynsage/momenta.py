r"""
Momenta as formal linear combinations, and their scalar products.

A momentum is a dict  {name: coefficient}  with rational coefficients, e.g.
{'l1': 1, 'p1': -1}  for  l1 - p1.  Names that are loop momenta and names that
are external momenta are kept apart by the caller.

Scalar products of external momenta are supplied by a Kinematics object; scalar
products that involve a loop momentum stay as polynomial variables.
"""
from sage.all import QQ, PolynomialRing


def mom(**kw):
    """Convenience constructor: mom(l1=1, p1=-1) is l1 - p1."""
    return {k: QQ(v) for k, v in kw.items() if v != 0}


def add(a, b, cb=1):
    out = dict(a)
    for k, v in b.items():
        out[k] = out.get(k, 0) + cb * v
        if out[k] == 0:
            del out[k]
    return out


class Kinematics:
    r"""
    External kinematics.

    INPUT:
    - ``externals`` -- list of independent external momentum names (a basis)
    - ``invariants`` -- list of names of symbolic invariants (e.g. ['s','t','m2'])
    - ``rules`` -- dict {(pi, pj): expression-string} giving pi.pj in terms of invariants.
      Missing pairs are an error when needed.
    - ``euclidean`` -- if True, squares are Euclidean (propagator k^2 + m^2);
      otherwise Minkowski (propagator k^2 - m^2).
    """

    def __init__(self, externals, invariants, rules, euclidean=False, extra_vars=()):
        self.externals = list(externals)
        self.invariants = list(invariants)
        self.euclidean = euclidean
        self.R = PolynomialRing(QQ, ['d'] + self.invariants + list(extra_vars))
        try:                         # print the invariants the way notation() or var(..., latex_name=) asks
            from sage.all import SR, latex
            self.R._latex_names = ['{%s}' % latex(SR.var(str(g))) for g in self.R.gens()]
        except Exception:
            pass
        self.K = self.R.fraction_field()
        self.rules = {}
        for (a, b), e in rules.items():
            val = self.R(e) if isinstance(e, str) else self.R(e)
            self.rules[(a, b)] = val
            self.rules[(b, a)] = val

    def ext_sp(self, a, b):
        if (a, b) not in self.rules:
            raise KeyError("scalar product %s.%s not given" % (a, b))
        return self.rules[(a, b)]

    def square_external(self, q):
        """q.q for a purely external momentum q (dict)."""
        tot = self.R(0)
        for a, ca in q.items():
            for b, cb in q.items():
                tot += ca * cb * self.ext_sp(a, b)
        return tot
