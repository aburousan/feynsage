r"""
Sector decomposition of Feynman-parameter integrals, with the expansion in eps.

    secs = sector_decompose([(U, 2*eps - 1), (F, -1 - eps)], [x1, x2, x3])
    integrate_sectors(secs, eps, order=0, values={s: -1})

The integral is  Int_0^oo prod dx_i delta(1 - sum x_i) prod_j P_j^(e_j)  (prod x_i^(a_i - 1) may be included
as polynomials x_i with their exponents).  Primary sectors: x_l is the largest parameter, the others are
x_i = x_l t_i with 0 <= t_i <= 1, and by homogeneity x_l drops out (x_l = 1).  Then, while some polynomial
vanishes when a set S of the t's is zero, the sector is split |S| ways: t_i = t_k t_i' for i in S, i != k,
which gives a Jacobian t_k^(|S| - 1) and a power of t_k that comes out of every polynomial.  At the end
every polynomial is non-zero at t = 0 and the singularities sit in the monomials t_i^(-1 + c_i eps).

Expansion: t^(-1 + c eps) = delta(t)/(c eps) + sum_k (c eps)^k/k! [log^k t/t]_+, so every sector becomes
a sum of finite integrals over the unit cube.  integrate_sectors does them numerically
(tanh-sinh quadrature, about 12 digits for each coefficient of eps), as programs like pySecDec do;
exact=True asks Sage for the integrals (it manages only simple cases, logarithm branches often defeat it).
"""
import itertools
from sage.all import SR, QQ, prod, factorial, log, integrate, numerical_integral


class Sector:
    """prod_i t_i^(mono_i) prod_j P_j^(e_j) over the unit cube in the variables ts."""

    def __init__(self, label, ts, mono, polys):
        self.label, self.ts, self.mono, self.polys = label, list(ts), dict(mono), list(polys)

    def _latex_(self):
        from sage.all import latex
        m = ' '.join('%s^{%s}' % (latex(t), latex(e)) for t, e in self.mono.items() if not SR(e).is_trivial_zero())
        p = ' '.join(r'\left(%s\right)^{%s}' % (latex(P), latex(e)) for P, e in self.polys if not (SR(P) - 1).is_trivial_zero())
        return r'\text{%s:}\ %s\, %s' % (self.label, m or '1', p)

    def __repr__(self):
        m = '*'.join('%s^(%s)' % (t, e) for t, e in self.mono.items() if not SR(e).is_trivial_zero())
        p = '*'.join('(%s)^(%s)' % (P, e) for P, e in self.polys)
        return 'Sector %s: %s * %s' % (self.label, m or '1', p)


def _lowest(P, t):
    """Lowest power of t in the expanded polynomial P."""
    P = SR(P).expand()
    terms = P.operands() if P.operator() is not None and 'add' in str(P.operator()) else [P]
    return min(int(SR(term).degree(t)) for term in terms)


def _zero_set(P, ts):
    """A smallest set S of the variables such that P vanishes when all t in S are 0 (or None)."""
    P = SR(P).expand()
    if not P.subs({t: 0 for t in ts}).is_trivial_zero():
        return None
    for k in range(1, len(ts) + 1):
        for S in itertools.combinations(ts, k):
            if P.subs({t: 0 for t in S}).expand().is_trivial_zero():
                return list(S)
    return list(ts)


def sector_decompose(polys, xs):
    """polys: [(P, exponent)], xs: the Feynman parameters.  Returns a list of Sectors."""
    xs = [SR(x) for x in xs]
    todo = []
    for l, xl in enumerate(xs):                                  # primary sectors
        ts = [x for x in xs if x is not xl]
        sub = {xl: 1}
        ps = [(SR(P).subs(sub).expand(), SR(e)) for P, e in polys]
        todo.append(Sector(str(l + 1), ts, {t: SR(0) for t in ts}, ps))
    done = []
    while todo:
        sec = todo.pop()
        S = None
        for P, e in sec.polys:
            S = _zero_set(P, sec.ts)
            if S:
                break
        if not S:
            done.append(sec)
            continue
        for k in S:                                              # t_i = t_k t_i' for i in S, i != k
            sub = {t: t * k for t in S if t is not k}
            mono = dict(sec.mono)
            mono[k] = mono[k] + (len(S) - 1) + sum(sec.mono[t] for t in S if t is not k)
            ps = []
            for P, e in sec.polys:
                Q = SR(P).subs(sub).expand()
                n = _lowest(Q, k)
                if n:
                    Q = (Q / k**n).expand()
                    mono[k] = mono[k] + n * e
                ps.append((Q, e))
            todo.append(Sector(sec.label + str(sec.ts.index(k) + 1), sec.ts, mono, ps))
    return done


def _expand_sector(sec, eps, order):
    """The sector as a sum of (coefficient, finite integrand) pairs up to eps^order."""
    poles, regs = [], []
    for t in sec.ts:
        e = SR(sec.mono[t]).expand()
        a0, c = e.subs({eps: 0}), e.coefficient(eps)
        if a0 == -1:
            poles.append((t, c))
        elif a0 < -1:
            raise ValueError("power of %s below -1 in sector %s: not a logarithmic singularity" % (t, sec.label))
        else:
            regs.append(t)
    g = prod(t**sec.mono[t] for t in regs) * prod(P**e for P, e in sec.polys)
    out = []
    for r in range(len(poles) + 1):
        for S in itertools.combinations(range(len(poles)), r):
            pref = prod(1 / (poles[i][1] * eps) for i in S)
            h = g.subs({poles[i][0]: 0 for i in S})
            for i in range(len(poles)):                          # plus prescription on the others
                if i in S:
                    continue
                t, c = poles[i]
                h = (h - h.subs({t: 0})) * t**(c * eps - 1)
            # expand in eps to the order needed (pref has |S| poles)
            ser = h.series(eps, order + r + 1).truncate()
            for k in range(-len(poles), order + r + 1):
                ck = ser.coefficient(eps, k)
                if not SR(ck).is_trivial_zero():
                    out.append((pref * eps**k, ck))
    return out


def integrate_sectors(sectors, eps, order=0, values=None, exact=False):
    """Sum over the sectors of the integral over the unit cube, expanded to eps^order.  values: numbers for
    the invariants.  exact=False integrates numerically (each coefficient to about 15 digits)."""
    values = values or {}
    total = SR(0)
    for sec in sectors:
        for pref, f in _expand_sector(sec, eps, order):
            f = SR(f).subs(values)
            ts = [t for t in sec.ts if f.has(t)]
            if exact:
                from sage.all import assume, forget
                facts = [t > 0 for t in sec.ts] + [t < 1 for t in sec.ts]
                assume(*facts)
                try:
                    val = f
                    for t in ts:
                        val = integrate(val, t, 0, 1)
                finally:
                    forget(*facts)
            else:
                val = _numint(f, ts)
            total += pref * val
    total = total.expand()
    return sum((c * eps**k for c, k in total.coefficients(eps)), SR(0))


def _numint(f, ts):
    import mpmath
    from sage.ext.fast_callable import fast_callable
    from sage.all import RDF
    if not ts:
        return SR(f).n()
    fc = fast_callable(SR(f), vars=ts, domain=RDF)
    mpmath.mp.dps = 15
    val = mpmath.quad(lambda *a: float(fc(*[float(x) for x in a])), *[[0, 1]] * len(ts))
    return SR(float(val))
