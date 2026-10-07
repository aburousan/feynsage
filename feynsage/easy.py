r"""
Short, learning-friendly front ends.  Everything here builds the objects of the other
modules from plain strings, so a first calculation takes two or three lines.

    from feynsage import *
    kite = family(["l1", "l1 + q", "l1 + l2", "l1 + l2 + q", "l2"], kin={"q^2": 1})
    kite.info()                                # what the family is
    res = ibp_reduce(kite, ["F(1,1,1,1,1)", "F(2,2,1,2,2)"])
    res                                        # the reductions, printed
    g = diagram("kite"); g.plot()              # a named diagram, ready to draw or use

Conventions: a propagator is written as a momentum, or as (momentum, mass); masses are
masses (the family stores mass^2).  kin maps "p^2", "p.q" to expressions or numbers.
Minkowski propagators are (k^2 - m^2), euclidean=True gives (k^2 + m^2).
"""
import re
from itertools import permutations
from sage.all import SR, QQ, ZZ, Graph, gcd, lcm

from .momenta import Kinematics, mom
from .family import IntegralFamily
from .graph import FeynmanGraph
from .laporta import Reducer, weight
from ._parse import sr


# ---------------------------------------------------------------------------- parsing helpers
def _names(text):
    return set(re.findall(r'[A-Za-z_]\w*', str(text)))


def _mom(text, names):
    e = sr(str(text))
    out = {}
    for n in names:
        c = e.coefficient(SR.var(n))
        if c != 0:
            out[n] = QQ(c)
    return out


def _kin_pairs(kin):
    out = {}
    for key, val in (kin or {}).items():
        k = key.replace(' ', '')
        a, b = (k[:-2], k[:-2]) if k.endswith('^2') else k.split('.')
        out[(a, b)] = sr(val)
    return out


def kinematics(externals, kin=None, euclidean=False, extra=()):
    """Kinematics from {"p^2": "s", "p1.p2": "t/2"}; products that are not given become symbols."""
    pairs = _kin_pairs(kin)
    rules = {}
    for i, a in enumerate(externals):
        for b in externals[i:]:
            v = pairs.get((a, b), pairs.get((b, a)))
            if v is None:
                v = SR.var(a + '2' if a == b else a + b)
            rules[(a, b)] = v
    invariants = sorted({str(x) for v in rules.values() for x in v.variables()} | {str(x) for e in extra for x in sr(e).variables()})
    return Kinematics(list(externals), invariants, {k: str(v) for k, v in rules.items()}, euclidean=euclidean)


# ---------------------------------------------------------------------------- families
class Family(IntegralFamily):
    """An IntegralFamily built by family(); it remembers how it was written."""

    def info(self):
        print("Integral family %r" % self.name)
        print("  loop momenta: %s   external momenta: %s   %s metric"
              % (", ".join(self.loops), ", ".join(self.kin.externals) or "none",
                 "Euclidean" if self.kin.euclidean else "Minkowski"))
        for i, (q, m2) in enumerate(self.props):
            mom_text = " + ".join("%s*%s" % (c, a) if c != 1 else a for a, c in q.items()).replace("+ -", "- ")
            mass_text = "" if m2 == 0 else " %s %s" % ("+" if self.kin.euclidean else "-", m2)
            print("  D%d = (%s)^2%s" % (i + 1, mom_text, mass_text))
        print("  %s(a1,...,at) = Int prod_r d^d l_r prod_i D_i^(-a_i), a_i <= 0 is a numerator."
              % self.name)

    def integral(self, *a):
        """The integral with powers a as a Sage symbol, e.g. T(2) or J(1,1)."""
        from sage.symbolic.function_factory import function as _fn
        return _fn(self.name)(*[SR(x) for x in a])

    def ibp_eq(self, a=None):
        """The IBP identities at the powers a, written as equations  sum_k c_k F(...) == 0.
        Without a: the identities for symbolic powers a1, ..., at."""
        if a is None:
            return self._ibp_symbolic()
        a = tuple(a) if isinstance(a, (tuple, list)) else (a,)
        out = []
        for rel in self.ibp(a):
            lhs = sum((SR(str(c)) * self.integral(*k) for k, c in rel.items()), SR(0))
            out.append(lhs == 0)
        return out

    def _ibp_symbolic(self):
        from .ff import ibp_templates
        gens = [SR(str(g)) for g in self.kin.R.gens()]
        avars = [SR.var('a%d' % (i + 1)) for i in range(self.t)]

        def poly(dct):
            tot = SR(0)
            for exps, c in dct.items():
                term = SR(c)
                for g, e in zip(gens, exps):
                    term *= g ** e
                tot += term
            return tot
        out = []
        for tpl in ibp_templates(self):
            lhs = SR(0)
            for shift, c0, cs in tpl:
                coef = poly(c0) + sum((avars[i] * poly(c) for i, c in enumerate(cs) if c), SR(0))
                lhs += coef * self.integral(*[avars[i] + shift[i] for i in range(self.t)])
            out.append(lhs.collect_common_factors() == 0)
        return out


def family(props, kin=None, loops=None, euclidean=False, name="F", explain=False):
    r"""
    An integral family from strings.

        family(["l", "l + p"], kin={"p^2": "s"})                        # massless bubble
        family([("l", "m"), ("l - p", "m")], kin={"p^2": "s"})          # equal-mass bubble
        family(["l1", "l1 + q", "l1 + l2", "l1 + l2 + q", "l2"], kin={"q^2": 1})   # kite

    loops: the loop momenta (default: every name that is not an external momentum in kin;
    if kin is empty, names starting with l or k).  Returns a Family (an IntegralFamily with .info()).
    """
    if explain:
        print(family.__doc__)
    items = [p if isinstance(p, (tuple, list)) else (p, 0) for p in props]
    names = set().union(*[_names(q) for q, _ in items])
    ext = sorted({a for pair in _kin_pairs(kin) for a in pair})
    if loops is None:
        loops = sorted(n for n in names if n not in ext) if ext else sorted(n for n in names if n[0] in 'lk')
    elif isinstance(loops, str):
        loops = loops.replace(',', ' ').split()
    externals = sorted((names - set(loops)) | set(ext))
    masses = [sr(m) for _, m in items]
    K = kinematics(externals, kin, euclidean, extra=masses)
    fam_props = [(_mom(q, list(loops) + externals), str(sr(m) ** 2)) for q, m in items]
    f = Family(name, list(loops), K, fam_props)
    return f


# ---------------------------------------------------------------------------- graphs
def graph(edges, legs, kin=None, euclidean=False, explain=False):
    r"""
    A Feynman graph from strings.

        graph("A-B, A-B", {"A": "p", "B": "-p"}, kin={"p^2": "s"})             # bubble
        graph("A-B:m, A-B:m", {"A": "p", "B": "-p"}, kin={"p^2": "s"})         # massive bubble
        graph("L-B, L-T, B-R, T-R, T-B", {"L": "p", "R": "-p"}, kin={"p^2": "s"})  # kite

    An edge is "U-V" or "U-V:mass".  legs gives the incoming momentum at each vertex.
    The line numbers (x_1, x_2, ...) follow the order of the edges.
    """
    if explain:
        print(graph.__doc__)
    if isinstance(edges, str):
        edges = [e.strip() for e in edges.split(',') if e.strip()]
    lines, masses = [], []
    for e in edges:
        uv, _, m = e.partition(':')
        u, v = [x.strip() for x in uv.split('-')]
        m = sr(m.strip()) if m.strip() else sr(0)
        lines.append((u, v, str(m ** 2)))
        masses.append(m)
    names = sorted(set().union(*[_names(q) for q in legs.values()]) if legs else set())
    K = kinematics(names, kin, euclidean, extra=masses)
    return FeynmanGraph(lines, {v: _mom(q, names) for v, q in legs.items()}, K)


_LIBRARY = {
    "tadpole": (lambda: graph("A-A:m", {}, kin={}), "one massive line from a vertex back to itself"),
    "bubble": (lambda: graph("A-B, A-B", {"A": "p", "B": "-p"}, kin={"p^2": "s"}), "massless bubble, p^2 = s"),
    "bubble_mass": (lambda: graph("A-B:m, A-B:m", {"A": "p", "B": "-p"}, kin={"p^2": "s"}), "equal-mass bubble"),
    "triangle": (lambda: graph("Z-X, X-Y, Y-Z", {"X": "p1", "Y": "p2", "Z": "-p1-p2"},
                               kin={"p1^2": "P1Sq", "p2^2": "P2Sq", "p1.p2": "(QSq - P1Sq - P2Sq)/2"}),
                 "massless triangle with three off-shell legs"),
    "box": (lambda: graph("W-X, X-Y, Y-Z, Z-W", {"X": "p1", "Y": "p2", "Z": "-p1-p2-p3", "W": "p3"},
                          kin={"p1^2": 0, "p2^2": 0, "p3^2": 0, "p1.p2": "s/2", "p1.p3": "t/2", "p2.p3": "-(s+t)/2"}),
            "massless box, on-shell legs, s and t"),
    "sunset": (lambda: graph("A-B:m, A-B:m, A-B:m", {"A": "p", "B": "-p"}, kin={"p^2": "s"}), "two-loop sunset, three equal masses"),
    "kite": (lambda: graph("L-B, L-T, B-R, T-R, T-B", {"L": "p", "R": "-p"}, kin={"p^2": "s"}), "two-loop massless kite"),
    "phi4_two": (lambda: graph("A-B, A-C, B-C, B-C", {"A": "p3+p4", "B": "-p3", "C": "-p4"},
                               kin={"p3^2": 0, "p4^2": 0, "p3.p4": "Q2/2"}), "two-loop phi^4 diagram of the notes"),
    "vertex2": (lambda: graph("B1-C1, A-B1, A-B2, B2-C2, C1-C2, B1-B2", {"A": "p1+p2", "C1": "-p1", "C2": "-p2"},
                              kin={"p1^2": 0, "p2^2": 0, "p1.p2": "q2/2"}),
                "planar two-loop vertex, on-shell massless legs (lines x1..x6 as in the notes)"),
    "banana3": (lambda: graph("A-B:m, A-B:m, A-B:m, A-B:m", {"A": "p", "B": "-p"}, kin={"p^2": "s"}),
                "three-loop banana, four equal masses"),
    "ladder3": (lambda: graph("L-A, A-B, B-R, L-C, C-D, D-R, A-C, B-D", {"L": "p", "R": "-p"}, kin={"p^2": "s"}),
                "three-loop massless ladder propagator"),
    "triplebox": (lambda: graph("B1-T1, T1-T2, T2-T3, T3-T4, T4-B4, B4-B3, B3-B2, B2-B1, T2-B2, T3-B3",
                                {"B1": "p1", "T1": "p2", "T4": "p3", "B4": "-p1-p2-p3"},
                                kin={"p1^2": 0, "p2^2": 0, "p3^2": 0, "p1.p2": "s/2", "p2.p3": "t/2", "p1.p3": "-(s+t)/2"}),
                  "three-loop massless planar triple box, on-shell legs, s and t"),
}


def diagram(name=None):
    """A ready-made FeynmanGraph by name; diagram() lists them."""
    if name is None:
        for k, (_, text) in _LIBRARY.items():
            print("%-12s %s" % (k, text))
        return None
    return _LIBRARY[name][0]()


# ---------------------------------------------------------------------------- symmetries
def symmetries(fam):
    r"""
    Permutations of the propagators that leave U + F unchanged (Pak's criterion), found as
    automorphisms of a coloured graph (variables and monomials) with Sage's graph code.
    Each one is a symmetry of every integral of the family.  Returns generators for Reducer.
    """
    U, F = fam.UF()
    G = U + F
    R = G.parent()
    t = fam.t
    gr = Graph(multiedges=False, loops=False)
    xs = [('x', i) for i in range(t)]
    gr.add_vertices(xs)
    colours = {}
    for j, (c, m) in enumerate(zip(G.coefficients(), G.monomials())):
        v = ('m', j)
        gr.add_vertex(v)
        colours.setdefault(str(c), []).append(v)
        exps = m.exponents()[0]
        for i, e in enumerate(exps if hasattr(exps, '__iter__') else (exps,)):   # one line: an int
            if e:
                gr.add_edge(v, ('x', i), e)
    partition = [xs] + list(colours.values())
    A = gr.automorphism_group(partition=partition, edge_labels=True)
    gens = []
    for g in A.gens():
        perm = tuple(g(('x', i))[1] for i in range(t))
        if perm != tuple(range(t)):
            gens.append(perm)
    return gens


# ---------------------------------------------------------------------------- reduction
def _tidy_coefficient(c):
    """A rational function with its common factors cancelled: numerator and denominator with
    integer coefficients and no common number (Sage's fraction fields of several variables keep
    factors like 2^22 in both)."""
    try:
        n, d = c.numerator(), c.denominator()
        g = n.gcd(d)
        n, d = n // g, d // g
        L = lcm([QQ(x).denominator() for x in list(n.coefficients()) + list(d.coefficients())])
        n, d = n * L, d * L
        G = gcd([ZZ(x) for x in list(n.coefficients()) + list(d.coefficients())])
        if d.lc() < 0:
            G = -G
        return c.parent()(n / G) / c.parent()(d / G)
    except (AttributeError, TypeError, ValueError, ArithmeticError):
        return c


class Reduction:
    """What ibp_reduce() returns: {target: {master: coefficient}} with printing and .info()."""

    def __init__(self, table, masters, method, seconds, fam):
        table = {t: {m: _tidy_coefficient(c) for m, c in row.items()} for t, row in table.items()}
        self.table, self.masters, self.method, self.seconds, self.fam = table, masters, method, seconds, fam

    def __getitem__(self, target):
        """r["J(2,1)"] or r[(2, 1)]."""
        return self.table[_target(target)]

    def _lab(self, a):
        """An integral as the family's name with its powers, e.g. T(2) or J(1,1)."""
        return "%s(%s)" % (getattr(self.fam, 'name', 'F'), ",".join(str(x) for x in a))

    def draw(self, graph, targets=None, rename=None, size=1.5, fontsize=12):
        """The reductions as equations of diagrams (see plotting.draw_reduction): `graph` is a
        FeynmanGraph whose lines are this family's propagators in the same order.  targets are
        given as "J(2,1)" or (2, 1)."""
        from .plotting import draw_reduction
        if targets is not None:
            targets = [_target(t) for t in targets]
        return draw_reduction(self.table, graph, targets=targets, rename=rename, size=size, fontsize=fontsize)

    def info(self):
        print("Reduction of %d integrals of family %r with the %s reducer (%.1f s)."
              % (len(self.table), self.fam.name, self.method, self.seconds))
        print("Each target %s is written as sum_M c_M(d, invariants) M over the master integrals"
              % self._lab(["a1", "a2", "..."][:max(1, min(3, getattr(self.fam, 't', 3)))]))
        print("M = %s; d is the space-time dimension." % ", ".join(self._lab(m) for m in self.masters))

    def _repr_latex_(self):
        from sage.all import latex
        rows = []
        for t, row in self.table.items():
            terms = []
            body = ''
            for k, (mm, c) in enumerate(row.items()):
                term = r'%s\; %s' % (_latex_ratio(c), self._lab(mm))
                if k == 0:
                    body = term
                elif term.startswith('-'):
                    body += r' \\ &\quad - ' + term[1:]
                else:
                    body += r' \\ &\quad + ' + term
            rows.append(r'%s &= %s' % (self._lab(t), body or '0'))
        return r'$$\begin{aligned} %s \end{aligned}$$' % r' \\[6pt] '.join(rows)

    def _latex_(self):
        # used by Sage's latex() and by %display latex in notebooks
        return self._repr_latex_()[2:-2]

    def __repr__(self):
        lines = []
        for t, row in self.table.items():
            terms = ["[%s] %s" % (c.factor() if hasattr(c, 'factor') else c, self._lab(m)) for m, c in row.items()]
            lines.append("%s = %s" % (self._lab(t), " + ".join(terms) if terms else "0"))
        return "\n".join(lines)


def _latex_ratio(c):
    """A rational function as  k N / D  in LaTeX with N and D factored and one rational constant k
    (factoring numerator and denominator apart leaves a constant in each, e.g. (-16) ... / (-16) ...)."""
    from sage.all import latex, QQ
    from sage.structure.factorization import Factorization
    if not (hasattr(c, 'numerator') and hasattr(c, 'denominator')):
        return latex(c.factor() if hasattr(c, 'factor') else c)
    try:
        fN, fD = c.numerator().factor(), c.denominator().factor()
        k = QQ(fN.unit()) / QQ(fD.unit())
    except (TypeError, ValueError, ArithmeticError, AttributeError):
        return latex(c)
    N, Dn = list(fN), list(fD)
    sign = '-' if k < 0 else ''
    k = abs(k)
    num = ' '.join(([latex(k.numerator())] if k.numerator() != 1 or not N else []) +
                   [latex(Factorization([(p, e)])) for p, e in N])
    den = ' '.join(([latex(k.denominator())] if k.denominator() != 1 else []) +
                   [latex(Factorization([(p, e)])) for p, e in Dn])
    return sign + (r'\frac{%s}{%s}' % (num, den) if den else num)


def _target(t):
    """'T(2,1)' or 'F(2,1)' (any family name) or (2, 1) -> (2, 1)."""
    if isinstance(t, str):
        inside = re.search(r'\(([^)]*)\)', t)
        return tuple(int(x) for x in re.findall(r'-?\d+', inside.group(1) if inside else t))
    return tuple(t)


def ibp_reduce(fam, targets, method="auto", symmetries_="auto", nproc=None, explain=False, verbose=False):
    r"""
    Reduce integrals of a family to master integrals.

        ibp_reduce(kite, ["F(1,1,1,1,1)", (2,2,1,2,2)])

    The seed range is chosen from the targets (enough dots and numerator powers), the
    symmetries are found automatically (symmetries(fam)), and the method is
    "ff" (finite fields, fast; at most one symbolic invariant besides d), "exact"
    (Laporta with exact rational functions on every seeded equation) or "trimmed" (exact
    Laporta on only the equations the targets need, found by one modular probe: exact
    arithmetic, much faster than "exact").  "auto" picks ff when it can.
    The finite-field samples run on every core when the system is big enough (nproc=None);
    nproc=1 keeps it serial.
    """
    import time
    if explain:
        print(reduce.__doc__)
    ts = [_target(t) for t in targets]
    rmax = max(sum(x - 1 for x in t if x > 0) for t in ts)
    smax = max(sum(-x for x in t if x < 0) for t in ts)
    if len(fam.loops) > 1:
        smax = max(smax, 1)
    syms = symmetries(fam) if symmetries_ == "auto" else (symmetries_ or [])
    red = Reducer(fam, symmetries=syms)
    red.set_top(ts)                                    # seeds and masters inside the targets' sectors
    if method == "auto":
        method = "ff" if fam.kin.R.ngens() <= 2 else "exact"
    t0 = time.time()
    if method == "ff":
        from .ff import reduce_ff
        table = reduce_ff(red, ts, rmax=rmax, smax=smax, verbose=verbose, nproc=nproc)
    elif method == "trimmed":
        from .ff import reduce_exact_trimmed
        table = reduce_exact_trimmed(red, ts, rmax=rmax, smax=smax, verbose=verbose)
    else:
        red.run(rmax=rmax, smax=smax)
        table = {t: red.reduce(t) for t in ts}
    masters = sorted({m for row in table.values() for m in row}, key=weight)
    return Reduction(table, masters, method, time.time() - t0, fam)
