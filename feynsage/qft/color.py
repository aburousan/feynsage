r"""
SU(N) colour algebra (QCD colour factors), reduced by FORM.

    a, b, c = color_indices("a b c")         # adjoint (gluon) indices, N^2 - 1 values
    i, j, k = quark_colors("i j k")          # fundamental (quark) indices, N values
    T_color(a, i, j)                         # (T^a)_{ij}, normalised Tr(T^a T^b) = delta^{ab}/2
    f_color(a, b, c)                         # f^{abc}, [T^a, T^b] = i f^{abc} T^c
    color_trace(a, b, c)                     # Tr(T^a T^b T^c) with fresh internal quark indices
    color_chain(a, b, i=i, j=j)              # (T^a T^b)_{ij}
    color_factor(expr, N=3)                  # sum every repeated colour index

Every repeated colour index is summed.  The reduction is the Fierz identity
    (T^a)_{ij} (T^a)_{kl} = 1/2 (delta_{il} delta_{kj} - delta_{ij} delta_{kl} / N)
with f^{abc} = -2 i Tr(T^a [T^b, T^c]) and Tr(T^a T^b) = delta^{ab}/2 for what is left, so the answer is
a polynomial in N and 1/N: C_F = (N^2 - 1)/(2 N), C_A = N.  Open indices come back as Kron(i, j)
(delta_{ij}) and Kron(a, b) (delta^{ab}).
"""
import itertools
from sage.all import SR, I, function, latex
from . import tensors as T

_COUNT = itertools.count(1)
_lx = latex

COLT = function('T_color', nargs=3, print_latex_func=lambda self, a, i, j: r'(T^{%s})_{%s %s}' % (_lx(a), _lx(i), _lx(j)))
COLF = function('f_color', nargs=3, print_latex_func=lambda self, a, b, c: r'f^{%s %s %s}' % (_lx(a), _lx(b), _lx(c)))
KRON = function('Kron', nargs=2, print_latex_func=lambda self, a, b: r'\delta_{%s %s}' % (_lx(a), _lx(b)))
NC = SR.var('N')


def color_indices(names):
    """Adjoint (gluon) colour indices a, b, ... = 1 .. N^2 - 1."""
    return T._declare(names, 'adjoint')


def quark_colors(names):
    """Fundamental (quark) colour indices i, j, ... = 1 .. N."""
    return T._declare(names, 'fundamental')


def is_adjoint(x):
    try:
        return x.is_symbol() and T._KIND.get(str(x)) == 'adjoint'
    except AttributeError:
        return False


def is_fundamental(x):
    try:
        return x.is_symbol() and T._KIND.get(str(x)) == 'fundamental'
    except AttributeError:
        return False


def _need(x, test, what, fn):
    x = SR(x)
    if not test(x):
        raise TypeError("%s: %s is not %s (declare it with %s)" % (fn, x, what,
                        'color_indices' if test is is_adjoint else 'quark_colors'))
    return x


def T_color(a, i, j):
    """(T^a)_{ij}: the SU(N) generator in the fundamental representation."""
    return COLT(_need(a, is_adjoint, "a gluon colour index", "T_color"),
                _need(i, is_fundamental, "a quark colour index", "T_color"),
                _need(j, is_fundamental, "a quark colour index", "T_color"))


def color_delta(x, y):
    """delta_{ij} for two quark indices or delta^{ab} for two gluon indices."""
    x, y = SR(x), SR(y)
    if not ((is_fundamental(x) and is_fundamental(y)) or (is_adjoint(x) and is_adjoint(y))):
        raise TypeError("color_delta: give two quark colour indices or two gluon colour indices")
    a, b = sorted((x, y), key=str)
    return KRON(a, b)


def f_color(a, b, c):
    """The structure constant f^{abc} (totally antisymmetric)."""
    args = [_need(x, is_adjoint, "a gluon colour index", "f_color") for x in (a, b, c)]
    return T._eps3_sorted(args)


def _fresh_quark():
    while True:
        n = 'qc%d' % next(_COUNT)
        if n not in T._KIND:
            return quark_colors(n)


def color_chain(*adj, i=None, j=None):
    """(T^{a1} T^{a2} ... )_{ij} with internal quark indices summed."""
    ks = [SR(i) if i is not None else _fresh_quark()] + [_fresh_quark() for _ in adj[:-1]] + \
         [SR(j) if j is not None else _fresh_quark()]
    out = SR(1)
    for n, a in enumerate(adj):
        out *= T_color(a, ks[n], ks[n + 1])
    return out


def color_trace(*adj):
    """Tr(T^{a1} ... T^{an}) (Tr 1 = N for no argument)."""
    if not adj:
        return NC
    k0 = _fresh_quark()
    return color_chain(*adj, i=k0, j=k0)


def color_factor(expr, N=None, rules=None, debug=False):
    """Sum all repeated colour indices of expr (a Sage expression) with FORM; N = 3 for QCD, or a
    symbol (the default N)."""
    from .operations import contract
    r = contract(expr, rules=rules, debug=debug)
    return r.subs({NC: N}) if N is not None else r


# ---------------------------------------------------------------------------- for the compiler
def expand_f(e):
    """f^{abc} -> -2 i [Tr(T^a T^b T^c) - Tr(T^a T^c T^b)] with fresh quark indices for every f."""
    e = SR(e)
    w = [SR.wild(n) for n in range(3)]
    fs = e.find(COLF(*w))
    if not fs:
        return e
    # each occurrence needs its own internal indices: expand monomial by monomial
    from .compiler import _opname
    op = _opname(e)
    if op in ('add_vararg', 'add'):
        return sum((expand_f(x) for x in e.operands()), SR(0))
    if op in ('mul_vararg', 'mul'):
        out = SR(1)
        for x in e.operands():
            out *= expand_f(x)
        return out
    if op == 'pow':
        b, n = e.operands()
        if n.is_integer() and n > 0:
            out = SR(1)
            for _ in range(int(n)):
                out *= expand_f(b)
            return out
        raise TypeError("f_color in a denominator")
    if op == 'f_color':
        a, b, c = e.operands()
        return -2 * I * (color_trace(a, b, c) - color_trace(a, c, b))
    raise TypeError("f_color inside the function %s" % op)


FORM_RULES = """repeat;
  id FSCT(fsi1?,fsi2?,fsa1?)*FSCT(fsi3?,fsi4?,fsa1?) = 1/2*(d_(fsi1,fsi4)*d_(fsi3,fsi2) - FSNC^-1*d_(fsi1,fsi2)*d_(fsi3,fsi4));
  id FSCT(fsi1?,fsi1?,fsa1?) = 0;
  id FSCT(fsi1?,fsi2?,fsa1?)*FSCT(fsi2?,fsi1?,fsa2?) = 1/2*d_(fsa1,fsa2);
endrepeat;
id FSNA = FSNC^2-1;"""
