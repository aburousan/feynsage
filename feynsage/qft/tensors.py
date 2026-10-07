r"""
Momenta, Lorentz indices and the commuting Lorentz tensors, as ordinary Sage symbolic expressions.

    p, q = momenta("p q")
    mu, nu = lorentz_indices("mu nu")
    dot(p, q)            p.q          (dot(p, p) is p^2)
    comp(p, mu)          p^mu
    metric(mu, nu)       g^{mu nu}    (delta_{mu nu} with euclidean=True)
    epsilon(mu, nu, p, q)             the Levi-Civita tensor

Momenta and indices are plain Sage symbols, so p + q, 2*p, x*k work as usual.  The tensors are Sage
symbolic functions with their arguments in a fixed order (dot and metric are symmetric, epsilon is
antisymmetric) and linear in each momentum, so dot(p + k, q) = dot(p, q) + dot(k, q) and results can
be substituted with .subs({dot(p, p): m^2}).

Conventions (see set_convention):
    "standard"  Peskin and Schroeder: metric (+,-,-,-), gamma5 = i g^0 g^1 g^2 g^3,
                tr(g^mu g^nu g^rho g^sigma g5) = -4 i eps^{mu nu rho sigma} (their eq. 5.5),
                so eps^{0123} = +1; the symbolic dimension is called d.
    "note"      the lecture note: the same Minkowski algebra, the dimension is called D = 4 - 2 eps.
    euclidean=True   {g_mu, g_nu} = 2 delta_{mu nu}, gamma5 = g_1 g_2 g_3 g_4,
                tr(g_mu g_nu g_rho g_sigma g5) = 4 eps_{mu nu rho sigma}, eps_{1234} = +1,
                p.q the Euclidean product.
"""
from sage.all import SR, function, latex, Integer

# ---------------------------------------------------------------------------- conventions
_CONV = {'name': 'standard', 'euclidean': False}


def set_convention(name=None, euclidean=None):
    """Choose "standard" (Peskin and Schroeder, dimension d) or "note" (the lecture note,
    dimension D = 4 - 2 eps), and Minkowski (default) or euclidean=True.  Returns the setting."""
    if name is not None:
        if name not in ('standard', 'note'):
            raise ValueError('convention must be "standard" or "note"')
        _CONV['name'] = name
    if euclidean is not None:
        _CONV['euclidean'] = bool(euclidean)
    return dict(_CONV)


def convention():
    """The current convention, e.g. {'name': 'standard', 'euclidean': False}."""
    return dict(_CONV)


def dimension():
    """The symbol for the space-time dimension: d ("standard") or D ("note")."""
    return SR.var('D') if _CONV['name'] == 'note' else SR.var('d')


def _euclid(euclidean):
    return _CONV['euclidean'] if euclidean is None else bool(euclidean)


# ---------------------------------------------------------------------------- registry
_KIND = {}          # symbol name -> 'vector' or 'index'


def _names(names):
    if isinstance(names, str):
        return names.replace(',', ' ').split()
    return [str(n) for n in names]


def _declare(names, kind, latex_name=None):
    out = []
    names = _names(names)
    if latex_name is not None:
        latex_name = [latex_name] if isinstance(latex_name, str) else list(latex_name)
        if len(latex_name) != len(names):
            raise ValueError("give one latex_name for each name")
    for i, n in enumerate(names):
        if not n.isidentifier():
            raise ValueError("%r is not a valid name (letters, digits and _ only)" % n)
        old = _KIND.get(n)
        if old is not None and old != kind:
            raise ValueError("%s is already a %s; choose another name" % (n, 'momentum' if old == 'vector' else 'Lorentz index'))
        _KIND[n] = kind
        out.append(SR.var(n, latex_name=latex_name[i]) if latex_name else SR.var(n))
    return out[0] if len(out) == 1 else tuple(out)


def momenta(names, latex_name=None):
    """momenta("p q k") -> the symbols p, q, k, marked as four-vectors.  latex_name gives how they
    print in LaTeX, for example momenta("pp kp", latex_name=["p'", "k'"])."""
    return _declare(names, 'vector', latex_name)


def lorentz_indices(names, latex_name=None):
    """lorentz_indices("mu nu") -> the symbols mu, nu, marked as Lorentz indices (latex_name as in momenta)."""
    return _declare(names, 'index', latex_name)


def is_momentum(x):
    try:
        return x.is_symbol() and _KIND.get(str(x)) == 'vector'
    except AttributeError:
        return False


def is_index(x):
    try:
        return x.is_symbol() and _KIND.get(str(x)) == 'index'
    except AttributeError:
        return False


def _fresh_index(base, taken):
    """A new index of the same kind as base, named base_c, base_cc, ... not in taken."""
    kind = _KIND.get(str(base), 'index')
    n = str(base) + '_c'
    while n in taken or (n in _KIND and _KIND[n] != kind):
        n += 'c'
    return _declare(n, kind)


# ---------------------------------------------------------------------------- the symbolic functions
def _lx(a):
    return latex(a)


def _latex_dot(self, a, b):
    E = r'_{E}' if _CONV['euclidean'] else ''           # the note writes Euclidean squares p_E^2
    if str(a) == str(b):
        return r'{{%s}%s^{2}}' % (_lx(a), E)                       # a group, so that (p^2)^2 is valid LaTeX
    return r'\left({%s}%s \cdot {%s}%s\right)' % (_lx(a), E, _lx(b), E)    # (p.k)(p'.k') as in the books


# Names differ from the plain function('dot'), function('g'), ... of form.py on purpose: two Sage
# functions with one name and different options disturb each other in Pynac's registry.
DOT = function('Dot', nargs=2, print_latex_func=_latex_dot)
COMP = function('Comp', nargs=2, print_latex_func=lambda self, p, mu: r'{%s}^{%s}' % (_lx(p), _lx(mu)))
G = function('Metric', nargs=2, print_latex_func=lambda self, a, b: r'g^{%s %s}' % (_lx(a), _lx(b)))
DELTA = function('Delta', nargs=2, print_latex_func=lambda self, a, b: r'\delta_{%s %s}' % (_lx(a), _lx(b)))
EPS = function('Epsilon', nargs=4, print_latex_func=lambda self, *a: r'\epsilon^{%s}' % ' '.join(_lx(x) for x in a))

# Sage function name -> what it is; the lower-case names are those of form.compute ('eps' is FORM's e_)
TENSOR_NAMES = {'Dot': 'dot', 'dot': 'dot', 'Comp': 'comp', 'comp': 'comp', 'Metric': 'g', 'g': 'g',
                'Delta': 'delta', 'delta': 'delta', 'Epsilon': 'Eps', 'Eps': 'Eps', 'eps': 'eps'}


# ---------------------------------------------------------------------------- linear algebra of momenta
def vector_parts(x, what="a momentum"):
    """x = sum c_i p_i (registered momenta p_i) -> [(p_i, c_i)]; raises if anything else is left."""
    x = SR(x)
    if is_momentum(x):
        return [(x, Integer(1))]
    if x.is_trivial_zero():
        return []
    vs = [v for v in x.variables() if is_momentum(v)]
    if not vs:
        raise TypeError("expected %s (a combination of momenta from momenta(...)), got %s" % (what, x))
    e = x.expand()
    parts, rest = [], e
    for v in sorted(vs, key=str):
        c = e.coefficient(v, 1)
        if any(c.has(w) for w in vs):
            raise TypeError("%s is not linear in the momenta" % x)
        parts.append((v, c))
        rest = rest - c * v
    rest = rest.expand()
    if not rest.is_trivial_zero():
        raise TypeError("%s: the part %s is not a momentum.  Write slash(p) + m, not slash(p + m)" % (x, rest))
    return parts


def _index_or_vector(x):
    x = SR(x)
    if is_index(x):
        return [(x, Integer(1))]
    if x.is_symbol() and not is_momentum(x):
        raise TypeError("%s is neither a Lorentz index (lorentz_indices) nor a momentum (momenta)" % x)
    return vector_parts(x)


def dot(a, b):
    """The scalar product a.b of two momenta (sums of momenta are expanded); dot(p, p) = p^2."""
    out = SR(0)
    for p, cp in vector_parts(a):
        for q, cq in vector_parts(b):
            x, y = sorted((p, q), key=str)
            out += cp * cq * DOT(x, y)
    return out


def comp(p, mu):
    """The component p^mu of a momentum (or of a sum of momenta)."""
    mu = SR(mu)
    if not is_index(mu):
        raise TypeError("comp(p, mu): %s is not a Lorentz index" % mu)
    return sum((c * COMP(v, mu) for v, c in vector_parts(p)), SR(0))


def metric(mu, nu, euclidean=None):
    """The metric g^{mu nu} (delta_{mu nu} in the Euclidean convention)."""
    mu, nu = SR(mu), SR(nu)
    for x in (mu, nu):
        if not is_index(x):
            raise TypeError("metric(mu, nu): %s is not a Lorentz index" % x)
    a, b = sorted((mu, nu), key=str)
    return (DELTA if _euclid(euclidean) else G)(a, b)


def _perm_sign(seq):
    seq, sign = list(seq), 1
    for i in range(len(seq)):
        for j in range(len(seq) - 1 - i):
            if seq[j] > seq[j + 1]:
                seq[j], seq[j + 1] = seq[j + 1], seq[j]
                sign = -sign
    return sign


def _eps_sorted(args, fn=None):
    names = [str(a) for a in args]
    if len(set(names)) < 4:
        return SR(0)
    order = sorted(range(4), key=lambda i: names[i])
    return _perm_sign(names) * (fn or EPS)(*[args[i] for i in order])


def _eps3_sorted(args):
    from .color import COLF
    names = [str(a) for a in args]
    if len(set(names)) < 3:
        return SR(0)
    order = sorted(range(3), key=lambda i: names[i])
    return _perm_sign(names) * COLF(*[args[i] for i in order])


def is_summable(x):
    """A Lorentz index or a colour index (summed when it appears twice)."""
    try:
        return x.is_symbol() and _KIND.get(str(x)) in ('index', 'adjoint', 'fundamental')
    except AttributeError:
        return False


def epsilon(a, b, c, d):
    """The Levi-Civita tensor with indices or momenta in its slots (eps^{0123} = +1 in Minkowski
    space, eps_{1234} = +1 in Euclidean space); it is linear in momentum slots and antisymmetric."""
    terms = [([], SR(1))]
    for x in (a, b, c, d):
        terms = [(args + [v], co * cv) for args, co in terms for v, cv in _index_or_vector(x)]
    return sum((co * _eps_sorted(args) for args, co in terms), SR(0))


# The printed names of results, usable as input too: 4*Dot(p, q)*Metric(mu, nu) can be pasted back.
Dot, Comp, Metric, Epsilon = dot, comp, metric, epsilon


def Delta(mu, nu):
    """The Euclidean metric delta_{mu nu}."""
    return metric(mu, nu, euclidean=True)


# ---------------------------------------------------------------------------- polarization vectors
_POL = {}           # name -> (momentum, name of the complex conjugate vector, mass)


def polarization(name, k, mass=0, latex_name=None):
    r"""A polarization vector eps(k) of a photon, gluon or massive vector boson with momentum k:
        e1 = polarization("e1", k1)
    It is a momentum-like vector (use slash(e1), dot(e1, p), comp(e1, mu)); its complex conjugate
    eps*(k) is made automatically (named e1_c, returned by conjugate_vector(e1)) and conjugate()
    exchanges the two.  Sum over the polarizations with polarization_sum.  latex_name, for example
    r"\epsilon'", is how it prints; the conjugate prints with a star."""
    k = SR(k)
    if not is_momentum(k):
        raise TypeError("polarization(%r, k): k must be a momentum from momenta()" % name)
    e = momenta(name, latex_name=latex_name)
    ec = momenta(str(name) + '_c', latex_name=(r'{%s}^{*}' % latex_name) if latex_name else None)
    _POL[str(e)] = (k, str(ec), SR(mass))
    _POL[str(ec)] = (k, str(e), SR(mass))
    return e


def conjugate_vector(x):
    """The complex conjugate partner of a polarization vector (the vector itself otherwise)."""
    hit = _POL.get(str(x))
    return SR.var(hit[1]) if hit else x


# ---------------------------------------------------------------------------- Minkowski <-> Euclidean
def _flip_dots(expr):
    expr = SR(expr)
    w0, w1 = SR.wild(0), SR.wild(1)
    if any(T_ in str(expr) for T_ in ('Comp(', 'Metric(', 'Delta(', 'Epsilon(')):
        raise TypeError("to_euclidean/to_minkowski work on scalar expressions (scalar products only); "
                        "contract free indices and remove epsilon first")
    dots = expr.find(DOT(w0, w1))
    return expr.subs({x: -x for x in dots}) if dots else expr


def to_euclidean(expr):
    """Minkowski scalar products -> Euclidean ones, p.q = -p_E.q_E (l^0 = i l_4; so s = p^2 = -p_E^2 as
    in the note).  For a scalar expression; the answer is in Euclidean Dot's."""
    return _flip_dots(expr)


def to_minkowski(expr):
    """The inverse of to_euclidean: p_E.q_E = -p.q."""
    return _flip_dots(expr)


def normalize_dim(dim):
    """4, a symbol, or a name such as "d" or "D" -> 4 or the Sage symbol."""
    if dim is None:
        return 4
    if isinstance(dim, str):
        return 4 if dim.strip() == '4' else SR.var(dim.strip())
    try:
        if SR(dim) == 4:
            return 4
    except TypeError:
        pass
    return SR(dim)
