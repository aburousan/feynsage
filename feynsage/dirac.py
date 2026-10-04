r"""
Dirac algebra on open fermion lines, in the style of Package-X's FermionLine and FermionLineExpand.

    from feynsage.dirac import line_expand
    line_expand(("u", "p2", "m"), ["mu", "p2 + m", "nu"], ("v", "p1", "m"), on_shell=True)
        # u-bar(p2) g^mu (p2/ + m) g^nu v(p1)  ->  {structure: coefficient}

A product of Dirac matrices is a list of factors, as in form.dirac_trace:
    'mu'              gamma^mu (any name that is not a momentum)
    'p - k + m'       slashed momenta plus a multiple of the unit matrix
    'g5', 'PL', 'PR'  gamma_5 and the chiral projectors (1 -+ g5)/2
The algebra is d-dimensional with an anticommuting gamma_5, as in Package-X.  gamma_5, PL and PR
are moved to the right end; the gamma matrices are put in a fixed order with {g^a, g^b} = 2 g^ab, so
a repeated index contracts to d and a repeated momentum to its square.  On a line the momentum of the
left spinor is ordered to the left end and that of the right spinor to the right end, where the
Dirac equation removes it:
    p/ u(p) = m u(p),  p/ v(p) = -m v(p),  u-bar(p) p/ = m u-bar(p),  v-bar(p) p/ = -m v-bar(p).
The result is written in the basis 1, g^a, sigma^{ab} = (i/2)[g^a, g^b] and ordered products of three
or more distinct gamma matrices, each times 1, g5 or PL, PR.
"""
import re
from sage.all import SR, I
from ._parse import sr

D = SR.var('d')


# ---------------------------------------------------------------------------- symbols
def dot(a, b):
    """The symbol for a.b (a^2 for a = b), as in loop(): p1.p2 -> p1p2, p.p -> p2 is avoided: 'p^2' -> p_sq."""
    a, b = sorted((a, b))
    return SR.var('%s_sq' % a) if a == b else SR.var('%s_dot_%s' % (a, b))


def comp(p, mu):
    """The component p^mu."""
    return SR.var('%s__%s' % (p, mu))


def metric(mu, nu):
    """g^{mu nu}."""
    mu, nu = sorted((mu, nu))
    return D if mu == nu else SR.var('g__%s__%s' % (mu, nu))


# ---------------------------------------------------------------------------- parsing
def _items(factors, momenta):
    """Factor strings -> terms [(coefficient, [items])], items ('i', mu), ('p', p) or ('c', '5'|'L'|'R')."""
    terms = [(SR(1), [])]
    for f in factors:
        f = str(f).strip()
        if f in ('g5', 'PL', 'PR'):
            terms = [(c, s + [('c', {'g5': '5', 'PL': 'L', 'PR': 'R'}[f])]) for c, s in terms]
            continue
        if re.fullmatch(r'[A-Za-z_]\w*', f) and f not in momenta:
            terms = [(c, s + [('i', f)]) for c, s in terms]
            continue
        e = sr(f).expand()
        parts, rest = [], e
        for p in momenta:
            cp = e.coefficient(SR.var(p))
            if not cp.is_zero():
                parts.append((cp, ('p', p)))
                rest -= cp * SR.var(p)
        rest = rest.expand()
        new = []
        for c, s in terms:
            for cp, it in parts:
                new.append((c * cp, s + [it]))
            if not rest.is_zero():
                new.append((c * rest, s))
        terms = new
    return terms


def _chiral_mul(a, b):
    """(sign, product) for two chiral matrices of {'1', '5', 'L', 'R'}; product '0' means zero."""
    if a == '1':
        return 1, b
    if b == '1':
        return 1, a
    if a == '5' and b == '5':
        return 1, '1'
    if a == '5':
        return (-1, 'L') if b == 'L' else (1, 'R')
    if b == '5':
        return (-1, 'L') if a == 'L' else (1, 'R')
    return (1, a) if a == b else (1, '0')


def _chiral_right(c, items):
    """Move the chiral matrices to the right end: C g^a = g^a C' with C' = -g5 for g5 and PL <-> PR."""
    slots, state = [], '1'
    for it in items:
        if it[0] == 'c':
            sg, state = _chiral_mul(state, it[1])
            c *= sg
            if state == '0':
                return SR(0), [], '1'
        else:
            slots.append(it)
            if state == '5':
                c = -c
            elif state in ('L', 'R'):
                state = 'R' if state == 'L' else 'L'
    return c, slots, state


# ---------------------------------------------------------------------------- ordering
class _Order:
    """Sort key: the left spinor's momentum first, then indices, other momenta, the right spinor's momentum."""

    def __init__(self, left=None, right=None, sp=None):
        self.left, self.right, self.sp = left, right, dict(sp or {})

    def key(self, it):
        kind, name = it
        if kind == 'p' and name == self.right:
            return (3, name)
        if kind == 'p' and name == self.left:
            return (0, name)
        return (1, name) if kind == 'i' else (2, name)

    def half_anti(self, a, b):
        """{g^a, g^b}/2."""
        if a[0] == 'i' and b[0] == 'i':
            return metric(a[1], b[1])
        if a[0] == 'p' and b[0] == 'p':
            k = tuple(sorted((a[1], b[1])))
            return self.sp.get(k, dot(*k))
        idx, mom = (a[1], b[1]) if a[0] == 'i' else (b[1], a[1])
        return comp(mom, idx)


def _contract_pair(a, b, rest, order):
    """{g^a, g^b}/2 with the rest of the product; an index that appears again in the rest is summed:
    g^{mu nu} turns the other gamma_mu into gamma^nu and p^mu turns it into p/."""
    idx = [x[1] for x in (a, b) if x[0] == 'i']
    for mu in idx:
        if ('i', mu) in rest:
            j = rest.index(('i', mu))
            other = b if a == ('i', mu) else a
            return 1, rest[:j] + [other] + rest[j + 1:]
    return order.half_anti(a, b), rest


def _normal_order(c, slots, order):
    """{tuple(slots): coefficient} with the slots sorted and no two equal."""
    out = {}
    todo = [(c, list(slots))]
    while todo:
        c, s = todo.pop()
        if c.is_zero():
            continue
        for i in range(len(s) - 1):
            a, b = s[i], s[i + 1]
            if a == b:
                todo.append((c * order.half_anti(a, b), s[:i] + s[i + 2:]))
                break
            if order.key(a) > order.key(b):
                todo.append((-c, s[:i] + [b, a] + s[i + 2:]))
                f, rest = _contract_pair(a, b, s[:i] + s[i + 2:], order)
                todo.append((2 * c * f, rest))
                break
        else:
            k = tuple(s)
            out[k] = out.get(k, SR(0)) + c
    return out


# ---------------------------------------------------------------------------- the line
def _sign(kind):
    if kind not in ('u', 'v'):
        raise ValueError("a spinor is 'u' or 'v'")
    return 1 if kind == 'u' else -1


def _collect(dct):
    return {k: v for k, v in ((k, v.expand()) for k, v in dct.items()) if not v.is_zero()}


def line_expand(left, factors, right, sp=None, momenta=None, gordon=True, on_shell=True):
    r"""
    Expand w-bar(p2, m2) Gamma w(p1, m1) (Package-X's FermionLineExpand):
        left  = ('u' or 'v', 'p2', m2)    the left spinor u-bar(p2) or v-bar(p2)
        right = ('u' or 'v', 'p1', m1)    the right spinor u(p1) or v(p1)
        factors   the Dirac matrices between them, as in form.dirac_trace
        sp        scalar products {('p1', 'p2'): value} (on shell p_i^2 = m_i^2 is added)
        gordon    rewrite (p1 + p2)^mu w-bar w and (p1 + p2)^mu w-bar g5 w with the Gordon identities
    Returns {(structure, chiral): coefficient}; structure is () for 1, (('i', mu),) for g^mu,
    ('sigma', a, b) for sigma^{ab} (a, b an index or a momentum), or an ordered tuple of slots, and
    chiral is '1', '5', 'L' or 'R'.  Symbols in the coefficients: d, p_sq, p1_dot_p2, p__mu (p^mu),
    g__mu__nu (g^{mu nu}).
    """
    lk, p2, m2 = left
    rk, p1, m1 = right
    s2, s1 = _sign(lk), _sign(rk)
    m1, m2 = sr(m1), sr(m2)
    moms = list(momenta or [])
    for p in (p1, p2):
        if p not in moms:
            moms.append(p)
    spd = dict(sp or {})
    if on_shell:
        spd[(p1, p1)] = m1 ** 2
        spd[(p2, p2)] = m2 ** 2
    out = _core(_items(factors, moms), p1, p2, s1 * m1, s2 * m2, spd)
    if gordon:
        out = _gordon(out, p1, p2, s1 * m1, s2 * m2)
    return _collect(out)


def _core(terms, p1, p2, M1, M2, spd):
    """Parsed terms -> {(structure, chiral): coefficient}, Dirac equation with M = s m at both ends."""
    order = _Order(left=p2, right=p1, sp={tuple(sorted(k)): sr(v) for k, v in spd.items()})
    out = {}
    for c, items in terms:
        c, slots, ch = _chiral_right(c, items)
        if c.is_zero():
            continue
        for s, cc in _normal_order(c, slots, order).items():
            s, chs = list(s), ch
            if s and s[-1] == ('p', p1):           # Dirac equation on the right, through the chiral end
                s.pop()
                cc *= M1
                if chs == '5':
                    cc = -cc
                elif chs in ('L', 'R'):
                    chs = 'R' if chs == 'L' else 'L'
            if s and s[0] == ('p', p2):            # and on the left
                s.pop(0)
                cc *= M2
            k = (tuple(s), chs)
            out[k] = out.get(k, SR(0)) + cc
    return _sigma_basis(_collect(out))


def _sigma_basis(dct):
    """Ordered pairs g^a g^b = g^ab - i sigma^{ab} (a before b in the canonical order)."""
    out = {}
    order = _Order()
    for (s, ch), c in dct.items():
        if len(s) == 2:
            a, b = s
            k0 = ((), ch)
            out[k0] = out.get(k0, SR(0)) + c * order.half_anti(a, b)
            k1 = (('sigma', a, b), ch)
            out[k1] = out.get(k1, SR(0)) - I * c
        else:
            k = (s, ch)
            out[k] = out.get(k, SR(0)) + c
    return out


def _gordon(dct, p1, p2, M1, M2):
    r"""
    (p1 + p2)^mu w-bar w     = (M1 + M2) w-bar g^mu w - i w-bar sigma^{mu nu} (p2 - p1)_nu w
    (p1 + p2)^mu w-bar g5 w  = (M2 - M1) w-bar g^mu g5 w - i w-bar sigma^{mu nu} (p2 - p1)_nu g5 w
    with M = s m (s = +1 for u, -1 for v).  Applied to the coefficient of p1^mu + p2^mu in front of
    the scalar and pseudoscalar structures; the difference (p2 - p1)^mu stays.
    """
    out = dict(dct)
    for ch in ('1', '5', 'L', 'R'):
        k = ((), ch)
        if k not in out:
            continue
        c = out[k]
        idx = {str(v).split('__')[1] for v in c.variables() if str(v).startswith(p1 + '__') or str(v).startswith(p2 + '__')}
        for mu in idx:
            a, b = comp(p1, mu), comp(p2, mu)
            ca, cb = c.coefficient(a), c.coefficient(b)
            sym = ((ca + cb) / 2).expand()                 # c = sym (p1 + p2)^mu + anti (p2 - p1)^mu + rest
            if sym.is_zero() or sym.has(a) or sym.has(b):
                continue
            c = (c - sym * (a + b)).expand()
            vec = ((('i', mu),),)
            if ch == '1':
                parts = [('1', M1 + M2)]
            elif ch == '5':
                parts = [('5', M2 - M1)]
            else:                                          # g^mu p1/ PL = g^mu PR p1/: the chirality flips
                other = 'R' if ch == 'L' else 'L'
                parts = [(ch, M2), (other, M1)]
            for chv, cv in parts:
                kv = vec + (chv,)
                out[kv] = out.get(kv, SR(0)) + sym * cv
            for mom, sgn in ((p2, 1), (p1, -1)):
                kt = (('sigma', ('i', mu), ('p', mom)), ch)
                out[kt] = out.get(kt, SR(0)) - I * sgn * sym
        out[k] = c
    return out


# ---------------------------------------------------------------------------- loop integrals on a line
def _kin_sp(kin):
    out = {}
    for key, val in (kin or {}).items():
        k = key.replace(' ', '')
        a, b = (k[:-2], k[:-2]) if k.endswith('^2') else k.split('.')
        out[tuple(sorted((a, b)))] = sr(val)
    return out


_LABEL = re.compile(r'(?:g|delta)\^\{(\w+) (\w+)\}|(\w+)\^(\w+)')


def _finalize_sr(dct, euclidean=False):
    """sum_M c_M(d) M with symbolic coefficients -> expression to eps^0 (as pv._finalize)."""
    from .pv import _master_formula, eps
    total = SR(0)
    for key, c in dct.items():
        full, (p2, p1) = _master_formula(key, euclidean)
        c = SR(c)
        c4 = c.subs({D: 4})
        dc4 = c.diff(D).subs({D: 4})
        d2c4 = c.diff(D, 2).subs({D: 4}) / 2
        total += c4 * full
        if not dc4.is_zero() and not (p2.is_trivial_zero() and p1.is_trivial_zero()):
            total += -2 * dc4 * (p2 / eps + p1)
        if not d2c4.is_zero() and not p2.is_trivial_zero():
            total += 2 * d2c4 * p2
    return total


class LineResult:
    """What loop_line() returns: {(structure, chiral): coefficient}, printed spinor by spinor."""

    def __init__(self, parts, left, right):
        self.parts, self.left, self.right = parts, left, right

    def coefficients(self):
        return dict(self.parts)

    def __getitem__(self, key):
        return self.parts[key]

    def __repr__(self):
        lk, p2, m2 = self.left
        rk, p1, m1 = self.right
        rows = []
        for (st, ch), c in sorted(self.parts.items(), key=lambda t: str(t[0])):
            if lk:
                rows.append("[%sbar(%s)| %s |%s(%s)] :  %s" % (lk, p2, structure_name(st, ch), rk, p1, c))
            else:
                rows.append("%s :  %s" % (structure_name(st, ch), c))
        return "\n".join(rows) if rows else "0"


def structure_name(st, ch):
    def one(x):
        return x[1] if x[0] == 'i' else x[1]
    if not st:
        body = "1"
    elif st[0] == 'sigma':
        body = "sigma^{%s %s}" % (one(st[1]), one(st[2]))
    elif st[0] == 'asym':
        body = "Gamma^{[%s]}" % " ".join(one(x) for x in st[1:])
    elif st[0] == 'eps':
        body = "eps^{%s %s %s s} g_s" % (one(st[1]), one(st[2]), one(st[3]))
    else:
        body = " ".join(("g^%s" % x[1]) if x[0] == 'i' else ("%s/" % x[1]) for x in st)
    return body + {'1': '', '5': ' g5', 'L': ' PL', 'R': ' PR'}[ch]


def loop_line(left, factors, right, *props, kin=None, gordon=True, on_shell=True):
    r"""
    A one-loop integral over a fermion line (Package-X's LoopIntegrate of a FermionLine, then LoopRefine):
        Int w-bar(p2, m2) Gamma(l) w(p1, m1) / prod D_i,     D_i = (l + q_i)^2 - m_i^2,
    in the Package-X normalisation mu^(2 eps) e^(eps gamma_E) Int d^d l/(i pi^(d/2)).

        loop_line(("u", "p2", "m"), ["rho", "p2 - l + m", "mu", "p1 - l + m", "rho"], ("u", "p1", "m"),
                  ["l", "0"], ["l - p2", "m"], ["l - p1", "m"], kin={"p1^2": "m^2", "p2^2": "m^2", "p1.p2": "m^2 - q2/2"})

    The factors are as in line_expand and may contain the loop momentum l.  Every l/ becomes
    gamma_a l^a; the tensor integral is reduced exactly in d, put back into the line, which is then
    simplified (contractions, Dirac equation), and only then expanded in eps, so the factors d of the
    contractions meet the poles correctly.  Returns a LineResult {(structure, chiral): coefficient}.
    """
    from .pv import loop, _SR
    lk, p2, m2 = left
    rk, p1, m1 = right
    s2, s1 = _sign(lk), _sign(rk)
    m1, m2 = sr(m1), sr(m2)
    spd = _kin_sp(kin)
    if on_shell:
        spd[(p1, p1)] = m1 ** 2
        spd[(p2, p2)] = m2 ** 2
    moms = sorted({n for q, _ in props for n in re.findall(r'[A-Za-z_]\w*', str(q))} | {p1, p2} |
                  {a for k in spd for a in k})
    if 'l' not in moms:
        moms.append('l')
    terms = _items(factors, moms)
    by_rank = {}
    for c, items in terms:
        new, n = [], 0
        for it in items:
            if it == ('p', 'l'):
                n += 1
                new.append(('i', 'fsL%d' % n))
            else:
                new.append(it)
        by_rank.setdefault(n, []).append((c, new))
    acc = {}
    for n, tl in by_rank.items():
        num = " ".join("l^fsL%d" % (j + 1) for j in range(n)) or "1"
        res = loop(num, *props, kin=kin)
        for label, dct in res.raw:
            ren, subs_p = {}, {}
            for mm in _LABEL.finditer(label):
                if mm.group(1):
                    ren[mm.group(2)] = mm.group(1)       # g^{a b}: b is summed with a
                else:
                    subs_p[mm.group(4)] = mm.group(3)    # v^a: gamma_a -> v/
            for c, items in tl:
                its = []
                for it in items:
                    if it[0] == 'i' and it[1] in subs_p:
                        its.append(('p', subs_p[it[1]]))
                    elif it[0] == 'i' and it[1] in ren:
                        its.append(('i', ren[it[1]]))
                    else:
                        its.append(it)
                for k, cc in _core([(c, its)], p1, p2, s1 * m1, s2 * m2, spd).items():
                    slot = acc.setdefault(k, {})
                    for mk, mc in dct.items():
                        slot[mk] = slot.get(mk, SR(0)) + cc * _SR(mc)
    out = {k: _finalize_sr(v) for k, v in acc.items()}
    if gordon:
        out = _gordon(out, p1, p2, s1 * m1, s2 * m2)
    return LineResult(_collect(out), left, right)


def loop_matrix(factors, *props, kin=None):
    r"""
    A one-loop integral of a product of Dirac matrices with no spinors (Package-X's LoopIntegrate of a
    DiracMatrix), for example a self-energy:
        loop_matrix(["mu", "p - l + m", "mu"], ["l", "0"], ["l - p", "m"], kin={"p^2": "s"})
    gives {(structure, chiral): coefficient} with the structures 1 and p/ (and sigma for more momenta).
    """
    return _loop_core(None, factors, None, props, kin, gordon=False, on_shell=False)


def _loop_core(left, factors, right, props, kin, gordon, on_shell):
    from .pv import loop, _SR
    spd = _kin_sp(kin)
    if left is not None:
        lk, p2, m2 = left
        rk, p1, m1 = right
        M2, M1 = _sign(lk) * sr(m2), _sign(rk) * sr(m1)
        if on_shell:
            spd[(p1, p1)] = sr(m1) ** 2
            spd[(p2, p2)] = sr(m2) ** 2
    else:
        p1 = p2 = None
        M1 = M2 = SR(0)
    moms = sorted({n for q, _ in props for n in re.findall(r'[A-Za-z_]\w*', str(q))} |
                  {x for x in (p1, p2) if x} | {a for k in spd for a in k} |
                  {n for f in factors for n in re.findall(r'[A-Za-z_]\w*', str(f))
                   if re.fullmatch(r'[A-Za-z_]\w*', str(f).strip()) is None})
    moms = [x for x in moms if x not in ('m', 'g5', 'PL', 'PR')]
    masses = {n for _, mm in props for n in re.findall(r'[A-Za-z_]\w*', str(mm))}
    moms = [x for x in moms if x not in masses]
    if 'l' not in moms:
        moms.append('l')
    by_rank = {}
    for c, items in _items(factors, moms):
        new, n = [], 0
        for it in items:
            if it == ('p', 'l'):
                n += 1
                new.append(('i', 'fsL%d' % n))
            else:
                new.append(it)
        by_rank.setdefault(n, []).append((c, new))
    acc = {}
    for n, tl in by_rank.items():
        num = " ".join("l^fsL%d" % (j + 1) for j in range(n)) or "1"
        res = loop(num, *props, kin=kin)
        for label, dct in res.raw:
            ren, subs_p = {}, {}
            for mm in _LABEL.finditer(label):
                if mm.group(1):
                    ren[mm.group(2)] = mm.group(1)
                else:
                    subs_p[mm.group(4)] = mm.group(3)
            for c, items in tl:
                its = [('p', subs_p[it[1]]) if it[0] == 'i' and it[1] in subs_p else
                       ('i', ren[it[1]]) if it[0] == 'i' and it[1] in ren else it for it in items]
                for k, cc in _core([(c, its)], p1, p2, M1, M2, spd).items():
                    slot = acc.setdefault(k, {})
                    for mk, mc in dct.items():
                        slot[mk] = slot.get(mk, SR(0)) + cc * _SR(mc)
    out = {k: _finalize_sr(v) for k, v in acc.items()}
    if gordon and left is not None:
        out = _gordon(out, p1, p2, M1, M2)
    return LineResult(_collect(out), left or ('', '', ''), right or ('', '', ''))


def transverse(T, v, vsq=None):
    r"""
    Package-X's Transverse: write a rank-2 tensor as T^{mu nu} = (g^{mu nu} - v^mu v^nu/v^2) A + (v^mu v^nu/v^2) B
    and return A.  T is a loop() result (or a dict of its coefficients) with the structures 'g^{mu nu}' and
    'v^mu v^nu'; vsq is v^2 (default: the symbol of loop(), e.g. p2 for p).
    """
    a, b = _rank2(T, v)
    return a


def longitudinal(T, v, vsq=None):
    """Package-X's Longitudinal: B in T = (g - v v/v^2) A + (v v/v^2) B, i.e. A + v^2 b."""
    a, b = _rank2(T, v)
    vsq = SR.var(v + '2') if vsq is None else sr(vsq)
    return (a + vsq * b).expand()


def _rank2(T, v):
    coeffs = T.coefficients() if hasattr(T, 'coefficients') else dict(T)
    a = b = SR(0)
    for k, c in coeffs.items():
        if re.fullmatch(r'(?:g|delta)\^\{\w+ \w+\}', k):
            a += c
        elif re.fullmatch(r'%s\^\w+ %s\^\w+' % (v, v), k):
            b += c
        else:
            raise ValueError("structure %s: only g^{mu nu} and %s^mu %s^nu are expected" % (k, v, v))
    return a, b


def to_chiral(result):
    """Write 1 and g5 through PL and PR: 1 = PL + PR, g5 = PR - PL (Package-X's ChiralBasis -> True)."""
    out = {}
    for (st, ch), c in (result.parts if hasattr(result, 'parts') else result).items():
        for nch, f in {'1': [('L', 1), ('R', 1)], '5': [('R', 1), ('L', -1)]}.get(ch, [(ch, 1)]):
            out[(st, nch)] = out.get((st, nch), SR(0)) + f * c
    return _collect(out)


def to_g5(result):
    """Write PL and PR through 1 and g5: PL = (1 - g5)/2, PR = (1 + g5)/2."""
    out = {}
    for (st, ch), c in (result.parts if hasattr(result, 'parts') else result).items():
        for nch, f in {'L': [('1', SR(1) / 2), ('5', -SR(1) / 2)], 'R': [('1', SR(1) / 2), ('5', SR(1) / 2)]}.get(ch, [(ch, 1)]):
            out[(st, nch)] = out.get((st, nch), SR(0)) + f * c
    return _collect(out)


# ---------------------------------------------------------------------------- traces and projectors
def _contract_coeff(c, items):
    """Sum an index carried by a vector component in the coefficient (p__mu) with a gamma^mu in the
    product (it becomes p/) or with another component (a dot product).  Returns [(c, items)]."""
    out = [(SR(c).expand(), list(items))]
    changed = True
    while changed:
        changed = False
        new = []
        for c, its in out:
            idxs = [x[1] for x in its if x[0] == 'i']
            comps = [v for v in c.variables() if '__' in str(v) and not str(v).startswith('g__')]
            hit = None
            for v in comps:
                vec, mu = str(v).split('__', 1)
                if mu in idxs:
                    hit = (v, vec, mu)
                    break
            if hit is None:
                new.append((c, its))
                continue
            v, vec, mu = hit
            cv = c.coefficient(v)
            rest = (c - cv * v).expand()
            j = its.index(('i', mu))
            new.append((cv, its[:j] + [('p', vec)] + its[j + 1:]))
            if not rest.is_zero():
                new.append((rest, its))
            changed = True
        out = new
    return out


def _trace_slots(s, order):
    """Trace of an ordered product of gamma slots (d dimensions, Tr 1 = 4)."""
    if not s:
        return SR(4)
    if len(s) % 2:
        return SR(0)
    a = s[0]
    tot = SR(0)
    for k in range(1, len(s)):
        tot += (-1) ** (k + 1) * order.half_anti(a, s[k]) * _trace_slots(s[1:k] + s[k + 1:], order)
    return tot


def trace(factors, momenta, sp=None):
    """Tr of a product of Dirac matrices (factors as in line_expand, also 'p^mu' components), d-dimensional,
    gamma_5 anticommuting; a single gamma_5 with fewer than four other matrices gives 0."""
    order = _Order(sp={tuple(sorted(k)): sr(v) for k, v in (sp or {}).items()})
    tot = SR(0)
    for c, items in _items(factors, momenta):
        for c2, its in _contract_coeff(c, items):
            c3, slots, ch = _chiral_right(c2, its)
            if c3.is_zero():
                continue
            for s, cc in _normal_order(c3, slots, order).items():
                t = _trace_slots(list(s), order)
                if ch == '1':
                    tot += cc * t
                elif ch in ('L', 'R'):
                    if len(s) >= 4:
                        raise NotImplementedError("trace with gamma_5 and four or more gamma matrices: use form.dirac_trace")
                    tot += cc * t / 2
                elif len(s) >= 4:
                    raise NotImplementedError("trace with gamma_5 and four or more gamma matrices: use form.dirac_trace")
    return tot.expand()


def _vertex_basis(mu, p1, m1, p2, m2, equal):
    """Package-X's on-shell vector/axial-vector vertex structures, as lists of (coefficient, factors);
    q = p2 - p1, written q in the factors and qsq = q^2."""
    qsq = SR.var('fs_qsq')
    M = 2 * m1 if equal else m1 + m2
    qmu = "(%s__%s - %s__%s)" % (p2, mu, p1, mu)
    qs = "%s - %s" % (p2, p1)
    trans = [(SR(1), [mu]), (-1 / qsq, [qmu, qs])] if not equal else [(SR(1), [mu])]
    # i sigma^{mu nu} q_nu = -(g^mu q/ - q/ g^mu)/2
    sig = [(-SR(1) / 2 / M, [mu, qs]), (SR(1) / 2 / M, [qs, mu])]
    lon = [((2 if not equal else 1) / (M if not equal else m1), [qmu])]
    g5 = lambda t: [(c, f + ['g5']) for c, f in t]
    trans5 = [(SR(1), [mu]), (-1 / qsq, [qmu, qs])]
    return {'F1': trans, 'F2': sig, 'F3': lon, 'G1': g5(trans5), 'G2': g5(sig), 'G3': g5(lon)}


_PROJ_CACHE = {}


def _inverse(rows):
    """Inverse of a matrix of symbolic rational functions, done in a fraction field (the imaginary unit
    as a formal symbol, put back at the end), which is much faster than in SR."""
    from sage.all import PolynomialRing, QQ, matrix
    Ig = SR.var('fs_Ig')
    rows = [[SR(x).subs({I: Ig}) for x in r] for r in rows]
    names = sorted({str(v) for r in rows for x in r for v in x.variables()} | {'fs_Ig'})
    F = PolynomialRing(QQ, names).fraction_field()
    M = matrix(F, [[F(str(x)) for x in r] for r in rows])
    Mi = M.inverse()
    return [[SR(Mi[i, j]).subs({Ig: I}) for j in range(len(rows))] for i in range(len(rows))]


_ALIAS = {'Kinetic': 'A', 'Mass': 'B', 'AxialKinetic': 'C', 'ImaginaryMass': 'E', 'Scalar': 'S',
          'Pseudoscalar': 'P', 'Dirac': 'F1', 'Pauli': 'F2', 'Anapole': 'G1', 'EDM': 'G2'}


def _chi(t, ch):
    return [(c, f + [ch]) for c, f in t]


def projector(name, mu=None, p1=None, m1=None, p2=None, m2=None, q2='q2', p=None, m=None):
    r"""
    Package-X's Projector: P with Tr[M P] = the form factor of the Dirac matrix M (spur(), form_factor()).

    Vertex, u-bar(p2) M^mu u(p1) on shell, q = p2 - p1, q2 = q^2:
        projector(name, 'mu', 'p1', m1, 'p2', m2)
        M^mu = (g^mu - q/ q^mu/q^2) F1 + i sigma^{mu nu} q_nu/(m1 + m2) F2 + 2 q^mu/(m1 + m2) F3
             + (g^mu - q/ q^mu/q^2) g5 G1 + i sigma^{mu nu} g5 q_nu/(m1 + m2) G2 + 2 q^mu/(m1 + m2) g5 G3
        (m1 = m2 = m: g^mu F1 + i sigma q/(2m) F2 + q^mu/m F3 + ...), names F1 'Dirac', F2 'Pauli', F3,
        G1 'Anapole', G2 'EDM', G3, 'SachsElectric' = F1 + q^2/(4m^2) F2, 'SachsMagnetic' = F1 + F2, and
        the chiral set AL, BL, CL, AR, BR, CR (the same three structures times PL and PR).
    Scalar density, u-bar(p2) M u(p1):  projector('S' or 'P', None, 'p1', m1, 'p2', m2),  M = GS + i g5 GP.
    Self-energy, off shell:  projector(name, p='p', m=m),  M = p/ A + m B + p/ g5 C + i m g5 E,
        names A 'Kinetic', B 'Mass', C 'AxialKinetic', E 'ImaginaryMass', and AL, BL, AR, BR for
        M = (p/ AL + m BL) PL + (p/ AR + m BR) PR.
    Built by linear algebra on the structures (Tr[B_j P_i] = delta_ij), so it is exact in d.
    """
    from sage.all import matrix
    name = _ALIAS.get(name, name)
    key = (name, mu, p1, str(m1), p2, str(m2), q2, p, str(m))
    if key in _PROJ_CACHE:
        return _PROJ_CACHE[key]
    if p is not None:                                       # self-energy
        mm = sr(m)
        sp = {(p, p): SR.var('%s2' % p)}
        mom = [p]
        if name in ('A', 'B', 'C', 'E'):
            B = {'A': [(SR(1), [p])], 'B': [(mm, [])], 'C': [(SR(1), [p, 'g5'])], 'E': [(I * mm, ['g5'])]}
            block = ['A', 'B', 'C', 'E']
        else:
            B = {'AL': [(SR(1), [p, 'PL'])], 'BL': [(mm, ['PL'])], 'AR': [(SR(1), [p, 'PR'])], 'BR': [(mm, ['PR'])]}
            block = ['AL', 'BL', 'AR', 'BR']
        wrap = lambda f: f
    else:
        m1, m2 = sr(m1), sr(m2)
        Q2 = sr(q2)
        sp = {(p1, p1): m1 ** 2, (p2, p2): m2 ** 2, (p1, p2): (m1 ** 2 + m2 ** 2 - Q2) / 2}
        mom = [p1, p2]
        wrap = lambda f: ["%s + %s" % (p1, m1)] + f + ["%s + %s" % (p2, m2)]
        if mu is None:
            B = {'S': [(SR(1), [])], 'P': [(I, ['g5'])]}
            block = ['S', 'P']
        else:
            V = _vertex_basis(mu, p1, m1, p2, m2, bool((m1 - m2).is_zero()))
            if name in ('SachsElectric', 'SachsMagnetic'):
                PF1 = projector('F1', mu, p1, m1, p2, m2, q2)
                PF2 = projector('F2', mu, p1, m1, p2, m2, q2)
                k2 = Q2 / (4 * m1 ** 2) if name == 'SachsElectric' else SR(1)
                out = PF1 + [(c * k2, f) for c, f in PF2]
                _PROJ_CACHE[key] = out
                return out
            if name in ('AL', 'BL', 'CL', 'AR', 'BR', 'CR'):
                base = _vertex_basis(mu, p1, m1, p2, m2, False)
                B = {}
                for ch in ('L', 'R'):
                    for a_, f_ in (('A', 'F1'), ('B', 'F2'), ('C', 'F3')):
                        B[a_ + ch] = _chi(base[f_], 'P' + ch)
                block = ['AL', 'BL', 'CL', 'AR', 'BR', 'CR']
            else:
                B = V
                block = ['F1', 'F2', 'F3'] if name in ('F1', 'F2', 'F3') else ['G1', 'G2', 'G3']
    if name not in block:
        raise ValueError("unknown projector %r" % name)
    qsq = SR.var('fs_qsq')

    def tr(Bj, Rk):
        tot = SR(0)
        for c1, f1 in Bj:
            for c2, f2 in Rk:
                tot += c1 * c2 * trace(f1 + wrap(f2), mom, sp)
        return tot.subs({qsq: sr(q2)}) if p is None else tot
    n = len(block)
    X = _inverse([[tr(B[block[j]], B[block[k]]) for j in range(n)] for k in range(n)])
    i = block.index(name)
    P = []
    for k in range(n):
        coef = X[i][k]
        if coef.is_zero():
            continue
        for c, fl in B[block[k]]:
            cc = (coef * c).subs({qsq: sr(q2)}) if p is None else coef * c
            P.append((cc, wrap(fl)))
    _PROJ_CACHE[key] = P
    return P


def form_factor(factors, P, *props, kin=None):
    r"""
    Int Tr[M(l) P] / prod D_i: a form factor of a one-loop vertex, the projection done before the
    integral (exact in d), as Package-X's LoopIntegrate of a Spur with a Projector.
        P = projector('F2', 'mu', 'p1', 'm', 'p2', 'm')
        form_factor(['rho', 'p2 - l + m', 'mu', 'p1 - l + m', 'rho'], P, ['l', '0'], ['l - p2', 'm'],
                    ['l - p1', 'm'], kin={'p1^2': 'm^2', 'p2^2': 'm^2', 'p1.p2': 'm^2 - q2/2'})
    """
    from .pv import loop
    spd = _kin_sp(kin)
    moms = sorted({n for q, _ in props for n in re.findall(r'[A-Za-z_]\w*', str(q))} | {a for k in spd for a in k} | {'l'})
    masses = {n for _, mm in props for n in re.findall(r'[A-Za-z_]\w*', str(mm))}
    moms = [x for x in moms if x not in masses]
    t = SR(0)
    for c, fl in P:
        t += c * trace(list(factors) + fl, moms, {})        # keep the l products symbolic
    t = t.expand()
    num = str(t)
    num = re.sub(r'\b(\w+)_dot_(\w+)\b', r'\1.\2', num)
    num = re.sub(r'\b(\w+)_sq\b', r'\1^2', num)
    return loop(num, *props, kin=kin)


def spur(M, P, momenta, sp=None):
    """Tr[M P] for two lists of (coefficient, factors), with the free indices of M and P summed
    (they must use the same names)."""
    tot = SR(0)
    for c1, f1 in M:
        for c2, f2 in P:
            tot += c1 * c2 * trace(f1 + f2, momenta, sp)
    return tot.expand()


def chisholm(result):
    r"""
    Four-dimensional Chisholm identity on every ordered product of three gamma matrices (Package-X's
    ChisholmExpand), with eps^{0123} = +1 and g5 = i g^0 g^1 g^2 g^3 as in Peskin and Schroeder and
    Package-X (checked on explicit matrices):
        g^a g^b g^c = g^{ab} g^c + g^{bc} g^a - g^{ac} g^b + i eps^{a b c s} g_s g5.
    The last term is the structure ('eps', a, b, c) = eps^{a b c s} g_s, with its chirality times g5.
    Only valid in four dimensions: use it on finite results or after the eps expansion.
    """
    order = _Order()
    out = {}
    parts = result.parts if hasattr(result, 'parts') else result
    for (st, ch), c in parts.items():
        if len(st) != 3 or st[0] == 'sigma':
            out[(st, ch)] = out.get((st, ch), SR(0)) + c
            continue
        a, b, cc = st
        for coef, single in ((order.half_anti(a, b), cc), (order.half_anti(b, cc), a), (-order.half_anti(a, cc), b)):
            k = ((single,), ch)
            out[k] = out.get(k, SR(0)) + c * coef
        sg, nch = _chiral_mul('5', ch)
        if nch != '0':
            k = (('eps', a, b, cc), nch)
            out[k] = out.get(k, SR(0)) + sg * I * c
    out = _collect(out)
    return LineResult(out, result.left, result.right) if hasattr(result, 'parts') else out


# ---------------------------------------------------------------------------- products of fermion lines
def _slot_key(s):
    return (0 if s[0] == 'i' else 1, str(s[1]))


def _asym_sort(slots):
    """Sort the slots of an antisymmetric Gamma^{[...]}: (sign, sorted tuple), sign 0 if a slot repeats."""
    sl = list(slots)
    if len(set(sl)) < len(sl):
        return 0, ()
    sign = 1
    for i in range(len(sl)):                                   # bubble sort, counting swaps
        for j in range(len(sl) - 1 - i):
            if _slot_key(sl[j]) > _slot_key(sl[j + 1]):
                sl[j], sl[j + 1] = sl[j + 1], sl[j]
                sign = -sign
    return sign, tuple(sl)


def _to_asym(slots):
    r"""An ordered product g^{a1} ... g^{an} in the d-dimensional antisymmetric basis:
    Gamma^{[A]} g^b = Gamma^{[A b]} + sum_i (-1)^(k-1-i) g^{a_i b} Gamma^{[A without a_i]}.
    Returns [(coefficient, [metric pairs], antisymmetric slot tuple)]."""
    terms = [(1, [], ())]
    for b in slots:
        new = []
        for c, gs, A in terms:
            if b not in A:
                new.append((c, gs, A + (b,)))
            k = len(A)
            for i, a in enumerate(A):
                new.append((c * (-1) ** (k - 1 - i), gs + [(a, b)], A[:i] + A[i + 1:]))
        terms = new
    return terms


def _coef_tensors(c):
    """A coefficient from line_expand -> [(scalar, [metric pairs])]: g__a__b -> (a, b), p__a -> (p, a)."""
    out = []
    c = SR(c).expand()
    terms = c.operands() if c.operator() is not None and 'add' in getattr(c.operator(), '__name__', '') else [c]
    for t in terms:
        gs, scal = [], t
        for v in t.variables():
            n = str(v)
            if '__' not in n:
                continue
            parts = n.split('__')
            deg = t.degree(v)
            if parts[0] == 'g' and len(parts) == 3:
                pair = (('i', parts[1]), ('i', parts[2]))
            elif len(parts) == 2:
                pair = (('p', parts[0]), ('i', parts[1]))
            else:
                continue
            gs += [pair] * int(deg)
            scal = scal / v ** deg
        out.append((scal, gs))
    return out


def _metric_value(a, b):
    """g between two slots that are not summed: a symbol."""
    if a[0] == 'p' and b[0] == 'p':
        return dot(a[1], b[1])
    if a[0] == 'p':
        return comp(a[1], b[1])
    if b[0] == 'p':
        return comp(b[1], a[1])
    return metric(a[1], b[1])


_DUMMY = [0]


def _fresh():
    _DUMMY[0] += 1
    return ('i', 'fsd%d' % _DUMMY[0])


def _chisholm_line(G, ch):
    """Gamma^{[A]} times the chiral factor, in four dimensions: rank 3 -> i eps^{abcs} g_s g5,
    rank 4 -> -i eps^{abce} g5, rank >= 5 -> 0.  Returns [(coef, eps slots or None, new A, new chiral)]."""
    n = len(G)
    if n <= 2:
        return [(1, None, G, ch)]
    if n >= 5:
        return []
    flip = {'1': ('5', 1), '5': ('1', 1), 'L': ('L', -1), 'R': ('R', 1)}       # g5 X
    ch2, s = flip[ch]
    if n == 3:
        dm = _fresh()
        return [(I * s, G + (dm,), (dm,), ch2)]
    return [(-I * s, G, (), ch2)]


def _eps_eps(A, B):
    """eps^{A} eps^{B} = -det[g(A_i, B_j)] (eps^0123 = +1): [(sign, metric pairs)]."""
    import itertools
    out = []
    for perm in itertools.permutations(range(4)):
        inv = sum(1 for i in range(4) for j in range(i + 1, 4) if perm[i] > perm[j])
        out.append((-(-1) ** inv, [(A[i], B[perm[i]]) for i in range(4)]))
    return out


def _contract(term):
    """Sum every index that appears twice.  term = (coef, metrics, epses, lines) with lines a list of
    [antisymmetric slots, chiral].  Returns a list of terms with the remaining contractions kept."""
    todo, done = [term], []
    while todo:
        c, gs, es, lines = todo.pop()
        if c == 0:
            continue
        where = {}
        for k, (a, b) in enumerate(gs):
            for s in (a, b):
                if s[0] == 'i':
                    where.setdefault(s, []).append(('g', k))
        for k, e in enumerate(es):
            for s in e:
                if s[0] == 'i':
                    where.setdefault(s, []).append(('e', k))
        for k, (A, ch) in enumerate(lines):
            for s in A:
                if s[0] == 'i':
                    where.setdefault(s, []).append(('l', k))
        act = None
        for s, occ in where.items():
            if any(o[0] == 'g' for o in occ) and len(occ) >= 2 or (len(occ) == 2 and occ[0] == occ[1] and occ[0][0] == 'g'):
                act = ('g', s, occ)
                break
        if act is None:
            for s, occ in where.items():
                kinds = [o[0] for o in occ]
                if len(occ) == 2 and kinds == ['e', 'e'] and occ[0][1] != occ[1][1]:
                    act = ('ee', s, occ)
                    break
        if act is None:
            for s, occ in where.items():
                if len(occ) == 2 and occ[0] == occ[1]:              # twice in one antisymmetric object: zero
                    act = ('zero', s, occ)
                    break
        if act is None:
            done.append((c, gs, es, lines))
            continue
        kind, s, occ = act
        if kind == 'zero':
            continue
        if kind == 'g':
            k = next(o[1] for o in occ if o[0] == 'g')
            a, b = gs[k]
            rest = gs[:k] + gs[k + 1:]
            if a == s and b == s:
                todo.append((c * D, rest, es, lines))
                continue
            other = b if a == s else a
            rep = lambda x: other if x == s else x
            if any(s in pair for pair in rest) or any(s in e for e in es) or any(s in A for A, _ in lines):
                rest = [(rep(x), rep(y)) for x, y in rest]
                es2 = [tuple(rep(x) for x in e) for e in es]
                lines2 = [(tuple(rep(x) for x in A), ch) for A, ch in lines]
                cc = c
                ok = True
                es3 = []
                for e in es2:
                    if len(set(e)) < 4:
                        ok = False
                    es3.append(e)
                lines3 = []
                for A, ch in lines2:
                    sg, A2 = _asym_sort(A)
                    if sg == 0:
                        ok = False
                    cc *= sg
                    lines3.append((A2, ch))
                if ok:
                    todo.append((cc, rest, es3, lines3))
            else:
                done.append((c, rest + [(a, b)], es, lines))                 # nothing to contract with
            continue
        k1, k2 = occ[0][1], occ[1][1]
        A, B = es[k1], es[k2]
        rest = [e for i, e in enumerate(es) if i not in (k1, k2)]
        for sg, pairs in _eps_eps(A, B):
            todo.append((c * sg, gs + pairs, rest, lines))
    return done


def _metric_scalar(gs):
    """Metric pairs with no summed index left -> (scalar, kept pairs whose index is still summed elsewhere)."""
    v = SR(1)
    for a, b in gs:
        v *= _metric_value(a, b)
    return v


def _canonical(lines, es):
    """Rename the summed indices to fs1, fs2, ... in the order that gives the smallest key."""
    import itertools
    count = {}
    for A, _ in lines:
        for s in A:
            if s[0] == 'i':
                count[s] = count.get(s, 0) + 1
    for e in es:
        for s in e:
            if s[0] == 'i':
                count[s] = count.get(s, 0) + 1
    dummies = sorted([s for s, n in count.items() if n >= 2], key=_slot_key)
    best = None
    perms = itertools.permutations(range(len(dummies))) if len(dummies) <= 5 else [tuple(range(len(dummies)))]
    for perm in perms:
        ren = {d: ('i', 'fs%d' % (perm[i] + 1)) for i, d in enumerate(dummies)}
        sign = 1
        newl = []
        for A, ch in lines:
            sg, A2 = _asym_sort(tuple(ren.get(x, x) for x in A))
            sign *= sg
            newl.append((A2, ch))
        newe = []
        for e in es:
            sg, e2 = _asym_sort(tuple(ren.get(x, x) for x in e))
            sign *= sg
            newe.append(e2)
        key = (tuple(newl), tuple(sorted(newe)))
        if best is None or str(key) < str(best[0]):
            best = (key, sign)
    return best


class LineProduct:
    """What line_product() returns: {((structure, chiral), (structure, chiral), ...): coefficient}."""

    def __init__(self, parts, lines):
        self.parts, self.lines = parts, lines

    def items(self):
        return self.parts.items()

    def coefficients(self):
        return dict(self.parts)

    def __getitem__(self, k):
        return self.parts[k]

    def __len__(self):
        return len(self.parts)

    def __repr__(self):
        rows = []
        for key, c in sorted(self.parts.items(), key=lambda t: str(t[0])):
            if key == 'eps':
                continue
            rows.append("%s  :  %s" % (" (x) ".join(structure_name(s, ch) for s, ch in key), c))
        return "\n".join(rows) if rows else "0"


def line_product(*lines, chisholm=True, gordon=True, on_shell=True, sp=None, momenta=None):
    r"""
    The product of fermion lines with Lorentz indices summed between them (Package-X's
    FermionLineProduct with FermionLineExpand), e.g. for four-fermion operators:

        line_product((("u", "p3", "m"), ["mu", "nu", "rho"], ("u", "p1", "m")),
                     (("u", "p4", "m"), ["mu", "nu", "rho"], ("u", "p2", "m")))

    Each line is (left spinor, factors, right spinor) as in line_expand, which reduces it first (Dirac
    equation, Gordon identities).  An index that appears in two lines is summed.  Every line is then
    written in the d-dimensional antisymmetric basis Gamma^{[a1...an]} and the indices are summed in
    d dimensions.  chisholm=True (Package-X's ChisholmExpand) then puts each line in the
    four-dimensional SVTAP basis: Gamma^{[abc]} = i eps^{abcs} g_s g5 and Gamma^{[abce]} = -i eps^{abce} g5
    (eps^0123 = +1), ranks five and up vanish, and two eps tensors contract with the d-dimensional
    metric, eps^{abcs} eps_{abct} = -(d-1)(d-2)(d-3) g^s_t.  At d = 4 the result is the usual
    four-dimensional identity; away from d = 4 it differs from Package-X's for some index orders
    (evanescent terms, a choice of scheme).  chisholm=False keeps the d-dimensional antisymmetric
    products, written ('asym', a, b, c, ...).

    Returns a LineProduct: {((structure, chiral) for each line): coefficient}, with the summed indices
    renamed fs1, fs2, ... (structures as in line_expand; sigma^{ab} = i Gamma^{[ab]}).
    """
    expanded = [line_expand(l, f, r, sp=sp, momenta=momenta, gordon=gordon, on_shell=on_shell) for l, f, r in lines]
    # every line: [(coefficient, metric pairs, antisymmetric slots, chiral)]
    per = []
    for res in expanded:
        items = []
        for (st, ch), c in res.items():
            if st == ():
                base = [(1, [], ())]
            elif st[0] == 'sigma':
                base = [(I, [], (st[1], st[2]))]                       # sigma^{ab} = i Gamma^{[ab]}
                sg, A = _asym_sort(base[0][2])
                base = [(I * sg, [], A)] if sg else []
            else:
                base = []
                for cc, gs, A in _to_asym(st):
                    sg, A2 = _asym_sort(A)
                    if sg:
                        base.append((cc * sg, gs, A2))
            for sc, gs0 in _coef_tensors(c):
                for cc, gs, A in base:
                    items.append((sc * cc, gs0 + gs, A, ch))
        per.append(items)
    import itertools
    out = {}
    for combo in itertools.product(*per):
        c = SR(1)
        gs, lines_ = [], []
        for sc, g_, A, ch in combo:
            c *= sc
            gs += g_
            lines_.append((A, ch))
        starts = [(c, gs, [], lines_)]
        if chisholm:
            new = []
            for c0, g0, e0, l0 in starts:
                options = [_chisholm_line(A, ch) for A, ch in l0]
                for pick in itertools.product(*options):
                    cc, es, ls = c0, list(e0), []
                    for f, eps, A, ch in pick:
                        cc *= f
                        if eps is not None:
                            es.append(eps)
                        ls.append((A, ch))
                    new.append((cc, list(g0), es, ls))
            starts = new
        for st in starts:
            for c2, g2, e2, l2 in _contract(st):
                scal = c2 * _metric_scalar(g2)
                key, sign = _canonical(l2, e2)
                lines_key, eps_key = key
                kk = tuple(_struct_of(A, ch) for A, ch in lines_key)
                factor = sign
                for A, ch in lines_key:
                    if len(A) == 2:
                        factor *= -I                                 # Gamma^{[ab]} = -i sigma^{ab}
                if eps_key:
                    from .form import EPS
                    for e in eps_key:
                        factor *= EPS(*[SR.var(str(x[1])) for x in e])
                full = kk
                out[full] = out.get(full, SR(0)) + (scal * factor).expand()
    out = {k: v for k, v in out.items() if not v.is_zero()}
    return LineProduct(out, lines)


def _struct_of(A, ch):
    if len(A) == 0:
        return ((), ch)
    if len(A) == 1:
        return ((A[0],), ch)
    if len(A) == 2:
        return (('sigma', A[0], A[1]), ch)
    return (('asym',) + tuple(A), ch)
