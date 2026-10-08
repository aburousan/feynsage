r"""
The operations a QFT calculation needs, done by FORM:

    dirac_trace(expr, dim=4, rules=None)     Tr[...] of a product of Dirac matrices
    contract(expr, rules=None)               sum repeated Lorentz indices, put in kinematics
    conjugate(M)                             complex conjugate of an amplitude (spinor bilinears)
    dirac_bar(G)                             gamma^0 G^dagger gamma^0
    spin_sum(M * conjugate(M))               sum over all spins: u ubar -> p/ + m, v vbar -> p/ - m

dim = 4 (default) or a symbol (d, D, ...).  In 4 dimensions gamma5 is traced by FORM (trace4).  In d
dimensions gamma5 anticommutes with every gamma matrix (NDR, the convention of feynsage.dirac and
Package-X): an even number of gamma5 drops out, an odd number with fewer than four other gamma
matrices gives 0, and anything else is ambiguous in NDR and raises an error (use dim=4).
Every repeated Lorentz index is summed.  rules = {dot(p, p): m^2, dot(p, k): s/2, p4: p1 + p2 - p3,
m: 0} are put in by FORM (momentum rules before the trace, the others after the contraction).
debug=True prints the FORM program and FORM's raw output.
"""
from sage.all import SR, I, Integer
from . import tensors as T
from .objects import DiracExpr, Factor, Chain, Spinor, DiracError, _combine
from .compiler import Program, CompileError, _opname
from .parser import Reader, expression_text
from .backend import backend


# ---------------------------------------------------------------------------- running a program
def _execute(prog, exprs, lines_per_expr=0, vec_ids=(), post_ids=(), trace4=True, contract=True,
             debug=False, threads=1):
    """exprs: list of FORM texts for Local E1, E2, ...; returns the printed bodies."""
    src = [prog.declarations(), 'Off statistics;']
    for i, e in enumerate(exprs, 1):
        src.append('Local FSE%d = %s;' % (i, e))
    src += list(vec_ids)
    if not trace4 and any('g5_' in e for e in exprs):
        raise DiracError("internal: gamma5 reached a d-dimensional trace")      # never send g5_ to tracen
    for n in range(1, lines_per_expr + 1):
        src.append(('trace4,%d;' if trace4 else 'tracen,%d;') % n)
    if contract:
        src.append('contract;')
    if prog.has_color:
        from .color import FORM_RULES
        src.append(FORM_RULES)
    if post_ids:
        src += list(post_ids)
        src.append('.sort')
        src += list(post_ids)          # a second pass for rules whose right side contains other rules
    src += ['Format nospaces;', 'Print;', '.end', '']
    code = '\n'.join(src)
    if debug:
        print('--- FORM program ---\n' + code)
    out = backend().run(code, threads=threads)
    if debug:
        print('--- FORM output ---\n' + out)
    return [expression_text(out, 'FSE%d' % i) for i in range(1, len(exprs) + 1)]


def _check(prog, dim, eps_in_d, gamma5_scheme="NDR"):
    if gamma5_scheme not in ("NDR", "NDR-even"):
        raise NotImplementedError('only gamma5_scheme="NDR" (anticommuting gamma5) and "NDR-even" are implemented; '
                                  'in 4 dimensions (dim=4) FORM traces gamma5 exactly and no scheme is needed')
    if dim != 4 and prog.has_eps and not eps_in_d:
        raise DiracError("a Levi-Civita tensor is contracted in %s dimensions.  epsilon is a 4-dimensional object; "
                         "contracting it in d dimensions is a scheme choice (FORM uses d-dimensional metrics in "
                         "eps.eps).  Use dim=4, or pass eps_in_d=True if this is what you want." % dim)


def _threads(threads):
    if threads is None:
        import os
        try:
            return int(os.environ.get('FEYNSAGE_FORM_THREADS', '1'))
        except ValueError:
            return 1
    return threads


# ---------------------------------------------------------------------------- gamma5 in d dimensions
def _ndr(coef, factors, drop_odd=False):
    """Anticommuting gamma5 in d dimensions: [(c, factors without gamma5)] equal to the trace of
    coef * product(factors), or raise when the trace with one gamma5 left over is ambiguous.
    drop_odd: leave out the traces with one gamma5 (gamma5_scheme="NDR-even", for quantities whose
    parity-odd part is known to vanish, e.g. a parity-even form factor with too few momenta for an
    epsilon tensor)."""
    branches = [(SR(coef), [], 0)]                 # (coefficient, factors so far, gamma5 parity)
    for F in factors:
        c5 = F.gamma5_coefficient()
        rest = F.without_gamma5()
        new = []
        for c, fs, par in branches:
            if not rest.is_zero():
                new.append((c, fs + [rest.flip() if par else rest], par))
            if not c5.is_trivial_zero():
                new.append((c * c5, fs, 1 - par))
        branches = new
    out = []
    for c, fs, par in branches:
        if par == 0:
            out.append((c, fs))
            continue
        ngam = sum(1 for F in fs if F.odd_even()[1].items)
        if ngam < 4 or drop_odd:
            continue                                # Tr[gamma5 x (fewer than 4 gammas)] = 0
        raise DiracError("a trace with gamma5 and %d other gamma matrices in d dimensions is ambiguous with an "
                         "anticommuting gamma5 (NDR).  Use dim=4, or do the gamma5 part in 4 dimensions." % ngam)
    return out



# ---------------------------------------------------------------------------- momentum rules
def _vector_rules(rules):
    """The rules that replace a momentum (p4: p1 + p2 - p3), and the others."""
    vr, other = {}, {}
    for k, val in (rules or {}).items():
        if T.is_momentum(SR(k)):
            T.vector_parts(val, "the right side of a momentum rule")
            vr[SR(k)] = SR(val)
        else:
            other[k] = val
    return vr, other


def _canonical_tensors(e):
    """Dot, Comp and Epsilon with sums of momenta in their slots -> expanded and in canonical order."""
    e = SR(e)
    w = [SR.wild(n) for n in range(4)]
    sub = {}
    for x in e.find(T.DOT(w[0], w[1])):
        sub[x] = T.dot(*x.operands())
    for x in e.find(T.COMP(w[0], w[1])):
        sub[x] = T.comp(*x.operands())
    for x in e.find(T.EPS(*w)):
        sub[x] = T.epsilon(*x.operands())
    return e.subs(sub) if sub else e


def _vr_factor(F, vr):
    if not vr:
        return F
    items = []
    for a, c in F.items:
        c = _canonical_tensors(c.subs(vr))
        if a[0] == 'p' and SR(a[1]) in vr:
            items += [(('p', v), c * cv) for v, cv in T.vector_parts(vr[SR(a[1])])]
        else:
            items.append((a, c))
    return Factor(items)


def _apply_vector_rules(expr, vr):
    """Put momentum rules in everywhere, also in denominators (FORM's id does not reach 1/p.k)."""
    if not vr:
        return expr
    if not isinstance(expr, DiracExpr):
        return _canonical_tensors(SR(expr).subs(vr))

    fac = lambda F: _vr_factor(F, vr)

    def chain(ch):
        return Chain(ch.left, tuple(fac(F) for F in ch.factors), ch.right) if ch is not None else None
    return DiracExpr([(_canonical_tensors(SR(c).subs(vr)), chain(o), tuple(chain(x) for x in cl))
                      for c, o, cl in expr.terms])

# ---------------------------------------------------------------------------- traces
def _line_texts(prog, coef, lines, dim4, drop_odd=False):
    """A term coef * Tr(line_1) * Tr(line_2) ... -> FORM text (NDR for gamma5 in d dims)."""
    alts = [(SR(coef), [])]
    for n, factors in enumerate(lines, 1):
        if dim4:
            alts = [(c, ch + [(n, list(factors))]) for c, ch in alts]
        else:
            alts = [(c * c2, ch + [(n, f2)]) for c, ch in alts for c2, f2 in _ndr(1, factors, drop_odd)]
    texts = []
    for c, ch in alts:
        parts = [prog.chain(f, n) for n, f in ch]
        cs = prog.scalar(c)
        texts.append('*'.join(([cs] if cs != '1' else []) + parts) or '1')
    return texts


def _trace_terms(expr):
    if isinstance(expr, DiracExpr):
        k = expr.kind()
        if k not in ('matrix', None):
            raise DiracError("dirac_trace needs a product of Dirac matrices without spinors; this is a %s "
                             "expression.  For |M|^2 summed over spins use spin_sum(M * conjugate(M))." % k)
        out = []
        for c, o, cl in expr.terms:
            if cl:
                raise DiracError("dirac_trace: the expression contains spinor bilinears; use spin_sum")
            out.append((c, [o.factors]))
        return out
    return [(SR(expr), [()])]                       # a number times the unit matrix


def dirac_trace(expr, dim=4, rules=None, contract=True, euclidean=None, gamma5_scheme="NDR", eps_in_d=False,
                debug=False, threads=None):
    r"""
    The Dirac trace Tr[expr] (Tr 1 = 4), computed by FORM, as a Sage expression in dot(p, q),
    comp(p, mu), metric g(mu, nu) and the Levi-Civita Eps.  A list of expressions is traced in one
    FORM run and gives a list.

        p, q = momenta("p q"); mu, nu = lorentz_indices("mu nu")
        dirac_trace(slash(p) * gamma(mu) * slash(q) * gamma(nu))
          -> 4*comp(p, mu)*comp(q, nu) + 4*comp(p, nu)*comp(q, mu) - 4*dot(p, q)*g(mu, nu)
    """
    dim = T.normalize_dim(dim)
    batch = isinstance(expr, (list, tuple))
    vr, rules = _vector_rules(rules)
    items = [_apply_vector_rules(e, vr) for e in (list(expr) if batch else [expr])]
    eu = T._euclid(euclidean)
    prog = Program(dim, eu)
    dim4 = (dim == 4)
    texts, nlines = [], 0
    for e in items:
        _guard_indices(e)
        parts = []
        for c, lines in _trace_terms(e):
            parts += _line_texts(prog, c, lines, dim4, gamma5_scheme == "NDR-even")
            nlines = max(nlines, len(lines))
        texts.append('+'.join(parts) if parts else '0')
    _check(prog, dim, eps_in_d, gamma5_scheme)
    vec, post = prog.rules(rules)
    bodies = _execute(prog, texts, nlines, vec, post, trace4=dim4, contract=contract, debug=debug,
                      threads=_threads(threads))
    rd = Reader(prog)
    res = [rd.scalar(b) for b in bodies]
    return res if batch else res[0]


# ---------------------------------------------------------------------------- contraction
def contract(expr, rules=None, dim=4, euclidean=None, eps_in_d=False, debug=False, threads=None):
    r"""
    Sum every repeated Lorentz index of a Sage expression made of dot, comp, metric and epsilon
    (also the output of dirac_trace and of form.compute) and put in the rules, with FORM.

        contract(metric(mu, nu) * comp(p, mu) * comp(q, nu))        -> dot(p, q)
        contract(epsilon(mu, nu, rho, sigma) * epsilon(mu, nu, rho, sigma))   -> -24

    A DiracExpr is first contracted into its gamma matrices (g^{mu nu} gamma_nu -> gamma^mu,
    p_mu gamma^mu -> slash(p)) and then simplified with simplify_dirac.  A list gives a list.
    """
    dim = T.normalize_dim(dim)
    batch = isinstance(expr, (list, tuple))
    vr, rules = _vector_rules(rules)
    items = [_apply_vector_rules(e, vr) for e in (list(expr) if batch else [expr])]
    if any(isinstance(e, DiracExpr) for e in items):
        from .simplify import simplify_dirac
        res = [simplify_dirac(e, dim=dim, rules=rules, euclidean=euclidean, eps_in_d=eps_in_d) if isinstance(e, DiracExpr)
               else contract(e, rules=rules, dim=dim, euclidean=euclidean, eps_in_d=eps_in_d, debug=debug, threads=threads) for e in items]
        return res if batch else res[0]
    eu = T._euclid(euclidean)
    prog = Program(dim, eu)
    for e in items:
        _guard_indices(e)
    texts = [prog.scalar(SR(e)) for e in items]
    _check(prog, dim, eps_in_d)
    vec, post = prog.rules(rules)
    bodies = _execute(prog, texts, 0, vec, post, contract=True, debug=debug, threads=_threads(threads))
    rd = Reader(prog)
    res = [rd.scalar(b) for b in bodies]
    return res if batch else res[0]


# ---------------------------------------------------------------------------- conjugation
def _conj_scalar(c, complex_symbols=()):
    """Complex conjugate with every symbol real except those in complex_symbols; (T^a)_{ij} -> (T^a)_{ji}."""
    from .._parse import swap_i
    c = swap_i(c, -I)
    from .color import COLT
    ts = c.find(COLT(SR.wild(0), SR.wild(1), SR.wild(2)))
    if ts:
        c = c.subs({t: COLT(t.operands()[0], t.operands()[2], t.operands()[1]) for t in ts})
    for s in complex_symbols:
        c = c.subs({SR(s): SR(s).conjugate()})
    return c


def _index_counts(e, acc):
    """Add to acc the number of times each Lorentz index occurs in a monomial of e (max over sums)."""
    e = SR(e)
    op = _opname(e)
    if op is None:
        if T.is_summable(e):
            acc[str(e)] = acc.get(str(e), 0) + 1
        return acc
    if op in ('add_vararg', 'add'):
        best = {}
        for a in e.operands():
            for k, n in _index_counts(a, {}).items():
                best[k] = max(best.get(k, 0), n)
        for k, n in best.items():
            acc[k] = acc.get(k, 0) + n
        return acc
    if op == 'pow':
        b, x = e.operands()
        sub = _index_counts(b, {})
        n = int(x) if x.is_integer() and x > 0 else 1
        for k, m in sub.items():
            acc[k] = acc.get(k, 0) + n * m
        return acc
    for a in e.operands():
        _index_counts(a, acc)
    return acc


def _factor_indices(F, acc):
    best = {}
    for a, c in F.items:
        sub = _index_counts(c, {})
        if a[0] == 'i':
            sub[str(a[1])] = sub.get(str(a[1]), 0) + 1
        for k, n in sub.items():
            best[k] = max(best.get(k, 0), n)
    for k, n in best.items():
        acc[k] = acc.get(k, 0) + n
    return acc


def _term_index_counts(t):
    c, o, cl = t
    acc = _index_counts(c, {})
    for ch in ([o] if o is not None else []) + list(cl):
        for F in ch.factors:
            _factor_indices(F, acc)
    return acc


def _rename_factor(F, sub):
    by_name = {str(k): v for k, v in sub.items()}
    return Factor([((a[0], by_name.get(str(a[1]), a[1])) if a[0] == 'i' else a, c.subs(sub) if sub else c)
                   for a, c in F.items])


def _swap_vectors(F, sub):
    by_name = {str(k): v for k, v in sub.items()}
    return Factor([((a[0], by_name.get(str(a[1]), a[1])) if a[0] == 'p' else a, c) for a, c in F.items])


def _guard_indices(expr):
    """Raise if an index appears more than twice in one term (Einstein notation has no meaning then)."""
    if not isinstance(expr, DiracExpr) and not any(T.is_summable(x) for x in SR(expr).variables()):
        return                                       # no index at all: skip the walk over a big result
    terms = expr.terms if isinstance(expr, DiracExpr) else [(SR(expr), None, ())]
    for t in terms:
        many = [k for k, n in _term_index_counts(t).items() if n > 2]
        if many:
            raise DiracError("the index %s appears more than twice in one term; an index is either free (once) "
                             "or summed (twice).  Rename one pair." % many[0])


def dummy_indices(expr):
    """The Lorentz indices that occur twice in some term (summed), as Sage symbols."""
    out = set()
    terms = expr.terms if isinstance(expr, DiracExpr) else [(SR(expr), None, ())]
    for t in terms:
        out |= {k for k, n in _term_index_counts(t).items() if n >= 2}
    return sorted(SR.var(k) for k in out)


def conjugate(M, complex_symbols=()):
    r"""
    The complex conjugate of an amplitude built from spinor bilinears:
        [ubar(p2) G u(p1)]^* = ubar(p1) Gbar u(p2),   Gbar = gamma^0 G^dagger gamma^0
    (the order of the matrices is reversed, gamma5 -> -gamma5, i -> -i).  Summed Lorentz indices are
    renamed (mu -> mu_c) so that M * conjugate(M) has no index four times.  All symbols are taken real
    except those listed in complex_symbols.  Anything that is not a DiracExpr goes to Sage's conjugate.
    """
    if not isinstance(M, DiracExpr):
        from sage.all import conjugate as _sage_conjugate
        return _sage_conjugate(M)
    k = M.kind()
    if k not in ('scalar', None):
        raise DiracError("conjugate() is for amplitudes (spinor bilinears and numbers); for a Dirac matrix G "
                         "use dirac_bar(G) = gamma^0 G^dagger gamma^0")
    conj = lambda c: _conj_scalar(c, complex_symbols)
    taken = set(T._KIND)
    out = []
    for t in M.terms:
        c, o, cl = t
        dummies = [x for x, n in _term_index_counts(t).items() if n >= 2]
        sub = {}
        for x in dummies:
            new = T._fresh_index(x, taken)
            taken.add(str(new))
            sub[SR.var(x)] = new
        for name in T._POL:                                   # eps(k) <-> eps*(k)
            sub[SR.var(name)] = T.conjugate_vector(SR.var(name))
        chains = []
        for ch in cl:
            fs = tuple(_swap_vectors(_rename_factor(F, sub), sub).bar(conj) for F in reversed(ch.factors))
            chains.append(Chain(ch.right.dagger_bar(), fs, ch.left.dagger_bar()))
        cc = conj(c.subs(sub) if sub else c)
        out.append((cc, None, tuple(chains)))
    return DiracExpr(_combine(out))


def dirac_bar(G):
    """gamma^0 G^dagger gamma^0 for a product of Dirac matrices: reversed order, gamma5 -> -gamma5, i -> -i."""
    if not isinstance(G, DiracExpr) or G.kind() not in ('matrix', None):
        raise DiracError("dirac_bar needs a Dirac matrix (no spinors)")
    conj = lambda c: _conj_scalar(c)
    return DiracExpr([(conj(c), Chain(None, tuple(F.bar(conj) for F in reversed(o.factors)), None), ())
                      for c, o, cl in G.terms])


# ---------------------------------------------------------------------------- spin sums
def _projector(s):
    """sum over spins of w(p) wbar(p): p/ + m for u, p/ - m for v."""
    m = s.m if s.particle == 'u' else -s.m
    return Factor([(('p', s.p), SR(1)), (('1',), m)])


def _cycles(chains):
    """Join the bilinears of one term into closed spin-sum loops: [[factors of loop 1], ...]."""
    by_left = {}
    for i, ch in enumerate(chains):
        k = ch.left.key()
        if k in by_left:
            raise DiracError("spin_sum: the spinor %r appears twice; each spinor must appear once in M and once in M*" % ch.left)
        by_left[k] = i
    used, loops = set(), []
    for start in range(len(chains)):
        if start in used:
            continue
        seq, i = [], start
        while True:
            used.add(i)
            ch = chains[i]
            seq += list(ch.factors)
            seq.append(_projector(ch.right))
            j = by_left.get(ch.right.partner_key())
            if j is None:
                raise DiracError("spin_sum: %r has no partner %s(%s%s) to be summed with" %
                                 (ch.right, ch.right.partner_key()[0], ch.right.p,
                                  '' if ch.right.m.is_trivial_zero() else ', %s' % ch.right.m))
            if j == start:
                break
            if j in used:
                raise DiracError("spin_sum: the spinors do not close into loops")
            i = j
        loops.append(seq)
    return loops


def spin_sum(expr, dim=4, rules=None, contract=True, euclidean=None, gamma5_scheme="NDR", eps_in_d=False,
             debug=False, threads=None):
    r"""
    Sum a product of spinor bilinears over all spins, with
        sum_s u(p,s) ubar(p,s) = slash(p) + m,     sum_s v(p,s) vbar(p,s) = slash(p) - m,
    and let FORM do the traces (one trace per closed loop of spinors).  Usual use:
        M = ubar(p2, m) * gamma(mu) * u(p1, m) * ...
        spin_sum(M * conjugate(M), rules=...)
    Averages over initial spins (1/4 for two fermions) are not included.
    """
    dim = T.normalize_dim(dim)
    batch = isinstance(expr, (list, tuple))
    vr, rules = _vector_rules(rules)
    items = [_apply_vector_rules(e, vr) for e in (list(expr) if batch else [expr])]
    eu = T._euclid(euclidean)
    if eu:
        raise DiracError("spin sums are written for Minkowski spinors; use euclidean=False")
    prog = Program(dim, eu)
    texts, nlines = [], 0
    for e in items:
        if not isinstance(e, DiracExpr) or e.kind() not in ('scalar', None):
            raise DiracError("spin_sum needs a product of spinor bilinears, such as M * conjugate(M)")
        _guard_indices(e)
        parts = []
        for c, o, cl in e.terms:
            loops = [[_vr_factor(F, vr) for F in lp] for lp in _cycles(list(cl))]   # the projectors p/ + m too
            nlines = max(nlines, len(loops))
            parts += _line_texts(prog, c, loops, dim == 4, gamma5_scheme == "NDR-even")
        texts.append('+'.join(parts) if parts else '0')
    _check(prog, dim, eps_in_d, gamma5_scheme)
    vec, post = prog.rules(rules)
    bodies = _execute(prog, texts, nlines, vec, post, trace4=(dim == 4), contract=contract, debug=debug,
                      threads=_threads(threads))
    rd = Reader(prog)
    res = [rd.scalar(b) for b in bodies]
    return res if batch else res[0]


# ---------------------------------------------------------------------------- polarization sums
def polarization_sum(expr, e, gauge="feynman", n=None, dim=4, rules=None, debug=False, threads=None):
    r"""
    Sum a Sage expression (for example the output of spin_sum) over the polarizations of the vector
    e = polarization("e", k).  It must be linear in e and in its conjugate e_c.  The replacement of
    e^a e_c^b is
        gauge="feynman"    -g^{ab}                                   (photons: fine by the Ward identity)
        gauge="physical"   -g^{ab} + (k^a n^b + n^a k^b)/(k.n) - n^2 k^a k^b/(k.n)^2   (massless; needs n)
        massive vector (polarization(..., mass=M))   -g^{ab} + k^a k^b / M^2      (gauge is ignored)
    """
    dim = T.normalize_dim(dim)
    if isinstance(e, (list, tuple)):
        for x in e:
            expr = polarization_sum(expr, x, gauge, n, dim, rules, debug, threads)
        return expr
    hit = T._POL.get(str(e))
    if hit is None:
        raise DiracError("%s is not a polarization vector; make it with polarization(name, k)" % e)
    k, ec, mass = hit
    prog = Program(dim, False)
    text = prog.scalar(SR(expr))
    ev, ecv, kv = prog.name('vector', SR.var(str(e))), prog.name('vector', SR.var(ec)), prog.name('vector', k)
    if not mass.is_trivial_zero():
        rhs = '-d_(fsa,fsb) + %s(fsa)*%s(fsb)*%s' % (kv, kv, prog.scalar(1 / mass ** 2))
    elif gauge == "feynman":
        rhs = '-d_(fsa,fsb)'
    elif gauge == "physical":
        if n is None:
            raise DiracError('gauge="physical" needs a reference vector n (n.k != 0), e.g. n = another momentum')
        nv = prog.name('vector', SR(n))
        den = prog.opaque(1 / T.dot(k, n))
        rhs = ('-d_(fsa,fsb) + (%s(fsa)*%s(fsb) + %s(fsa)*%s(fsb))*%s - %s.%s*%s(fsa)*%s(fsb)*%s^2'
               % (kv, nv, nv, kv, den, nv, nv, kv, kv, den))
    else:
        raise ValueError('gauge must be "feynman" or "physical"')
    _check(prog, dim, False)
    vec, post = prog.rules(rules)
    src = [prog.declarations(), 'Indices fsa,fsb;', 'Tensors FSTA,FSTB;', 'Off statistics;',
           'Local FSE1 = %s;' % text, 'totensor %s,FSTA;' % ev, 'totensor %s,FSTB;' % ecv, '.sort',
           'id FSTA(fsa?)*FSTB(fsb?) = %s;' % rhs, 'contract;']
    if post:
        src += post + ['.sort'] + post
    src += ['if (count(FSTA,1,FSTB,1) > 0) exit "feynsage: polarization vector not linear";',
            'Format nospaces;', 'Print;', '.end', '']
    code = '\n'.join(src)
    if debug:
        print('--- FORM program ---\n' + code)
    from .backend import FormError
    try:
        out = backend().run(code, threads=_threads(threads))
    except FormError as err:
        if 'not linear' in str(err.output):
            raise DiracError("polarization_sum: the expression is not linear in %s and %s (each must appear once "
                             "per term, as in M * conjugate(M))" % (e, ec))
        raise
    if debug:
        print('--- FORM output ---\n' + out)
    return Reader(prog).scalar(expression_text(out, 'FSE1'))
