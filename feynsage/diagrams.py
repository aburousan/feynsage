r"""
Feynman diagrams from topologies and a model, and their amplitudes in textbook notation
(FeynArts' InsertFields and CreateFeynAmp, then FeynCalc's FCFAConvert, in one place).

    from feynsage.models import SM
    tops = topologies(0, 2, 2)
    diags = insert_fields(tops, ["e", "e~"], ["mu", "mu~"], SM())     # e- e+ -> mu- mu+
    diags.counts()           # [0, 4, 0, 0] per topology, as FeynArts prints
    diags.draw()             # the diagrams
    p1, p2, p3, p4 = momenta("p1 p2 p3 p4")
    M = diags.amplitude([p1, p2], [p3, p4])       # sum of all diagrams, M (not i M)

Field names are those of the model ('e', 'e~', 'mu', 'A', 'Z', 'W+', 'H', 'G0', 'g', ...; also
'e-', 'e+', 'mu-', 'mu+', 'gamma').  Incoming particles are listed as they are; outgoing ones too
(the conversion to all-incoming fields is done here).

Conventions (those of FeynArts' Lorentz.gen): vertex factors are the model couplings (with their i),
propagators i(k/ + m)/(k^2 - m^2) along the fermion flow, -i[g - (1 - xi) k k/(k^2 - xi M^2)]/(k^2 - M^2)
for vectors, i/(k^2 - xi M^2) for scalars and ghosts.  i M is the product; M = -i * (i M).  Relative
signs of diagrams with several fermion lines come from the order of the external fermions.
"""
import itertools
from sage.all import SR, I, Graph, Permutation, sqrt

from .topologies import Topology, topologies


# ---------------------------------------------------------------------------- diagrams
class Diagram:
    """A topology with a field on every line.  fields[k] = (name, a, b): the field `name` flows
    from end point a to end point b of line k (for an external line: from the leg into the vertex,
    with the all-incoming field)."""

    def __init__(self, top, fields, model, ext_in, ext_out):
        self.top, self.fields, self.model = top, fields, model
        self.ext_in, self.ext_out = list(ext_in), list(ext_out)

    @property
    def internal(self):
        return [(k, f) for k, f in enumerate(self.fields) if f[1][0] == 'v' and f[2][0] == 'v']

    def propagator_fields(self):
        return [f[0] for k, f in self.internal]

    def __repr__(self):
        inner = ', '.join('%s(v%d->v%d)' % (f[0], f[1][1], f[2][1]) for k, f in self.internal)
        return 'Diagram(%s -> %s: %s)' % (' '.join(self.ext_in), ' '.join(self.ext_out), inner or 'contact')

    def amplitude(self, p_in, p_out, **kw):
        return _amplitude(self, p_in, p_out, **kw)


class DiagramList(list):
    def __init__(self, items, tops, per_top):
        super().__init__(items)
        self.tops, self.per_top = tops, per_top

    def counts(self):
        """Number of insertions per topology (as FeynArts prints 'Top. n: k insertions')."""
        return list(self.per_top)

    def draw(self, ncols=4, size=2.2):
        from .diagram_plot import draw_diagrams
        return draw_diagrams(self, ncols=ncols, size=size)

    def amplitudes(self, p_in, p_out, **kw):
        """The amplitude M of every diagram (None for a diagram absent in unitary gauge)."""
        sh = kw.pop('shared', None) or _Shared()
        return [d.amplitude(p_in, p_out, prefix='d%d' % i, shared=sh, **kw) for i, d in enumerate(self)]

    def squared(self, p_in, p_out, **kw):
        """M * conj(M) summed over the spins of external fermions (spin_sum); polarization and
        colour sums are left to polarization_sum and color_factor."""
        from .qft.operations import spin_sum
        sh = _Shared()
        M = _sum_amplitudes([a for a in self.amplitudes(p_in, p_out, shared=sh, **kw) if a is not None])
        Mc = conjugate_amplitude(M, keep=sh.col.values())
        prod = M * Mc
        return prod if isinstance(prod, type(SR(1))) else spin_sum(prod)

    def loop_amplitude(self, p_in, p_out, project, kin, xi=1, nproc=None):
        """One-loop M of every diagram (see diagrams.loop_amplitude), diagrams in parallel."""
        return loop_amplitude(self, p_in, p_out, project, kin, xi=xi, nproc=nproc)

    def amplitude(self, p_in, p_out, **kw):
        """Sum of the amplitudes M of all diagrams."""
        return _sum_amplitudes([a for a in self.amplitudes(p_in, p_out, **kw) if a is not None])


def conjugate_amplitude(M, keep=()):
    """
    The complex conjugate of an amplitude: i -> -i, every polarization vector -> its conjugate and
    back, and every Lorentz or colour index renamed (they are summed inside M, so M * conj(M) must not
    reuse them).  `keep`: indices summed between M and conj(M) (the external colour indices).
    Amplitudes with spinors go through feynsage's own conjugate().
    """
    from .qft import tensors as T
    from .qft.operations import conjugate as dconj
    if not isinstance(M, type(SR(1))):
        return dconj(M)
    from ._parse import swap_i
    keep = set(str(k) for k in keep)
    sub = {}
    for v in M.variables():
        n = str(v)
        kind = T._KIND.get(n)
        if kind in ('index', 'adjoint', 'fundamental') and n not in keep:
            sub[v] = _one(T._declare(n + 'cc', kind))
        elif n in T._POL:
            sub[v] = T.conjugate_vector(v)
        else:
            for pn, (k, ec, m) in T._POL.items():
                if ec == n:
                    sub[v] = SR.var(pn)
    return swap_i(M, -I).subs(sub)


def _sum_amplitudes(amps):
    tot = amps[0]
    for a in amps[1:]:
        tot = tot + a
    return tot


def insert_fields(tops, fields_in, fields_out, model, exclude_fields=()):
    """
    Put fields of `model` on the lines of the topologies so that every vertex is a vertex of the
    model, like FeynArts' InsertFields[tops, in -> out] (at the level of particles).
    exclude_fields: names of fields that may not appear on internal lines (e.g. ["G0", "G+"]).
    Returns a DiagramList; .counts() gives the number of diagrams per topology.
    """
    if isinstance(tops, Topology):
        tops = [tops]
    fin = [model.aliases.get(f, f) for f in fields_in]
    fout = [model.aliases.get(f, f) for f in fields_out]
    for f in fin + fout:
        model.field(f)
    allin = fin + [model.anti(f) for f in fout]            # all incoming
    excl = set(model.aliases.get(f, f) for f in exclude_fields)
    excl |= set(model.anti(f) for f in excl)
    # candidate propagating fields: one per particle/antiparticle pair
    cands = []
    for name, p in sorted(model.particles.items()):
        if name in excl:
            continue
        if p.self_conjugate or name < p.anti:
            cands.append(name)
    out, per_top = [], []
    for top in tops:
        if top.E != len(allin):
            raise ValueError("topology has %d legs, the process %d" % (top.E, len(allin)))
        found = _insert_one(top, allin, model, cands)
        per_top.append(len(found))
        out += [Diagram(top, f, model, fin, fout) for f in found]
    return DiagramList(out, tops, per_top)


def _insert_one(top, allin, model, cands):
    lines = top.lines
    # vertex -> list of (line index, end) ; incoming field at the vertex depends on the line's field
    vlines = {}
    for k, (a, b) in enumerate(lines):
        for end in (a, b):
            if end[0] == 'v':
                vlines.setdefault(end[1], []).append(k)
    degree = {v: len(ls) for v, ls in vlines.items()}
    # admissible vertex field multisets by degree, for pruning
    by_deg = {}
    for vx in model.vertices:
        by_deg.setdefault(len(vx.fields), set()).add(tuple(sorted(vx.fields)))
    def possible(partial, deg):
        ps = sorted(partial)
        for s in by_deg.get(deg, ()):
            rest = list(s)
            ok = True
            for x in ps:
                if x in rest:
                    rest.remove(x)
                else:
                    ok = False
                    break
            if ok:
                return True
        return False
    assign = [None] * len(lines)
    # external lines are fixed: field flows from the leg into the vertex
    internal = []
    for k, (a, b) in enumerate(lines):
        if a[0] == 'x' or b[0] == 'x':
            leg, v = (a, b) if a[0] == 'x' else (b, a)
            assign[k] = (allin[leg[1] - 1], leg, v)
        else:
            internal.append(k)
    def incoming_at(v):
        fs = []
        for k in dict.fromkeys(vlines[v]):              # a self-loop line is listed twice: once
            if assign[k] is None:
                continue
            name, s, t = assign[k]
            a, b = lines[k]
            # a self-loop line enters the vertex once and leaves once
            if a == b:
                fs += [name, model.anti(name)]
                continue
            fs.append(name if t == ('v', v) else model.anti(name))
        return fs
    def vertex_ok(v):
        fs = incoming_at(v)
        if len(fs) == degree[v]:
            return bool(model.find(fs))
        return possible(fs, degree[v])
    # initial pruning with external lines only
    for v in vlines:
        if not vertex_ok(v):
            return []
    results, seen = [], set()
    def rec(i):
        if i == len(internal):
            key = _diagram_key(top, assign, model)
            if key not in seen:
                seen.add(key)
                results.append(list(assign))
            return
        k = internal[i]
        a, b = lines[k]
        for name in cands:
            p = model.particles[name]
            orients = [(a, b)] if p.self_conjugate or a == b else [(a, b), (b, a)]
            for (s, t) in orients:
                assign[k] = (name, s, t)
                if vertex_ok(a[1]) and vertex_ok(b[1]):
                    rec(i + 1)
                assign[k] = None
    rec(0)
    return results


def _diagram_key(top, assign, model):
    """Canonical form of a topology with fields: line k becomes a - e0 - e1 - b with e0 coloured by
    (field, tail) and e1 by (field, head); a self-conjugate field colours both ends the same."""
    G = Graph()
    ext = [('x', i) for i in range(1, top.E + 1)]
    verts = sorted({p for l in top.lines for p in l if p[0] == 'v'})
    G.add_vertices(ext + verts)
    colours = {}
    for k, (name, s, t) in enumerate(assign):
        p = model.particles[name]
        if not p.self_conjugate and name > p.anti:          # store the particle's direction
            name, s, t = p.anti, t, s
        e0, e1 = ('e', k, 0), ('e', k, 1)
        G.add_vertices([e0, e1])
        G.add_edges([(s, e0), (e0, e1), (e1, t)])
        tag0 = (name, 'sc') if p.self_conjugate else (name, 'tail')
        tag1 = (name, 'sc') if p.self_conjugate else (name, 'head')
        colours.setdefault(tag0, []).append(e0)
        colours.setdefault(tag1, []).append(e1)
    part = [[x] for x in ext] + ([verts] if verts else []) + [colours[c] for c in sorted(colours)]
    C = G.canonical_label(partition=part)
    return tuple(sorted(tuple(sorted(e[:2])) for e in C.edges(labels=False, sort=False))) + \
        tuple(len(colours[c]) for c in sorted(colours)) + tuple(sorted(colours))


# ---------------------------------------------------------------------------- amplitudes
class _Shared:
    """Symbols that must be the same in every diagram of a process: polarization vectors and colour
    indices of the external legs."""
    def __init__(self):
        self.pol, self.col = {}, {}


def _amplitude(diag, p_in, p_out, gauge=None, colour=True, prefix='d', shared=None, widths=None):
    """M for one tree diagram, in feynsage's textbook notation (u, vbar, gamma, slash, dot, ...).
    gauge: None (the model's xi), 'feynman' (xi = 1) or 'unitary' (xi -> oo: diagrams with Goldstone
    bosons or ghosts give None, vector propagators get k k / M^2).  widths: {field name: Gamma}
    puts q^2 - M^2 + i M Gamma in the propagators of those fields (Breit-Wigner)."""
    from .qft import objects as O, tensors as T
    from .qft.color import quark_colors, color_indices
    top, model = diag.top, diag.model
    shared = shared or _Shared()
    if top.loops != 0:
        raise NotImplementedError("amplitudes are built for tree diagrams; for loops use loop() and form_factor()")
    E, n_in = top.E, top.n_in
    pext = list(p_in) + list(p_out)
    if len(pext) != E:
        raise ValueError("give %d incoming and %d outgoing momenta" % (n_in, top.n_out))
    if gauge == 'unitary':
        xi = None
    elif gauge == 'feynman':
        xi = SR(1)
    else:
        xi = getattr(model, 'xi', SR(1))
    lines = top.lines
    if xi is None and any(model.particles[f[0]].kind in ('S', 'U') and model.particles[f[0]].xi is not None
                          for f in diag.fields):
        return None                                     # Goldstone bosons and ghosts are absent
    qin = [pext[i] if i < n_in else -pext[i] for i in range(E)]      # all-incoming momenta
    def behind(k, s):
        adj = {}
        for j, (a, b) in enumerate(lines):
            if j != k:
                adj.setdefault(a, []).append(b)
                adj.setdefault(b, []).append(a)
        seen, todo = set(), [s]
        while todo:
            x = todo.pop()
            if x not in seen:
                seen.add(x)
                todo += adj.get(x, [])
        return seen
    mom = {}
    for k, (name, s, t) in enumerate(diag.fields):
        if s[0] == 'x':
            mom[k] = qin[s[1] - 1]
        else:
            mom[k] = sum((qin[x[1] - 1] for x in behind(k, s) if x[0] == 'x'), 0 * pext[0])
    lidx, cidx = {}, {}
    def index(name):
        return _one(T._declare(name, 'index'))
    for k, (name, s, t) in enumerate(diag.fields):
        p = model.particles[name]
        external = s[0] == 'x'
        if p.kind == 'V':
            if external or (xi is not None and xi == 1):
                mu = index('%sl%d' % (prefix, k))
                lidx[(k, s)] = lidx[(k, t)] = mu
            else:
                lidx[(k, s)], lidx[(k, t)] = index('%sl%da' % (prefix, k)), index('%sl%db' % (prefix, k))
        if colour and p.colour is not None:
            if external:
                leg = s[1]
                if leg not in shared.col:
                    shared.col[leg] = (_one(quark_colors('i%d' % leg)) if p.colour in ('3', '3b') else _one(color_indices('a%d' % leg)))
                cidx[k] = shared.col[leg]
            else:
                cidx[k] = (_one(quark_colors('%si%d' % (prefix, k))) if p.colour in ('3', '3b') else _one(color_indices('%sa%d' % (prefix, k))))
    vfac = {}
    for v in range(top.V):
        legs = []
        for k, (name, s, t) in enumerate(diag.fields):
            if t == ('v', v):
                legs.append((name, mom[k], k, t))
            if s == ('v', v):
                legs.append((model.anti(name), -mom[k], k, s))
        vx = model.find([l[0] for l in legs])
        if not vx:
            raise ValueError("no vertex %s in the model" % [l[0] for l in legs])
        perm = _match(vx[0].fields, legs)
        vfac[v] = _vertex_value(vx[0], [legs[i] for i in perm], lidx, cidx, model)
    factor = SR(1)
    for k, (name, s, t) in enumerate(diag.fields):
        p = model.particles[name]
        if s[0] == 'x':
            if p.kind == 'V':
                leg = s[1]
                if leg not in shared.pol:
                    shared.pol[leg] = T.polarization('eps%d' % leg, pext[leg - 1], mass=p.mass)
                e = shared.pol[leg] if leg <= n_in else T.conjugate_vector(shared.pol[leg])
                factor *= T.comp(e, lidx[(k, s)])
            continue
        q = mom[k]
        q2 = T.dot(q, q)
        if widths and (name in widths or p.anti in widths):
            G = widths.get(name, widths.get(p.anti))
            q2 = q2 + I * p.mass * G                     # q^2 - M^2 + i M Gamma in every denominator
        if p.kind == 'V':
            M2 = p.mass ** 2
            mu, nu = lidx[(k, s)], lidx[(k, t)]
            if xi is None and M2.is_trivial_zero():
                factor *= (-I / q2) * T.metric(mu, nu)          # massless: no unitary gauge, use xi = 1
            elif xi is None:
                factor *= (-I / (q2 - M2)) * (T.metric(mu, nu) - T.comp(q, mu) * T.comp(q, nu) / M2)
            elif xi == 1:
                factor *= -I / (q2 - M2)
            else:
                factor *= (-I / (q2 - M2)) * (T.metric(mu, nu) - (1 - xi) * T.comp(q, mu) * T.comp(q, nu) / (q2 - xi * M2))
        elif p.kind in ('S', 'U'):
            m2 = p.mass ** 2 * (xi if p.xi is not None else 1)
            factor *= I / (q2 - m2)
    for v in range(top.V):
        if vfac[v][1] is None:
            factor *= vfac[v][0]
    chains, order = _fermion_chains(diag, mom, pext, vfac, model)
    if order:
        ranks = sorted(order)
        factor *= Permutation([ranks.index(x) + 1 for x in order]).sign()
    total = factor
    for ch in reversed(chains):
        total = ch * total
    return -I * total


def _one(x):
    return x[0] if isinstance(x, (list, tuple)) else x


def _match(model_fields, legs):
    """A permutation perm with legs[perm[i]] carrying model_fields[i]."""
    fs = [l[0] for l in legs]
    for perm in itertools.permutations(range(len(legs))):
        if all(fs[perm[i]] == model_fields[i] for i in range(len(legs))):
            return perm
    raise ValueError("vertex does not match")


def _vertex_value(vx, legs, lidx, cidx, model, colour_out=None):
    """(factor, dirac part or None) of a vertex; legs in model order: (field, k_in, line, end)."""
    from .qft import objects as O, tensors as T
    from .qft.color import T_color, f_color, color_delta
    kind, c = vx.kind, vx.coeffs
    idx = [lidx.get((l[2], l[3])) for l in legs]
    ks = [l[1] for l in legs]
    col = SR(1)
    if vx.colour and not all(l[2] in cidx for l in legs if model.particles[l[0]].colour):
        vx_colour = None
    else:
        vx_colour = vx.colour
    if vx_colour == 'delta':
        col = color_delta(cidx[legs[0][2]], cidx[legs[1][2]])
    elif vx_colour == 'T':
        col = T_color(cidx[legs[2][2]], cidx[legs[0][2]], cidx[legs[1][2]])
    elif vx_colour == 'f':
        col = f_color(cidx[legs[0][2]], cidx[legs[1][2]], cidx[legs[2][2]])
    if colour_out is not None:
        if vx_colour == 'ffff':
            # the colour of the four-gluon vertex does not factor out: each of its three Lorentz
            # structures gets a placeholder weight, replaced by its own colour factor after the integral
            g = T.metric
            m1, m2, m3, m4 = idx
            st = [g(m1, m2) * g(m3, m4), g(m1, m3) * g(m2, m4), g(m1, m4) * g(m3, m2)]
            a = [cidx[l[2]] for l in legs]
            cc = [_ff(a[0], a[2], a[1], a[3]) - _ff(a[0], a[3], a[2], a[1]),
                  _ff(a[0], a[1], a[2], a[3]) + _ff(a[0], a[3], a[2], a[1]),
                  -_ff(a[0], a[1], a[2], a[3]) - _ff(a[0], a[2], a[1], a[3])]
            n = sum(1 for x in colour_out if isinstance(x, tuple))
            ws = [SR.var('gluon4w%d_%d' % (n, k)) for k in range(3)]
            colour_out.append(('ffff', ws, cc))
            return (vx.factor * sum(w * x for w, x in zip(ws, st)), None, legs)
        colour_out.append(col)
        col = SR(1)
    if kind == 'FFV':
        gam = O.gamma(idx[2])
        dirac = gam * (c[0] * O.PL() + c[1] * O.PR())
        return (vx.factor * col, dirac, legs)
    if kind == 'FFS':
        return (vx.factor * col, c[0] * O.PL() + c[1] * O.PR(), legs)
    g, comp = T.metric, T.comp
    if kind == 'VVV':
        m1, m2, m3 = idx
        k1, k2, k3 = ks
        val = g(m1, m2) * comp(k2 - k1, m3) + g(m2, m3) * comp(k3 - k2, m1) + g(m3, m1) * comp(k1 - k3, m2)
        return (vx.factor * c[0] * col * val, None, legs)
    if kind == 'VVVV':
        m1, m2, m3, m4 = idx
        st = [g(m1, m2) * g(m3, m4), g(m1, m3) * g(m2, m4), g(m1, m4) * g(m3, m2)]
        if vx_colour == 'ffff':
            a = [cidx[l[2]] for l in legs]
            F = lambda w, x, y, z: _ff(w, x, y, z)
            cc = [F(a[0], a[2], a[1], a[3]) - F(a[0], a[3], a[2], a[1]),
                  F(a[0], a[1], a[2], a[3]) + F(a[0], a[3], a[2], a[1]),
                  -F(a[0], a[1], a[2], a[3]) - F(a[0], a[2], a[1], a[3])]
            return (vx.factor * sum(x * y for x, y in zip(cc, st)), None, legs)
        return (vx.factor * col * sum(x * y for x, y in zip(c, st)), None, legs)
    if kind in ('SSS', 'SSSS', 'SUU'):
        return (vx.factor * c[0] * col, None, legs)
    if kind == 'SSV':
        return (vx.factor * c[0] * col * comp(ks[0] - ks[1], idx[2]), None, legs)
    if kind == 'SVV':
        return (vx.factor * c[0] * col * g(idx[1], idx[2]), None, legs)
    if kind == 'SSVV':
        return (vx.factor * c[0] * col * g(idx[2], idx[3]), None, legs)
    if kind == 'UUV':
        return (vx.factor * col * (c[0] * comp(ks[0], idx[2]) + c[1] * comp(ks[1], idx[2])), None, legs)
    raise NotImplementedError(kind)


def _ff(a, b, c, d):
    """sum_e f^{a b e} f^{e c d} with a fresh summed adjoint index."""
    from .qft.color import color_indices, f_color
    e = _one(color_indices('ffs%d' % (abs(hash((str(a), str(b), str(c), str(d)))) % 10 ** 8)))
    return f_color(a, b, e) * f_color(e, c, d)


def _fermion_chains(diag, mom, pext, vfac, model):
    """The spinor chains [bar spinor] Gamma S Gamma ... [spinor] (against the fermion-number arrow)
    and the order of external fermion legs (end, start, end, start, ...) for the relative sign."""
    from .qft import objects as O, tensors as T
    top = diag.top
    n_in = top.n_in
    F = lambda name: model.particles[name].kind == 'F'
    flines = [k for k, f in enumerate(diag.fields) if F(f[0])]
    if not flines:
        return [], []
    # arrow direction of each fermion line: the field stored is the all-incoming one on external
    # lines, a particle or antiparticle on internal ones; the arrow follows the particle
    arrow = {}
    for k in flines:
        name, s, t = diag.fields[k]
        p = model.particles[name]
        is_particle = not name.endswith('~')
        arrow[k] = (s, t) if is_particle else (t, s)
    # at each vertex: the fermion line whose arrow enters and the one whose arrow leaves
    enter, leave = {}, {}
    for k in flines:
        a, b = arrow[k]
        if b[0] == 'v':
            enter.setdefault(b[1], []).append(k)
        if a[0] == 'v':
            leave.setdefault(a[1], []).append(k)
    def spinor(k, which):
        name, s, t = diag.fields[k]
        leg = s[1] if s[0] == 'x' else t[1]
        p = model.particles[name]
        m = p.mass
        q = pext[leg - 1]
        incoming = leg <= n_in
        # all-incoming field `name` on the leg
        if which == 'start':      # arrow into the diagram at this leg
            return O.u(q, m) if incoming else O.v(q, m)
        return O.vbar(q, m) if incoming else O.ubar(q, m)
    ends = [k for k in flines if diag.fields[k][1][0] == 'x' or diag.fields[k][2][0] == 'x']
    used, chains, order = set(), [], []
    for k0 in ends:
        a, b = arrow[k0]
        if b[0] != 'x' or k0 in used:
            continue                       # a chain starts where the arrow leaves the diagram (b = leg)
        used.add(k0)
        expr = spinor(k0, 'end')
        leg_end = b[1]
        k = k0
        while True:
            v = arrow[k][0][1]             # the vertex the arrow leaves from
            fac, dirac, legs = vfac[v]
            expr = expr * fac * dirac
            (kin,) = enter[v]
            if diag.fields[kin][1][0] == 'x' or diag.fields[kin][2][0] == 'x':
                used.add(kin)
                expr = expr * spinor(kin, 'start')
                leg_start = arrow[kin][0][1]
                break
            q = mom[kin] if arrow[kin] == (diag.fields[kin][1], diag.fields[kin][2]) else -mom[kin]
            m = model.particles[diag.fields[kin][0]].mass
            prop = I * (O.slash(q) + m) / (T.dot(q, q) - m ** 2)
            expr = expr * prop
            k = kin
        chains.append(expr)
        order += [leg_end, leg_start]
    return chains, order


# ---------------------------------------------------------------------------- one-loop amplitudes
def _line_momenta(top, fields, qin, lmom):
    """Momentum of every line along its field direction.  At one loop one line of the loop carries
    lmom; the others follow from momentum conservation (solved exactly)."""
    from .topologies import _components_without
    lines = top.lines
    mom, unknown = {}, {}
    internal = [k for k, (n, s, t) in enumerate(fields) if s[0] == 'v' and t[0] == 'v']
    for k, (n, s, t) in enumerate(fields):
        if s[0] == 'x':
            mom[k] = qin[s[1] - 1]
    chord = None
    if top.loops == 1:
        for k in internal:
            a, b = lines[k]
            if a == b or len(_components_without(lines, k)) == 1:
                chord = k
                break
        mom[chord] = lmom
    elif top.loops > 1:
        raise NotImplementedError("amplitudes with more than one loop")
    for k in internal:
        if k not in mom:
            unknown[k] = SR.var('fsmom%d' % k)
            mom[k] = unknown[k]
    eqs = []
    for v in range(top.V):
        tot = SR(0)
        for k, (n, s, t) in enumerate(fields):
            if t == ('v', v):
                tot += mom[k]
            if s == ('v', v):
                tot -= mom[k]
        eqs.append(tot == 0)
    if unknown:
        from sage.all import solve
        sol = solve(eqs, list(unknown.values()), solution_dict=True)
        if not sol:
            raise ValueError("momentum conservation has no solution")
        for k, x in unknown.items():
            mom[k] = sol[0][x].expand()
    return mom


def _symmetry_factor(top, fields, model):
    """1/(number of automorphisms of the diagram that keep the external legs fixed)."""
    G = Graph()
    ext = [('x', i) for i in range(1, top.E + 1)]
    verts = sorted({p for l in top.lines for p in l if p[0] == 'v'})
    G.add_vertices(ext + verts)
    colours = {}
    for k, (name, s, t) in enumerate(fields):
        p = model.particles[name]
        if not p.self_conjugate and name > p.anti:
            name, s, t = p.anti, t, s
        e0, e1 = ('e', k, 0), ('e', k, 1)
        G.add_vertices([e0, e1])
        G.add_edges([(s, e0), (e0, e1), (e1, t)])
        colours.setdefault((name, 'sc') if p.self_conjugate else (name, 'tail'), []).append(e0)
        colours.setdefault((name, 'sc') if p.self_conjugate else (name, 'head'), []).append(e1)
    part = [[x] for x in ext] + ([verts] if verts else []) + [colours[c] for c in sorted(colours)]
    return SR(1) / G.automorphism_group(partition=part).order()


def _loop_integrand(diag, p_in, p_out, lmom, xi=SR(1), colour=True, prefix='d', shared=None, spinors=False,
                    conserve=False):
    """
    (numerator, propagators, external indices) of a one-loop (or tree) diagram:
        i M = int d^d l/(2 pi)^d  numerator / prod (q_i^2 - m_i^2)
    numerator: a Sage expression (Dirac traces of closed fermion loops already done, in d dimensions)
    with a free Lorentz index on every external vector (no polarization vectors);
    propagators: [(q_i, m_i)] of the lines that carry the loop momentum.
    """
    from .qft import objects as O, tensors as T
    from .qft.operations import dirac_trace
    from .qft.color import quark_colors, color_indices
    top, model = diag.top, diag.model
    shared = shared or _Shared()
    E, n_in = top.E, top.n_in
    pext = list(p_in) + list(p_out)
    qin = [pext[i] if i < n_in else -pext[i] for i in range(E)]
    if conserve:                        # the last momentum from the others, for the internal lines only
        qin[E - 1] = -sum(qin[:E - 1], 0 * pext[0])
    mom = _line_momenta(top, diag.fields, qin, lmom)
    lidx, cidx = {}, {}
    ext_index = {}
    def index(name):
        return _one(T._declare(name, 'index'))
    for k, (name, s, t) in enumerate(diag.fields):
        p = model.particles[name]
        external = s[0] == 'x'
        if p.kind == 'V':
            if external or xi == 1:
                mu = index('%sl%d' % (prefix, k) if not external else 'mu%d' % s[1])
                lidx[(k, 'S')] = lidx[(k, 'T')] = mu
                if external:
                    ext_index[s[1]] = mu
            else:
                # the two ends of the line by tail S / head T (a line from a vertex to itself has both there)
                lidx[(k, 'S')], lidx[(k, 'T')] = index('%sl%da' % (prefix, k)), index('%sl%db' % (prefix, k))
        if colour and p.colour is not None:
            if external:
                leg = s[1]
                if leg not in shared.col:
                    shared.col[leg] = (_one(quark_colors('i%d' % leg)) if p.colour in ('3', '3b') else _one(color_indices('a%d' % leg)))
                cidx[k] = shared.col[leg]
            else:
                cidx[k] = (_one(quark_colors('%si%d' % (prefix, k))) if p.colour in ('3', '3b') else _one(color_indices('%sa%d' % (prefix, k))))
    vfac, colours = {}, []
    for v in range(top.V):
        legs = []
        for k, (name, s, t) in enumerate(diag.fields):
            if t == ('v', v):
                legs.append((name, mom[k], k, 'T'))
            if s == ('v', v):
                legs.append((model.anti(name), -mom[k], k, 'S'))
        vx = model.find([l[0] for l in legs])
        if not vx:
            raise ValueError("no vertex %s in the model" % [l[0] for l in legs])
        vfac[v] = _vertex_value(vx[0], [legs[i] for i in _match(vx[0].fields, legs)], lidx, cidx, model,
                                colour_out=colours)
    l_free = lambda q: not SR(q).has(lmom)
    factor, line_opts, loop_k = SR(1), [], []   # line_opts: per loop line, a list of (numerator, [masses])
    for k, (name, s, t) in enumerate(diag.fields):
        if s[0] == 'x':
            continue
        p = model.particles[name]
        q = mom[k]
        q2 = T.dot(q, q)
        if p.kind == 'V':
            M2 = p.mass ** 2
            mu, nu = lidx[(k, 'S')], lidx[(k, 'T')]
            kk = T.comp(q, mu) * T.comp(q, nu)
            if xi == 1:
                opts = [(-I, [p.mass])]
            elif l_free(q):
                opts = [(-I * (T.metric(mu, nu) - (1 - xi) * kk / (q2 - xi * M2)), [p.mass])]
            elif M2.is_trivial_zero():
                # -i [g - (1 - xi) k k/k^2]/k^2: the second piece has the propagator twice
                opts = [(-I * T.metric(mu, nu), [0]), (I * (1 - xi) * kk, [0, 0])]
            else:
                # partial fractions: -i[g - k k/M^2]/(k^2 - M^2) - i (k k/M^2)/(k^2 - xi M^2)
                opts = [(-I * (T.metric(mu, nu) - kk / M2), [p.mass]), (-I * kk / M2, [sqrt(xi) * p.mass])]
        elif p.kind in ('S', 'U'):
            opts = [(I, [p.mass * (sqrt(xi) if p.xi is not None else 1)])]
        else:
            opts = [(I, [p.mass])]                  # fermion: i (k/ + m) goes into the Dirac string
        if l_free(q):
            num, (mass,) = opts[0]
            factor *= num / (q2 - mass ** 2)
        else:
            line_opts.append((q, opts))
            loop_k.append(k)
    for v in range(top.V):
        if vfac[v][1] is None:
            factor *= vfac[v][0]
    # fermion lines: closed loops become traces; open lines are spinor chains (spinors=True, any
    # number of them, with the sign of the external fermion order) or one bare Dirac string
    F = [k for k, f in enumerate(diag.fields) if model.particles[f[0]].kind == 'F']
    open_ends = [k for k in F if diag.fields[k][1][0] == 'x']
    if spinors and open_ends:
        chains, order, used = _loop_fermion_chains(diag, mom, pext, vfac, model)
        if order:
            ranks = sorted(order)
            factor *= Permutation([ranks.index(x) + 1 for x in order]).sign()
        rest = [k for k in F if k not in used]
        if rest:
            factor *= _closed_fermion_loop(diag, rest, mom, vfac, model)
        for ch in reversed(chains):
            factor = ch * factor
    elif open_ends:
        # one open fermion line: the Dirac string X of  ubar/vbar(end) X u/v(start)  (spinors left out,
        # the projector of loop_amplitude traces it, e.g. with vertex_projector)
        if len(open_ends) != 2:
            raise NotImplementedError("one-loop amplitudes with more than one open fermion line")
        arrow = {}
        for k in F:
            name, s0, t0 = diag.fields[k]
            arrow[k] = (s0, t0) if not name.endswith('~') else (t0, s0)
        enter = {}
        for k in F:
            if arrow[k][1][0] == 'v':
                enter.setdefault(arrow[k][1][1], []).append(k)
        k0 = [k for k in open_ends if arrow[k][1][0] == 'x'][0]       # the arrow leaves the diagram here
        X, k = None, k0
        while True:
            v = arrow[k][0][1]
            fac, dirac, legs = vfac[v]
            X = fac * dirac if X is None else X * (fac * dirac)
            (kin_,) = enter[v]
            if arrow[kin_][0][0] == 'x':
                break
            q, m = mom[kin_], model.particles[diag.fields[kin_][0]].mass
            X = X * (O.slash(q) + m)
            k = kin_
        factor = factor * X
    elif F:
        factor *= _closed_fermion_loop(diag, F, mom, vfac, model)
    loop_lines = [k for k, f in enumerate(diag.fields) if f[1][0] == 'v' and SR(mom[k]).has(lmom)]
    if loop_lines and all(model.particles[diag.fields[k][0]].kind == 'U' for k in loop_lines):
        factor *= -1                     # a closed ghost loop (an open ghost line through the loop has no sign)
    factor *= _symmetry_factor(top, diag.fields, model)
    ext = [ext_index.get(i) for i in range(1, E + 1)]
    from sage.all import prod
    plain = [c for c in colours if not isinstance(c, tuple)]
    marks = [c for c in colours if isinstance(c, tuple)]
    colour_factor_ = prod(plain) if plain else SR(1)
    if marks:
        colour_factor_ = (colour_factor_, marks)
    # one term for every choice of the propagator pieces of the loop lines
    terms = []
    # the loop lines in their order around the loop (boxes: the IR classification needs it)
    if len(loop_k) > 2:
        ends = {k: (diag.fields[k][1], diag.fields[k][2]) for k in loop_k}
        order, at = [loop_k[0]], ends[loop_k[0]][1]
        while len(order) < len(loop_k):
            nxt = [k for k in loop_k if k not in order and at in ends[k]]
            if not nxt:
                break
            k = nxt[0]
            at = ends[k][1] if ends[k][0] == at else ends[k][0]
            order.append(k)
        if len(order) == len(loop_k):
            line_opts = [line_opts[loop_k.index(k)] for k in order]
    for choice in itertools.product(*[opts for q, opts in line_opts]):
        num = factor
        props = []
        for (q, opts), (n, masses) in zip(line_opts, choice):
            num = num * n
            props += [(q, m_) for m_ in masses]
        terms.append((num, props))
    return terms, ext, colour_factor_


def _closed_fermion_loop(diag, F, mom, vfac, model):
    """-Tr[...] of the closed fermion loop made of the internal fermion lines F (in d dimensions)."""
    from .qft import objects as O
    from .qft.operations import dirac_trace
    arrow = {k: (diag.fields[k][1], diag.fields[k][2]) for k in F}     # internal: particle direction
    enter = {}
    for k in F:
        enter.setdefault(arrow[k][1][1], []).append(k)
    v0 = arrow[F[0]][1][1]
    v, chain = v0, None
    while True:
        fac, dirac, legs = vfac[v]
        piece = fac * dirac
        chain = piece if chain is None else chain * piece
        (kin,) = [k for k in enter[v] if k in F]
        q, m = mom[kin], model.particles[diag.fields[kin][0]].mass
        chain = chain * (O.slash(q) + m)
        v = arrow[kin][0][1]
        if v == v0:
            break
    return -dirac_trace(chain, dim='d')


def _loop_fermion_chains(diag, mom, pext, vfac, model):
    """The open fermion lines of a loop diagram as spinor chains [bar spinor] G (k/ + m) G ... [spinor]:
    the factor i and the denominator of each propagator are kept with the line (loop or tree part).
    Returns (chains, order of the external legs for the sign, the fermion lines used)."""
    from .qft import objects as O
    top = diag.top
    n_in = top.n_in
    flines = [k for k, f in enumerate(diag.fields) if model.particles[f[0]].kind == 'F']
    arrow = {}
    for k in flines:
        name, s0, t0 = diag.fields[k]
        arrow[k] = (s0, t0) if not name.endswith('~') else (t0, s0)
    enter = {}
    for k in flines:
        a, b = arrow[k]
        if b[0] == 'v':
            enter.setdefault(b[1], []).append(k)

    def spinor(k, which):
        name, s0, t0 = diag.fields[k]
        leg = s0[1] if s0[0] == 'x' else t0[1]
        m = model.particles[name].mass
        q = pext[leg - 1]
        incoming = leg <= n_in
        if which == 'start':
            return O.u(q, m) if incoming else O.v(q, m)
        return O.vbar(q, m) if incoming else O.ubar(q, m)
    ends = [k for k in flines if diag.fields[k][1][0] == 'x' or diag.fields[k][2][0] == 'x']
    used, chains, order = set(), [], []
    for k0 in ends:
        a, b = arrow[k0]
        if b[0] != 'x' or k0 in used:
            continue
        used.add(k0)
        expr = spinor(k0, 'end')
        leg_end = b[1]
        k = k0
        while True:
            v = arrow[k][0][1]
            fac, dirac, legs = vfac[v]
            expr = expr * (fac * dirac)
            (kin,) = enter[v]
            used.add(kin)
            if diag.fields[kin][1][0] == 'x' or diag.fields[kin][2][0] == 'x':
                expr = expr * spinor(kin, 'start')
                leg_start = arrow[kin][0][1]
                break
            q = mom[kin] if arrow[kin] == (diag.fields[kin][1], diag.fields[kin][2]) else -mom[kin]
            expr = expr * (O.slash(q) + model.particles[diag.fields[kin][0]].mass)
            k = kin
        chains.append(expr)
        order += [leg_end, leg_start]
    return chains, order, used


def _to_loop_string(expr, lname):
    """A contracted Sage expression in dot products -> the numerator string of loop()."""
    import re
    s = str(SR(expr))
    s = re.sub(r'Dot\((\w+), (\w+)\)', lambda m: ('%s^2' % m.group(1)) if m.group(1) == m.group(2)
               else '%s.%s' % (m.group(1), m.group(2)), s)
    s = re.sub(r'\bd\b', '(4 - 2*eps)', s)
    return s


def loop_amplitude(diags, p_in, p_out, project, kin, xi=1, lname='l', nproc=None):
    """
    M (not i M) of every one-loop diagram, projected to a scalar and reduced to A0, B0, C0, D0:
        project(numerator, external_indices) -> a Sage expression without free indices
    (for example the contraction with a form-factor projector), kin: the scalar products of the
    external momenta for loop() ({"p1^2": 0, "p1.p2": "MH^2/2"}).  Returns the list of results (Sage
    expressions with 1/eps poles); finite_part(sum) is the finite amplitude when the poles cancel.
    Closed fermion loops are traced in d dimensions (NDR); Feynman-'t Hooft gauge by default.
    The diagrams are projected and reduced in parallel (nproc: see feynsage.parallel; 1 = serial).
    """
    from .qft import tensors as T
    from .qft.operations import contract
    from .pv import loop
    lmom = _one(T._declare(lname, 'vector'))
    sh = _Shared()
    # the external kinematics as substitutions in the numerator (so no p.q is left in a denominator)
    import re
    from ._parse import sr as _sr
    kin_sr = {}
    for key, val in (kin or {}).items():
        m2 = re.fullmatch(r'\s*(\w+)\s*\^\s*2\s*', key)
        m1 = re.fullmatch(r'\s*(\w+)\s*\.\s*(\w+)\s*', key)
        a, b = (m2.group(1), m2.group(1)) if m2 else (m1.group(1), m1.group(2))
        va, vb = _one(T._declare(a, 'vector')), _one(T._declare(b, 'vector'))
        kin_sr[T.dot(va, vb)] = _sr(val) if isinstance(val, str) else SR(val)
    # the integrands one after the other (they share the index names), then the expensive part,
    # projection, traces and loop reduction, diagram by diagram in parallel
    from .qft.color import color_factor
    from .parallel import pmap
    jobs = []
    for i, dg in enumerate(diags):
        terms, ext, col = _loop_integrand(dg, p_in, p_out, lmom, xi=SR(xi), prefix='d%d' % i, shared=sh)
        # the colour factor does not depend on l: summed here, multiplied back after the integral.
        # With four-gluon vertices: one colour factor for each choice of their Lorentz structures.
        if isinstance(col, tuple):
            plain, marks = col
            combos = []
            for ks in itertools.product(range(3), repeat=len(marks)):
                c = plain
                pick = {}
                for (tag, ws, cc), k in zip(marks, ks):
                    c = c * cc[k]
                    pick.update({w: (1 if j == k else 0) for j, w in enumerate(ws)})
                combos.append((pick, color_factor(c, N=3)))
            col = combos
        else:
            col = color_factor(col, N=3) if not SR(col).is_numeric() else col
        jobs.append((terms, ext, col))

    def reduce_one(job):
        terms, ext, col = job
        total = SR(0)
        for num, props in terms:
            num = SR(contract(project(num, ext), dim='d')).subs(kin_sr).expand()
            if not props:
                total += -I * num
                continue
            pr, back = [], {}
            for q, m in props:
                qs = SR(q)
                if qs.coefficient(lmom) == -1:
                    qs = -qs
                m = SR(m)
                if not _rational_poly(m):            # e.g. sqrt(xi) M_W: a placeholder symbol in loop()
                    if m not in back.values():
                        back[SR.var('fsmass%d' % len(back))] = m
                    m = [k for k, v in back.items() if v is m or (v - m).is_trivial_zero()][0]
                pr.append([str(qs), str(m)])
            res = loop(_to_loop_string(num, lname), *pr, kin=kin)
            val = res.coefficients().get('1', SR(0)) if hasattr(res, 'coefficients') else SR(res)
            if back:
                val = SR(val).subs(back)
            total += val / (16 * SR.pi() ** 2)          # i M = i/(16 pi^2) * loop(...)  ->  M
        if isinstance(col, list):
            return sum((SR(total).subs(pick) * c for pick, c in col), SR(0))
        return col * total
    return pmap(reduce_one, jobs, nproc=nproc)


def _rational_poly(m):
    """True if m is a polynomial in its symbols with rational coefficients (what loop() accepts)."""
    try:
        from sage.all import QQ
        m.polynomial(QQ)
        return True
    except (TypeError, ValueError, AttributeError):
        return False


def vertex_projector(p1, p2, m, mu, part="F2", dim="d"):
    r"""
    The Dirac matrix P with  Tr[(p1/ + m) P (p2/ + m) X] = F  for an on-shell vertex  ubar(p2) X^mu u(p1)
    (both masses m), X^mu = A gamma^mu + B (p1 + p2)^mu + C q^mu, q = p2 - p1:
        part="F2": the Pauli form factor F2 = -2 m B  (Gordon: i sigma^{mu nu} q_nu/(2m) = gamma^mu - (p1 + p2)^mu/(2m)),
        part="F1": F1 = A + 2 m B,  part="A", "B", "C": the coefficients themselves.
    Built exactly in `dim` dimensions; parity-odd (gamma5) structures do not contribute.  Use the same
    Lorentz index mu as the free index of X (the external photon).
    """
    from .qft import objects as O, tensors as T
    from .qft.operations import dirac_trace, contract
    from sage.all import matrix, vector
    S = [O.gamma(mu), T.comp(p1 + p2, mu), T.comp(p2 - p1, mu)]
    kin = {T.dot(p1, p1): m ** 2, T.dot(p2, p2): m ** 2}
    M = matrix(SR, 3, 3)
    for i in range(3):
        for j in range(3):
            M[i, j] = SR(contract(dirac_trace((O.slash(p1) + m) * S[i] * (O.slash(p2) + m) * S[j], dim=dim),
                                  dim=dim)).subs(kin)
    target = {"A": [1, 0, 0], "B": [0, 1, 0], "C": [0, 0, 1],
              "F2": [0, -2 * m, 0], "F1": [1, 2 * m, 0]}[part]
    c = (M.inverse() * vector(SR, target)).apply_map(lambda e: SR(e).simplify_rational())
    return c[0] * S[0] + c[1] * S[1] + c[2] * S[2]
