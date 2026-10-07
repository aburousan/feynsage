r"""
Topologies of Feynman diagrams, as FeynArts' CreateTopologies makes them.

    tops = topologies(0, 2, 2)            # tree level, 2 -> 2: 4 topologies
    tops = topologies(1, 2, 2, exclude=("tadpoles", "wf"))

A topology has E external legs, numbered 1..E (the first n_in are incoming), and internal vertices
of degree 3 or 4 (`degrees`).  Lines are pairs of end points; an end point is ('x', i) for external
leg i or ('v', j) for internal vertex j.  Several lines may join the same two vertices and a line may
start and end at the same vertex (a tadpole loop at a quartic vertex).

How they are made.  Tree topologies: start from three legs on one vertex and add the legs one at a
time, either on a new vertex that splits an existing line or on an existing vertex of degree 3.  This
makes every tree with labelled legs exactly once.  L loops: make the trees with E + 2L legs and join
the extra legs in pairs (cutting one line of every loop of a graph gives such a tree).  Copies are
removed by a canonical labelling of the graph with every external leg kept in place.
"""
from itertools import combinations
from sage.all import Graph

DEFAULT_DEGREES = (3, 4)


class Topology:
    """External legs 1..E (first n_in incoming), internal vertices 0..V-1, lines [(a, b)]."""

    def __init__(self, n_in, n_out, lines):
        self.n_in, self.n_out = n_in, n_out
        self.E = n_in + n_out
        # renumber internal vertices 0..V-1 in order of appearance
        ren, out = {}, []
        for a, b in lines:
            pair = []
            for p in (a, b):
                if p[0] == 'v':
                    if p[1] not in ren:
                        ren[p[1]] = len(ren)
                    p = ('v', ren[p[1]])
                pair.append(p)
            out.append(tuple(pair))
        self.lines = out
        self.V = len(ren)

    @property
    def loops(self):
        return len(self.internal_lines()) - self.V + 1

    def external_lines(self):
        return [l for l in self.lines if l[0][0] == 'x' or l[1][0] == 'x']

    def internal_lines(self):
        return [l for l in self.lines if l[0][0] == 'v' and l[1][0] == 'v']

    def degree(self, v):
        return sum((a == ('v', v)) + (b == ('v', v)) for a, b in self.lines)

    def key(self):
        return _canonical(self.lines, self.E)

    def __repr__(self):
        def name(p):
            return str(p[1]) if p[0] == 'x' else 'v%d' % p[1]
        ext = ', '.join('%s-%s' % (name(a), name(b)) for a, b in self.external_lines())
        inn = ', '.join('%s-%s' % (name(a), name(b)) for a, b in self.internal_lines())
        return 'Topology(%d -> %d, %d loop%s: %s | %s)' % (self.n_in, self.n_out, self.loops,
                                                         '' if self.loops == 1 else 's', ext, inn)


def _canonical(lines, E):
    """A hashable form, the same for graphs that differ only by renaming internal vertices."""
    G = Graph(multiedges=False, loops=False)
    ext = [('x', i) for i in range(1, E + 1)]
    verts = sorted({p for l in lines for p in l if p[0] == 'v'})
    G.add_vertices(ext + verts)
    enodes = []
    for k, (a, b) in enumerate(lines):
        e0, e1 = ('e', k, 0), ('e', k, 1)            # every line becomes a path a - e0 - e1 - b
        G.add_vertices([e0, e1])
        G.add_edges([(a, e0), (e0, e1), (e1, b)])
        enodes += [e0, e1]
    part = [[x] for x in ext] + [verts] + [enodes]
    part = [c for c in part if c]
    C = G.canonical_label(partition=part)
    return tuple(sorted(tuple(sorted(e[:2])) for e in C.edges(labels=False, sort=False)))


def _trees(E, degrees):
    """All trees with labelled legs 1..E and internal degrees in `degrees` (each exactly once)."""
    if E < 3:
        raise ValueError("trees need at least three legs")
    start = [(('x', 1), ('v', 0)), (('x', 2), ('v', 0)), (('x', 3), ('v', 0))]
    level = [(start, 1)]                                  # (lines, next vertex number)
    for leg in range(4, E + 1):
        new = []
        for lines, nv in level:
            for k, (a, b) in enumerate(lines):            # split line k by a new cubic vertex
                if 3 not in degrees:
                    break
                w = ('v', nv)
                nl = lines[:k] + lines[k + 1:] + [(a, w), (w, b), (('x', leg), w)]
                new.append((nl, nv + 1))
            if 4 in degrees:                               # put the leg on a vertex of degree 3
                verts = sorted({p for l in lines for p in l if p[0] == 'v'})
                for v in verts:
                    if sum((a == v) + (b == v) for a, b in lines) == 3:
                        new.append((lines + [(('x', leg), v)], nv))
        level = new
    out = []
    for lines, _ in level:
        if all(sum((a == v) + (b == v) for a, b in lines) in degrees
               for v in {p for l in lines for p in l if p[0] == 'v'}):
            out.append(lines)
    return out


def _join(lines, a, b):
    """Join external legs a and b of a graph into one internal line."""
    na = [q for (p, q) in lines if p == ('x', a)] + [p for (p, q) in lines if q == ('x', a)]
    nb = [q for (p, q) in lines if p == ('x', b)] + [p for (p, q) in lines if q == ('x', b)]
    rest = [l for l in lines if ('x', a) not in l and ('x', b) not in l]
    return rest + [(na[0], nb[0])]


def _components_without(lines, drop):
    """Connected components (as sets of end points) of the graph without line number `drop`."""
    adj = {}
    for k, (a, b) in enumerate(lines):
        adj.setdefault(a, set())
        adj.setdefault(b, set())
        if k == drop:
            continue
        adj[a].add(b)
        adj[b].add(a)
    seen, comps = set(), []
    for s in adj:
        if s in seen:
            continue
        comp, todo = set(), [s]
        while todo:
            x = todo.pop()
            if x in comp:
                continue
            comp.add(x)
            todo += list(adj[x] - comp)
        seen |= comp
        comps.append(comp)
    return comps


def is_tadpole(top):
    """A part with no external leg hangs on a single line (FeynArts' Tadpoles)."""
    for k, (a, b) in enumerate(top.lines):
        if a[0] == 'x' or b[0] == 'x':
            continue
        comps = _components_without(top.lines, k)
        if len(comps) == 2 and any(not any(p[0] == 'x' for p in c) for c in comps):
            return True
    return False


def is_wf(top):
    """A self-energy on an external leg: one line cuts off exactly one external leg together with
    internal vertices (FeynArts' WFCorrections)."""
    for k, (a, b) in enumerate(top.lines):
        if a[0] == 'x' or b[0] == 'x':
            continue
        comps = _components_without(top.lines, k)
        if len(comps) == 2 and any(sum(p[0] == 'x' for p in c) == 1 for c in comps):
            return True
    return False


EXCLUDE = {'tadpoles': is_tadpole, 'wf': is_wf}


def topologies(loops, n_in, n_out, degrees=DEFAULT_DEGREES, exclude=()):
    """
    All topologies with `loops` loops for n_in -> n_out external legs and internal vertices of the
    given degrees, like FeynArts' CreateTopologies[loops, n_in -> n_out].  exclude: any of
    "tadpoles" (a part without external legs hanging on one line) and "wf" (self-energy insertions
    on external legs), like ExcludeTopologies -> {Tadpoles, WFCorrections}.
    """
    E = n_in + n_out
    seen, out = set(), []
    for tree in _trees(E + 2 * loops, degrees):
        lines = tree
        for j in range(loops):
            lines = _join(lines, E + 2 * j + 1, E + 2 * j + 2)
        t = Topology(n_in, n_out, lines)
        k = t.key()
        if k in seen:
            continue
        seen.add(k)
        if any(EXCLUDE[x](t) for x in exclude):
            continue
        out.append(t)
    # order like FeynArts: fewer internal vertices of high degree first is not essential; sort by
    # number of internal vertices so that contact topologies come first
    out.sort(key=_order_key)
    return out


def _order_key(t):
    """FeynArts' order at tree level: contact vertices first, then the channels by the legs that sit
    with leg 1 (s: {1, 2}, t: {1, 3}, u: {1, 4} for 2 -> 2)."""
    sides = []
    for k, (a, b) in enumerate(t.lines):
        if a[0] == 'x' or b[0] == 'x':
            continue
        comps = _components_without(t.lines, k)
        if len(comps) == 2:
            c = [c for c in comps if ('x', 1) in c][0]
            sides.append(tuple(sorted(p[1] for p in c if p[0] == 'x')))
    return (t.V, sorted(-t.degree(v) for v in range(t.V)), sorted(sides))
