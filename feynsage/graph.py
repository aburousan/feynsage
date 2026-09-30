r"""
Feynman graphs and their graph polynomials, read off with graph theory.

A FeynmanGraph has internal lines  (u, v, mass2)  between vertices, and external
momenta attached to vertices,  {vertex: momentum dict}  (incoming momenta; their sum
must vanish).  Multiple lines between the same pair of vertices are allowed.

- spanning_trees() / two_forests() enumerate edge subsets directly from the
  definitions (connected, loop-free, touching every vertex).
- U() and F() are then read off exactly as in the lecture:
    U   = sum over spanning trees of the product of x_e of the removed lines,
    F_0 = sum over 2-forests of (product of x_e of the removed lines) * P^2,
          P = total external momentum entering one of the two trees,
    F   = F_0 + U sum_e x_e m_e^2        (Euclidean signs).
- U_kirchhoff() computes U independently from the matrix-tree theorem.
- family() assigns loop momenta through a cycle basis and returns an
  IntegralFamily, so that the matrix method can be compared with the tree rules.
"""
from itertools import combinations
from sage.all import Graph, PolynomialRing, matrix, QQ, prod

from .momenta import add


class FeynmanGraph:

    def __init__(self, lines, externals, kinematics):
        self.lines = [(u, v, m2) for (u, v, m2) in lines]
        self.ext = {k: dict(v) for k, v in externals.items()}
        self.kin = kinematics
        verts = set()
        for u, v, _ in self.lines:
            verts.update((u, v))
        verts.update(self.ext)
        self.vertices = sorted(verts)
        self.N = len(self.lines)
        self.V = len(self.vertices)
        self.L = self.N - self.V + 1
        tot = {}
        for q in self.ext.values():
            tot = add(tot, q)
        if tot:
            raise ValueError("external momenta do not add up to zero: %s" % tot)
        self.R = PolynomialRing(self.kin.R, ['x%d' % (i + 1) for i in range(self.N)])
        self.x = self.R.gens()

    # ------------------------------------------------------------------ graph theory
    def _graph(self, keep):
        G = Graph(multiedges=True, loops=True)
        G.add_vertices(self.vertices)
        for i in keep:
            u, v, _ = self.lines[i]
            G.add_edge(u, v, i)
        return G

    def _is_forest(self, G):
        # a multigraph is a forest iff edges = vertices - components (no cycles, no multi-edges used twice)
        return G.size() == G.order() - len(G.connected_components())

    def spanning_trees(self):
        """Removed-line sets of all spanning trees (remove exactly L lines)."""
        out = []
        for rem in combinations(range(self.N), self.L):
            keep = [i for i in range(self.N) if i not in rem]
            G = self._graph(keep)
            if G.is_connected() and self._is_forest(G):
                out.append(rem)
        return out

    def two_forests(self):
        """(removed lines, one of the two components) for all spanning 2-forests."""
        out = []
        for rem in combinations(range(self.N), self.L + 1):
            keep = [i for i in range(self.N) if i not in rem]
            G = self._graph(keep)
            comps = G.connected_components()
            if len(comps) == 2 and self._is_forest(G):
                out.append((rem, comps[0]))
        return out

    # ------------------------------------------------------------------ polynomials
    def U(self):
        return sum(prod(self.x[i] for i in rem) for rem in self.spanning_trees())

    def F0(self):
        tot = self.R(0)
        for rem, comp in self.two_forests():
            P = {}
            for vtx in comp:
                P = add(P, self.ext.get(vtx, {}))
            if not P:
                continue
            tot += prod(self.x[i] for i in rem) * self.kin.square_external(P)
        return tot

    def F(self):
        """Euclidean F = F_0 + U sum x m^2.  (Minkowski: F = -F_0 + U sum x m^2, with
        Minkowski invariants; see IntegralFamily.UF.)"""
        U = self.U()
        mass = sum(self.x[i] * self.kin.R(m2) for i, (_, _, m2) in enumerate(self.lines))
        if self.kin.euclidean:
            return self.F0() + U * mass
        return -self.F0() + U * mass

    def U_kirchhoff(self):
        r"""
        Matrix-tree theorem: with edge weights 1/x_e, the reduced weighted Laplacian has
        determinant  sum over spanning trees of prod_{e in tree} 1/x_e,  so
        U = (prod_e x_e) * det(reduced Laplacian).
        """
        F = self.R.fraction_field()
        idx = {v: n for n, v in enumerate(self.vertices)}
        Lap = [[F(0)] * self.V for _ in range(self.V)]
        for i, (u, v, _) in enumerate(self.lines):
            if u == v:
                continue
            w = 1 / F(self.x[i])
            a, b = idx[u], idx[v]
            Lap[a][a] += w; Lap[b][b] += w
            Lap[a][b] -= w; Lap[b][a] -= w
        red = matrix(F, [row[1:] for row in Lap[1:]])
        val = red.det() * prod(self.x)
        return self.R(val.numerator()) / self.R(val.denominator()) if val.denominator() != 1 else self.R(val)

    # ------------------------------------------------------------------ momentum routing
    def contract(self, lines):
        """The graph with the given lines (numbered from 1) shrunk to points.  This is what a
        sector with index 0 on those lines means: the propagator is absent, so its two end
        vertices become one.  Returns (graph, original numbers of the remaining lines)."""
        lines = set(lines)
        parent = {v: v for v in self.vertices}

        def find(v):
            while parent[v] != v:
                parent[v] = parent[parent[v]]
                v = parent[v]
            return v
        for i, (u, v, _) in enumerate(self.lines):
            if (i + 1) in lines:
                ru, rv = find(u), find(v)
                if ru != rv:
                    parent[rv] = ru
        new_lines, kept = [], []
        for i, (u, v, m2) in enumerate(self.lines):
            if (i + 1) not in lines:
                new_lines.append((find(u), find(v), m2))
                kept.append(i + 1)
        ext = {}
        for v, q in self.ext.items():
            ext[find(v)] = add(ext.get(find(v), {}), q)
        ext = {v: q for v, q in ext.items() if q}
        return FeynmanGraph(new_lines, ext, self.kin), kept

    def plot(self, **kw):
        """The diagram in the style of the notes (see plotting.draw_graph); returns a matplotlib Figure."""
        from .plotting import draw_graph
        return draw_graph(self, **kw)

    def family(self, name='fam', loop_names=None):
        r"""
        Route momenta: pick a spanning tree T; every line not in T (a chord) gets a loop
        momentum; each chord's loop momentum flows back to its start along the unique
        tree path; external momenta flow through the tree by momentum conservation.
        Returns (propagators, loop names), ready for IntegralFamily.
        """
        from .family import IntegralFamily
        rem = self.spanning_trees()[0]
        chords = list(rem)
        tree = [i for i in range(self.N) if i not in rem]
        loops = loop_names or ['l%d' % (r + 1) for r in range(len(chords))]
        mom = [dict() for _ in range(self.N)]
        T = self._graph(tree)
        # loops: chord i from u to v carries l_r (u -> v); close it through the tree v -> u
        for r, i in enumerate(chords):
            u, v, _ = self.lines[i]
            mom[i] = add(mom[i], {loops[r]: QQ(1)})
            path = T.shortest_path(v, u, by_weight=False)
            for a, b in zip(path[:-1], path[1:]):
                e = [lab for (p, q, lab) in T.edges_incident(a) if {p, q} == {a, b}][0]
                eu, ev, _ = self.lines[e]
                s = 1 if (eu, ev) == (a, b) else -1
                mom[e] = add(mom[e], {loops[r]: QQ(s)})
        # externals: tree line e (eu -> ev) carries the external momentum entering the component of eu
        for e in tree:
            eu, ev, _ = self.lines[e]
            G2 = self._graph([f for f in tree if f != e])
            comp = G2.connected_component_containing_vertex(eu)
            P = {}
            for vtx in comp:
                P = add(P, self.ext.get(vtx, {}))
            mom[e] = add(mom[e], P)
        props = [(mom[i], self.lines[i][2]) for i in range(self.N)]
        return props, loops
