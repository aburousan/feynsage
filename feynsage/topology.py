r"""
Loop integrals as graphs, with Sage's graph and matroid theory (FeynCalc's FCLoopIntegralToGraph,
FCLoopToPakForm, FCLoopFindTopologyMappings, FCLoopFindMomentumShifts, FCLoopFindSubtopologies,
FCLoopScalelessQ and their neighbours).

    integral_to_graph(fam)          the Feynman graph of a family, its lines in the order of the propagators
    canonical_form(fam, sector)     a canonical form of U + F (Pak), from a canonical labelling of a coloured graph
    topology_mappings(fams)         which families (or sectors) are the same integrals after a shift of the loop momenta,
                                    with the propagator map and the shift
    momentum_shift(fam1, fam2, perm)    the loop-momentum change that turns the propagators of fam1 into those of fam2

How integral_to_graph works.  The loop-momentum part A_i of each propagator momentum is a vector in Q^L.  At a
vertex the loop momenta must cancel, so a vertex is a set of lines with signs s_i = +-1 and sum s_i A_i = 0.
The minimal such sets are the bonds of the graph: the circuits of the vector matroid of the A_i (the bond
matroid of the graph).  The vertices are a choice of bonds (or unions of bonds, at a cut vertex) in which every
line occurs exactly twice, once leaving and once entering: an exact cover, found with dancing links.  Lines
in no bond are self-loops.  The external momentum entering a vertex is what is left of sum s_i q_i.  The
result is checked: its U and F (spanning trees and 2-forests) must equal U and F of the family (completing
the square).  Integrals that factorize (a cut vertex) are graphs as well.
"""
from itertools import combinations

from sage.all import QQ, matrix, Matroid

from .graph import FeynmanGraph
from .momenta import add


class NotAGraphError(ValueError):
    """The propagators are not the lines of a Feynman graph (or their loop momenta are not independent)."""


def _loop_vectors(fam):
    A = [[QQ(q.get(l, 0)) for l in fam.loops] for q, _ in fam.props]
    for i, a in enumerate(A):
        if not any(a):
            raise NotAGraphError("propagator %d carries no loop momentum (a tree line): give it to the "
                                 "external kinematics instead" % (i + 1))
    if matrix(QQ, A).rank() != len(fam.loops):
        raise NotAGraphError("the loop momenta are not independent in these propagators")
    return A


def _signed(A, C):
    """The relation sum_{i in C} s_i A_i = 0 of a circuit C with s_i = +-1, or None if it is not of that form."""
    C = sorted(C)
    K = matrix(QQ, [A[i] for i in C]).left_kernel().basis()
    if len(K) != 1:
        return None
    k = K[0]
    a = abs(k[0])
    if a == 0 or any(abs(x) != a for x in k):
        return None
    return {i: (1 if x > 0 else -1) for i, x in zip(C, k)}


def _stars(A, N, ext_parts, limit=5000):
    """Vertex stars {line: +1 out, -1 in}: every line that lies in a bond occurs in two stars, once with each
    sign (an exact cover, dancing links), and a connected graph has N - L + 1 vertices.  Several covers can
    exist (graphs with the same bonds, Whitney's 2-isomorphisms; they have the same U and F); the one kept
    has the fewest vertices with external momenta and then the simplest external momenta."""
    from sage.combinat.matrices.dancing_links import dlx_solver
    L = len(A[0])
    M = Matroid(matrix=matrix(QQ, A).transpose())
    bonds = [s for s in (_signed(A, C) for C in M.circuits()) if s is not None]
    inbond = sorted({i for s in bonds for i in s})
    selfloops = [i for i in range(N) if i not in inbond]
    if not inbond:
        if N - L + 1 == 1:
            return [], selfloops                           # only self-loops: one vertex
        raise NotAGraphError("these propagators are not a graph")
    col = {}
    for i in inbond:
        col[(i, 1)], col[(i, -1)] = len(col), len(col) + 1
    want = N - L + 1

    def score(stars):
        ext = []
        for star in stars:
            P = {}
            for i, sg in star.items():
                P = add(P, {a: sg * c for a, c in ext_parts[i].items()})
            if P:
                ext.append(sum(abs(c) for c in P.values()))
        return (len(ext), sum(ext))

    cands = []
    for depth in (1, 2, 3):                                # unions of disjoint bonds: a cut vertex
        for combo in combinations(range(len(bonds)), depth):
            sets = [bonds[j] for j in combo]
            if len(set().union(*[set(s) for s in sets])) != sum(len(s) for s in sets):
                continue
            for signs in _sign_choices(depth):
                star = {}
                for s, e in zip(sets, signs):
                    star.update({i: e * v for i, v in s.items()})
                cands.append(star)
        rows = [sorted(col[(i, v)] for i, v in star.items()) for star in cands]
        best = None
        for n, sol in enumerate(dlx_solver(rows).solutions_iterator()):
            if len(sol) == want and _connected([cands[r] for r in sol]):
                stars = [cands[r] for r in sol]
                sc = score(stars)
                if best is None or sc < best[0]:
                    best = (sc, stars)
            if n >= limit:
                break
        if best is not None:
            return best[1], selfloops
    raise NotAGraphError("no set of vertices conserves the loop momenta: these propagators are not a graph")


def _connected(stars):
    """Do the stars (vertices) form one connected graph?  (A cover can also give several components that
    share a loop momentum: momentum is conserved, but it is not a routing.)"""
    ends = {}
    for k, star in enumerate(stars):
        for i in star:
            ends.setdefault(i, []).append(k)
    parent = list(range(len(stars)))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a
    for ks in ends.values():
        for k in ks[1:]:
            parent[find(k)] = find(ks[0])
    return len({find(k) for k in range(len(stars))}) <= 1


def _sign_choices(n):
    if n == 0:
        return [()]
    return [(1,) + rest for rest in _sign_choices(n - 1)] + [(-1,) + rest for rest in _sign_choices(n - 1)]


def integral_to_graph(fam, sector=None, check=True):
    """
    The Feynman graph of an integral family: a FeynmanGraph whose line i is propagator i (same momentum up to
    sign, same mass), with the external momenta attached where momentum conservation needs them.  Vertices
    are named v1, v2, ...; g.plot() draws it.  sector: a tuple of 0/1 or of powers (only the lines with a
    positive entry are drawn: numerators and missing propagators are not lines); the returned graph then has
    those lines only, in their order.  Raises NotAGraphError if the propagators are not a graph (or the
    sector is scaleless because a loop momentum is in none of its lines).  With check=True the graph's U and
    F (from its trees) are compared with the family's (from its momenta).
    """
    if sector is not None:
        from .family import IntegralFamily
        keep = [i for i in range(fam.t) if sector[i] > 0]
        fam = IntegralFamily(fam.name, fam.loops, fam.kin, [(fam.props[i][0], str(fam.props[i][1])) for i in keep])
    A = _loop_vectors(fam)
    N = len(A)
    ext_parts = [{a: c for a, c in q.items() if a not in fam.loops} for q, _ in fam.props]
    stars, selfloops = _stars(A, N, ext_parts)
    names = ["v%d" % (k + 1) for k in range(max(1, len(stars)))]
    tail, head = {}, {}
    for v, star in zip(names, stars):
        for i, s in star.items():
            (tail if s > 0 else head)[i] = v
    for i in selfloops:
        tail[i] = head[i] = names[0]
    lines = [(tail[i], head[i], fam.props[i][1]) for i in range(N)]
    ext = {}
    for v, star in zip(names, stars):
        P = {}
        for i, s in star.items():
            q = {a: c for a, c in fam.props[i][0].items() if a not in fam.loops}
            P = add(P, {a: s * c for a, c in q.items()})
        if P:
            ext[v] = P
    g = FeynmanGraph(lines, ext, fam.kin)
    if check:
        Um, Fm = fam.UF()
        if g.R(str(Um)) != g.U() or g.R(str(Fm)) != g.F():
            raise RuntimeError("integral_to_graph: the graph's U and F differ from the family's; please report")
    return g


# ---------------------------------------------------------------------------- canonical forms (Pak)
def _sector_lines(fam, sector, powers):
    if powers is not None:
        powers = tuple(powers)
        if any(a < 0 for a in powers):
            raise ValueError("canonical forms are for integrals without numerators (powers >= 0)")
        sector = tuple(1 if a > 0 else 0 for a in powers)
    if sector is None:
        sector = (1,) * fam.t
    return tuple(sector), powers


def canonical_form(fam, sector=None, powers=None):
    """
    A canonical form of the integrals of a sector (default: all lines), as a hashable key, and the map
    {line i: position in the canonical order}.  Two sectors (of the same or of different families) with equal
    keys are the same integrals up to a renaming of the Feynman parameters (Pak's criterion), so their
    integrals are equal.  With powers (all >= 0) the propagator powers are part of the form, so equal keys
    mean equal integrals.  The key comes from a canonical labelling (Sage's graph canonical_label) of a
    coloured graph: one vertex per Feynman parameter, one per monomial of U + F coloured by its coefficient,
    an edge for each power of a parameter in a monomial.  This is FeynCalc's FCLoopToPakForm, done by graph
    canonisation instead of Pak's sorting.
    """
    from sage.all import Graph
    sector, powers = _sector_lines(fam, sector, powers)
    U, F = fam.UF(sector)
    G = U + F
    lines = [i for i in range(fam.t) if sector[i]]
    gr = Graph(multiedges=False, loops=False)
    xs = [("x", i) for i in lines]
    gr.add_vertices(xs)
    classes = {}
    for i in lines:
        classes.setdefault(("x", str(powers[i]) if powers else "1"), []).append(("x", i))
    for j, (e, c) in enumerate(sorted(G.dict().items(), key=lambda t: str(t[0]))):
        v = ("m", j)
        gr.add_vertex(v)
        classes.setdefault(("m", str(c)), []).append(v)
        e = e if hasattr(e, "__iter__") else (e,)
        for i in lines:
            if e[i]:
                gr.add_edge(v, ("x", i), e[i])
    order = sorted(classes)
    partition = [classes[k] for k in order]
    C, cert = gr.canonical_label(partition=partition, edge_labels=True, certificate=True)
    key = (tuple((k, len(classes[k])) for k in order), tuple(sorted(C.edges(sort=False, labels=True))))
    return key, {i: cert[("x", i)] for i in lines}


def _match(fam1, s1, p1, fam2, s2, p2):
    """The line map i -> j with the integrals of (fam1, s1) equal to those of (fam2, s2), or None."""
    k1, c1 = canonical_form(fam1, s1, p1)
    k2, c2 = canonical_form(fam2, s2, p2)
    if k1 != k2:
        return None
    inv = {pos: j for j, pos in c2.items()}
    return {i: inv[pos] for i, pos in c1.items()}


def momentum_shift(fam1, fam2, perm):
    """
    The change of loop momenta  l -> M l' + c  (c a combination of the external momenta, det M = +-1) that turns
    every propagator i of fam1 in perm into +- propagator perm[i] of fam2, as {loop of fam1: momentum dict in the
    loops of fam2 and the external momenta}; None if no such linear change exists (the integrals can still be
    equal: canonical_form compares U and F).  FeynCalc's FCLoopFindMomentumShifts.
    """
    import itertools
    from sage.all import vector
    L1, L2 = list(fam1.loops), list(fam2.loops)
    if len(L1) != len(L2):
        return None
    ext = sorted(set(fam1.kin.externals) | set(fam2.kin.externals))
    lines = sorted(perm)

    def parts(fam, loops, i):
        q = fam.props[i][0]
        return vector(QQ, [q.get(l, 0) for l in loops]), vector(QQ, [q.get(a, 0) for a in ext])
    A1 = {i: parts(fam1, L1, i) for i in lines}
    A2 = {i: parts(fam2, L2, perm[i]) for i in lines}
    basis = []
    for i in lines:                                       # L independent lines of fam1
        if matrix(QQ, [A1[j][0] for j in basis + [i]]).rank() == len(basis) + 1:
            basis.append(i)
        if len(basis) == len(L1):
            break
    if len(basis) < len(L1):
        return None
    B = matrix(QQ, [A1[i][0] for i in basis])
    for signs in itertools.product((1, -1), repeat=len(basis)):
        M = B.inverse() * matrix(QQ, [s * A2[i][0] for s, i in zip(signs, basis)])
        c = B.inverse() * matrix(QQ, [s * A2[i][1] - A1[i][1] for s, i in zip(signs, basis)])
        if abs(M.det()) != 1:
            continue
        ok = True
        for i in lines:                                   # q1_i(M l' + c) = +- q2_perm(i)(l')
            a, b = A1[i]
            new_loop, new_ext = a * M, a * c + b
            if not any(new_loop == s * A2[i][0] and new_ext == s * A2[i][1] for s in (1, -1)):
                ok = False
                break
        if ok:
            return {l: {**{k: M[r, j] for j, k in enumerate(L2) if M[r, j]},
                        **{a: c[r, j] for j, a in enumerate(ext) if c[r, j]}} for r, l in enumerate(L1)}
    return None


def topology_mappings(families, preferred=()):
    """
    Which families are the same integrals as another family or as a sector of a larger one (FeynCalc's
    FCLoopFindTopologyMappings, with PreferredTopologies -> preferred).  Families are compared by canonical
    forms of U + F; the preferred ones, then the larger ones, are tried first as targets.  Returns
    (mappings, independent): mappings is a list of dicts {"family", "target", "sector" (of the target),
    "lines" {line of family: line of target}, "shift" (momentum_shift, or None)}; independent lists the
    topologies that remain: the preferred ones that are used, then the families that map to nothing.
    """
    from itertools import combinations
    cache = {}

    def keys(fam, n):
        """canonical keys of all sectors of fam with n lines"""
        if (id(fam), n) not in cache:
            d = {}
            for keep in combinations(range(fam.t), n):
                sec = tuple(1 if i in keep else 0 for i in range(fam.t))
                d.setdefault(canonical_form(fam, sec)[0], sec)
            cache[(id(fam), n)] = d
        return cache[(id(fam), n)]

    targets = list(preferred)
    rest = sorted((f for f in families if all(f is not p for p in preferred)), key=lambda f: -f.t)
    mappings, independent = [], list(preferred)
    for f in rest:
        key = canonical_form(f)[0]
        hit = None
        for tf in [t for t in targets if t.t == f.t] + [t for t in targets if t.t > f.t]:
            sec = keys(tf, f.t).get(key)                  # a family of the same size first, then sectors
            if sec is not None:
                hit = (tf, sec)
                break
        if hit is None:
            targets.append(f)
            independent.append(f)
            continue
        tf, sec = hit
        perm = _match(f, None, None, tf, sec, None)
        mappings.append({"family": f, "target": tf, "sector": sec, "lines": perm,
                         "shift": momentum_shift(f, tf, perm)})
    used = [t for t in independent if any(t is m["target"] for m in mappings) or all(t is not p for p in preferred)]
    return mappings, used


def map_integral(fam, powers, targets):
    """
    An integral fam(powers) (powers >= 0) written as an integral of one of the target families with permuted
    powers: (target family, its powers), or None.  FeynCalc's FCLoopFindIntegralMappings / FCLoopApplyTopologyMappings.
    """
    for tf in targets:                                    # first the map of the whole family (as FeynCalc)
        perm = _match(fam, None, None, tf, None, None) if tf.t == fam.t else None
        if perm is not None:
            new = [0] * tf.t
            for i, j in perm.items():
                new[j] = powers[i]
            if canonical_form(fam, powers=powers)[0] == canonical_form(tf, powers=new)[0]:
                return tf, tuple(new)
    for tf in targets:
        for perm_powers in _candidate_powers(tf, powers):
            perm = _match(fam, None, powers, tf, None, perm_powers)
            if perm is not None:
                new = [0] * tf.t
                for i, j in perm.items():
                    new[j] = powers[i]
                return tf, tuple(new)
    return None


def _candidate_powers(tf, powers):
    """Power vectors of the target with the same multiset of non-zero powers: the target sectors to try."""
    from itertools import combinations, permutations
    nz = [a for a in powers if a > 0]
    seen = set()
    for lines in combinations(range(tf.t), len(nz)):
        for arr in set(permutations(nz)):
            p = [0] * tf.t
            for i, a in zip(lines, arr):
                p[i] = a
            p = tuple(p)
            if p not in seen:
                seen.add(p)
                yield p


# ---------------------------------------------------------------------------- sectors
def subtopologies(fam, top=None):
    """
    The non-scaleless sectors below top (default: all lines), one per class of equal integrals (canonical
    forms), from the largest down: a list of sectors (tuples of 0/1).  FeynCalc's FCLoopFindSubtopologies,
    with FCLoopScalelessQ done by Lee's criterion (IntegralFamily.is_zero_sector).
    """
    from itertools import combinations
    top = tuple(top) if top is not None else (1,) * fam.t
    lines = [i for i in range(fam.t) if top[i]]
    seen, out = set(), []
    for n in range(len(lines), 0, -1):
        for keep in combinations(lines, n):
            sec = tuple(1 if i in keep else 0 for i in range(fam.t))
            if fam.is_zero_sector(sec):
                continue
            key, _ = canonical_form(fam, sec)
            if key in seen:
                continue
            seen.add(key)
            out.append(sec)
    return out


def scaleless(fam, sector=None):
    """True if the integrals of the sector vanish in dimensional regularisation (Lee's criterion)."""
    return fam.is_zero_sector(tuple(sector) if sector is not None else (1,) * fam.t)


# ---------------------------------------------------------------------------- the propagator basis
def _rows(fam):
    """D_i = sum_j C_ij sp_j + c_i over the loop scalar products of the family."""
    sign = 1 if fam.kin.euclidean else -1
    C, c = [], []
    for k, m2 in fam.props:
        co, ext = fam.dot(k, k)
        C.append(list(co))
        c.append(ext + sign * m2)
    return matrix(fam.kin.K, C), c


def basis_status(fam):
    """(rank, number of loop scalar products, number of propagators): the propagators are a basis when all three
    are equal; rank < scalar products is incomplete (FCLoopBasisIncompleteQ), rank < propagators is
    overdetermined, a linear relation among them (FCLoopBasisOverdeterminedQ)."""
    C, _ = _rows(fam)
    return C.rank(), fam.nsp, fam.t


def complete_basis(fam, candidates=None, name=None):
    """
    A family with propagators added (irreducible numerators for IBP) so that they are a basis of the loop
    scalar products; the added ones come last.  candidates: momenta (strings such as "l1 + p") to choose from,
    default l_i, l_i - l_j, l_i + p, l_i - p, l_i + l_j, ...  FeynCalc's FCLoopBasisFindCompletion.
    """
    from .family import IntegralFamily
    from ._parse import sr
    from sage.all import SR
    C, _ = _rows(fam)
    if C.rank() < fam.t:
        raise ValueError("the propagators are linearly dependent: partial fraction them first (apart)")
    loops, ext = fam.loops, list(fam.kin.externals)
    if candidates is None:
        cand = [{l: 1} for l in loops]
        cand += [{a: 1, b: -1} for a, b in combinations(loops, 2)]
        cand += [{l: 1, p: s} for l in loops for p in ext for s in (1, -1)]
        cand += [{a: 1, b: 1} for a, b in combinations(loops, 2)]
        cand += [{l: 1, p: 1, q: 1} for l in loops for p, q in combinations(ext, 2)]
        cand += [{l: 1, p: 1, q: -1} for l in loops for p, q in combinations(ext, 2)]
    else:
        names = loops + ext
        cand = []
        for c in candidates:
            e = sr(c) if isinstance(c, str) else SR(c)
            cand.append({str(v): QQ(e.coefficient(SR.var(v))) for v in names if e.coefficient(SR.var(v)) != 0})
    props = list(fam.props)
    rows = [list(r) for r in C.rows()]
    for q in cand:
        if len(rows) == fam.nsp:
            break
        co, _ = fam.dot(q, q)
        if matrix(fam.kin.K, rows + [list(co)]).rank() > len(rows):
            rows.append(list(co))
            props.append((q, fam.kin.R(0)))
    if len(rows) < fam.nsp:
        raise ValueError("the candidates do not complete the basis; give more candidates")
    return IntegralFamily(name or fam.name, fam.loops, fam.kin, [(k, str(m2)) for k, m2 in props])


def apart(fam, powers, check=True):
    """
    Partial fractions for propagators that are linearly dependent: fam(powers) as {powers: coefficient}, every
    term with linearly independent positive-power propagators.  A relation sum c_i D_i = c_0 among the positive
    ones gives 1 = sum c_i D_i / c_0 if c_0 != 0, else D_j = -sum_{i != j} c_i D_i / c_j (Pak's algorithm).
    FeynCalc's ApartFF / FCLoopCreatePartialFractioningRules.  check=True compares both sides at random
    rational values of the scalar products and invariants.
    """
    K = fam.kin.K
    C, c0 = _rows(fam)
    memo = {}

    def go(a):
        a = tuple(a)
        if a in memo:
            return memo[a]
        S = [i for i in range(fam.t) if a[i] > 0]
        kern = matrix(K, [C.row(i) for i in S]).left_kernel().basis() if S else []
        if not kern:
            memo[a] = {a: K(1)}
            return memo[a]
        rel = {S[k]: x for k, x in enumerate(kern[0]) if x != 0}
        const = sum(x * c0[i] for i, x in rel.items())
        out = {}
        if const != 0:                                    # 1 = sum c_i D_i / c_0: each term loses one power
            terms = [(i, x / const) for i, x in rel.items()]
            for i, w in terms:
                b = list(a); b[i] -= 1
                for k, v in go(b).items():
                    out[k] = out.get(k, K(0)) + w * v
        else:                                             # D_j = -sum c_i D_i / c_j
            j = max(rel, key=lambda i: (a[i], -i))
            for i, x in rel.items():
                if i == j:
                    continue
                b = list(a); b[i] -= 1; b[j] += 1
                for k, v in go(b).items():
                    out[k] = out.get(k, K(0)) - x / rel[j] * v
        out = {k: v for k, v in out.items() if v != 0}
        memo[a] = out
        return out

    res = go(tuple(powers))
    if check:
        _check_apart(fam, powers, res, C, c0)
    return res


def _check_apart(fam, powers, res, C, c0):
    import random
    rnd = random.Random(3)
    K = fam.kin.K
    for _ in range(2):
        sp = [QQ(rnd.randint(-97, 97)) / rnd.randint(1, 13) for _ in range(fam.nsp)]
        vals = {g: QQ(rnd.randint(2, 50)) / 7 for g in fam.kin.R.gens()}
        sub = lambda e: K(e).numerator().subs(vals) / K(e).denominator().subs(vals) if hasattr(K(e), "numerator") else e
        D = [sum(sub(C[i, j]) * sp[j] for j in range(fam.nsp)) + sub(c0[i]) for i in range(fam.t)]

        def val(a):
            out = QQ(1)
            for i, x in enumerate(a):
                out *= D[i] ** (-x)
            return out
        lhs = val(powers)
        rhs = sum(sub(c) * val(a) for a, c in res.items())
        if lhs != rhs:
            raise RuntimeError("apart: the partial fractions do not reproduce the integrand; please report")


# ---------------------------------------------------------------------------- dimension shifts
def dimension_lower(fam, powers):
    """
    The integral in d - 2 dimensions as d-dimensional integrals (Tarasov): {powers: coefficient}.  In feynsage's
    normalisation (Euclidean measure d^d l / pi^(d/2) per loop) J^(d-2)(a) = sum over the monomials prod x_i^e_i
    of U of prod (a_i)_(e_i) J^(d)(a + e); Minkowski adds (-1)^L.  FeynCalc's FCLoopGLILowerDimension.
    """
    from sage.all import rising_factorial
    U, _ = fam.UF()
    sign = 1 if fam.kin.euclidean else (-1) ** len(fam.loops)
    out = {}
    for e, coef in U.dict().items():
        e = e if hasattr(e, "__iter__") else (e,)
        w = sign * coef
        for ai, ei in zip(powers, e):
            w *= rising_factorial(ai, ei)
        if w != 0:
            b = tuple(ai + ei for ai, ei in zip(powers, e))
            out[b] = out.get(b, 0) + w
    return {k: v for k, v in out.items() if v != 0}


def dimension_raise(fam, powers):
    """
    The integral in d + 2 dimensions as d-dimensional integrals (Lee's raising relation, Baikov form):
    J^(d+2)(a) = s 2^L / ((d - E - L + 1) ... (d - E)) / det G(p) * P(D_1, ..., D_t) J^(d)(a), where P is the Gram
    determinant of the loop and external momenta written in the propagators, each D_i lowering the power a_i,
    G(p) the Gram matrix of the external momenta and s = 1 (Euclidean) or (-1)^L (Minkowski).  Needs a complete
    basis (complete_basis).  FeynCalc's FCLoopGLIRaiseDimension.
    """
    from sage.all import PolynomialRing
    if not fam.complete:
        raise ValueError("dimension_raise needs a complete basis of propagators: use complete_basis(fam)")
    K = fam.kin.K
    d = K(fam.kin.R.gen(0))
    L, ext = len(fam.loops), list(fam.kin.externals)
    E = len(ext)
    R = PolynomialRing(K, ["D%d" % (i + 1) for i in range(fam.t)])
    Dv = R.gens()

    def sp(a, b):
        n = fam._sp_index(a, b)
        if n is None:
            return R(fam.kin.ext_sp(a, b))
        co, c = fam.sp_in_dens(n)
        return sum(R(x) * Dv[k] for k, x in enumerate(co)) + R(c)
    vecs = fam.loops + ext
    P = matrix(R, [[sp(a, b) for b in vecs] for a in vecs]).det()
    Gp = matrix(K, [[K(fam.kin.ext_sp(a, b)) for b in ext] for a in ext]).det() if E else K(1)
    s = 1 if fam.kin.euclidean else (-1) ** L
    pref = s * 2 ** L / Gp
    for j in range(1, L + 1):
        pref /= (d - E - L + j)
    out = {}
    for e, coef in P.dict().items():
        e = e if hasattr(e, "__iter__") else (e,)
        b = tuple(ai - ei for ai, ei in zip(powers, e))
        out[b] = out.get(b, K(0)) + pref * coef
    return {k: v for k, v in out.items() if v != 0}


# ---------------------------------------------------------------------------- cuts
def cuts(g, separate=None):
    """
    The cuts of a FeynmanGraph: sets of lines whose removal splits the graph into exactly two connected parts
    with external momenta on both sides.  separate: (vertices that must be on one side, vertices that must be on
    the other), default any split of the external vertices.  Each cut is (lines, (side A, side B), threshold),
    threshold the square of the sum of the masses of the cut lines (the normal threshold of that channel).
    """
    from sage.all import Graph, sqrt, SR
    verts = list(g.vertices)
    extv = [v for v in verts if g.ext.get(v)]
    out, seen = [], set()
    for n in range(1, len(verts)):
        for S in combinations(verts, n):
            S = set(S)
            T = set(verts) - S
            if separate is not None:
                a, b = set(separate[0]), set(separate[1])
                if not ((a <= S and b <= T) or (a <= T and b <= S)):
                    continue
            elif not (S & set(extv) and T & set(extv)):
                continue
            key = frozenset([frozenset(S), frozenset(T)])
            if key in seen:
                continue
            seen.add(key)
            ok = True
            for part in (S, T):
                G = Graph(multiedges=True, loops=True)
                G.add_vertices(list(part))
                G.add_edges([(u, v) for (u, v, _) in g.lines if u in part and v in part])
                ok = ok and G.is_connected()
            if not ok:
                continue
            lines = [i for i, (u, v, _) in enumerate(g.lines) if (u in S) != (v in S)]
            thr = (sum(sqrt(SR(str(g.lines[i][2]))) for i in lines) ** 2).canonicalize_radical()
            out.append((lines, (sorted(S), sorted(T)), thr))
    return out


def cuttable(g, masses, separate=None):
    """Can the graph be cut (cuts()) through lines none of which has one of the masses?  masses: the squared
    masses (or their symbols) that must not go on shell.  FeynCalc's FCGraphCuttableQ."""
    from sage.all import SR
    bad = {str(SR(m)) for m in masses} | {str(SR(m) ** 2) for m in masses}
    return any(all(str(g.lines[i][2]) not in bad for i in lines) for lines, _, _ in cuts(g, separate))


# ---------------------------------------------------------------------------- topologies of an amplitude
def find_topologies(integrals, kin=None, loops=None, complete=True, euclidean=False, merge_sectors=True):
    """
    The families (topologies) behind a list of loop integrals, as FeynCalc's FCLoopFindTopologies.  integrals:
    one entry per integral, the list of its propagators as strings or (momentum, mass) pairs (a propagator
    written twice is a squared propagator; tree propagators without a loop momentum are ignored).  Integrals
    whose propagators are the same as those of another integral, or of a sector of a larger one, after a shift
    of the loop momenta, share its topology (merge_sectors=False: only integrals with as many propagators, so a
    smaller integral keeps its own topology, as FeynCalc does).  Scaleless integrals are dropped (they vanish).
    Returns (topologies, rows): topologies is the list of families (with irreducible numerators added when
    complete=True); rows[k] is (index of the topology, powers, momentum shift) for integral k, or None if it is
    scaleless.  Linearly dependent propagators must be partial-fractioned first (apart).
    """
    from .easy import family
    from collections import Counter
    from ._parse import sr
    from sage.all import SR
    items = []
    for props in integrals:
        cnt = Counter()
        for pr in props:
            q, m = (pr if isinstance(pr, (tuple, list)) else (pr, 0))
            cnt[(str(sr(q)), str(sr(m)))] += 1
        items.append(cnt)

    lnames = (loops.replace(",", " ").split() if isinstance(loops, str) else list(loops)) if loops else None

    def is_loop(q):
        names = {str(v) for v in SR(sr(q)).variables()}
        return bool(names & set(lnames)) if lnames else any(n[0] in "lkq" for n in names)

    def fam_of(cnt, name):
        ps = [(q, m) for (q, m) in cnt if is_loop(q)]      # tree propagators are not part of the integral
        if not ps:
            return None, []
        f = family(ps, kin=kin, loops=loops, euclidean=euclidean, name=name)
        return f, [cnt[p] for p in ps]
    cands = []
    for k, cnt in enumerate(items):
        f, pw = fam_of(cnt, "t%d" % (k + 1))
        if f is None or scaleless(f):
            cands.append(None)
            continue
        cands.append((f, tuple(pw)))
    order = sorted([k for k in range(len(cands)) if cands[k] is not None], key=lambda k: -cands[k][0].t)
    tops, rows = [], [None] * len(cands)
    cache = {}
    for k in order:
        f, pw = cands[k]
        key = canonical_form(f)[0]
        hit = None
        for ti, tf in enumerate(tops):
            if tf.t < f.t or (tf.t > f.t and not merge_sectors):
                continue
            if (ti, f.t) not in cache:
                d = {}
                for keep in combinations(range(tf.t), f.t):
                    sec = tuple(1 if i in keep else 0 for i in range(tf.t))
                    if not tf.is_zero_sector(sec):
                        d.setdefault(canonical_form(tf, sec)[0], sec)
                cache[(ti, f.t)] = d
            sec = cache[(ti, f.t)].get(key)
            if sec is not None:
                hit = (ti, sec)
                break
        if hit is None:
            tops.append(f)
            ti, perm = len(tops) - 1, {i: i for i in range(f.t)}
        else:
            ti, sec = hit
            perm = _match(f, None, None, tops[ti], sec, None)
        new = [0] * tops[ti].t
        for i, j in perm.items():
            new[j] = pw[i]
        rows[k] = (ti, tuple(new), momentum_shift(f, tops[ti], perm))
    names = []
    for ti, tf in enumerate(tops):
        tf.name = "topology%d" % (ti + 1)
    if complete:
        full = []
        for tf in tops:
            if basis_status(tf)[0] < tf.nsp:
                full.append(complete_basis(tf, name=tf.name))
            else:
                full.append(tf)
        rows = [None if r is None else (r[0], r[1] + (0,) * (full[r[0]].t - len(r[1])), r[2]) for r in rows]
        tops = full
    return tops, rows
