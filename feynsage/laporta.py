r"""
Laporta reduction of an IntegralFamily.

Integrals are ordered by complexity (number of lines present, then dots, then
numerator powers, then lexicographically).  Identities generated from seeds are
solved one by one for their most complicated integral; the rules are kept fully
reduced (back substitution), so at the end every integral that has a rule is
written in terms of integrals without rules: the master integrals.

Symmetries are given as permutations of the propagator list that map the family
to itself (for example the reflection l -> p - l of an equal-mass bubble).
Zero sectors are found with Lee's criterion (IntegralFamily.is_zero_sector).

Sector symmetries (sector_symmetries=True, the default) go further, as in LiteRed and
Kira: an integral without numerators depends only on G = U + F with the x_j of its
missing lines set to zero.  Two sectors whose polynomials agree after renaming the
x_i hold the same integrals (Pak's criterion), also when no symmetry of the whole
family relates them, for example the three tadpole-product sectors of an equal-mass
sunset or the automorphisms of the sunset sector itself.  Each sector polynomial is
put in canonical form with Sage's canonical graph labelling and every such integral is
mapped to one representative sector.  Integrals with numerators cannot be renamed like
this: for them the loop momenta are shifted (a unimodular linear map, external momenta
fixed or all reversed) so that the propagators of the sector go into those of the
representative, the numerators are written again through the propagators and the
result is added to the system as one more equation (symmetry_relations), as LiteRed does.
"""
from itertools import product as iproduct, combinations
from sage.all import prod


def weight(a):
    pos = [x for x in a if x > 0]
    neg = [-x for x in a if x < 0]
    return (len(pos), sum(x - 1 for x in pos), sum(neg), tuple(x > 0 for x in a), tuple(a))


class Reducer:

    def __init__(self, family, symmetries=(), sector_symmetries=True, numerator_relations="auto", top=None):
        self.fam = family
        self.t = family.t
        self.syms = _closure([tuple(s) for s in symmetries], self.t)
        self._zero = {}
        self.rules = {}
        self.seen = set()
        self.K = family.kin.K
        self.sector_symmetries = sector_symmetries
        # top sector (0/1 per line): seeds stay inside it and representatives are chosen inside it
        # when a symmetry allows, as Kira does with its top-level sectors.  None: no restriction.
        self.top = tuple(top) if top is not None else None
        # "auto": reduce without the momentum-shift relations first and add them only if a
        # master with numerators sits in a sector where they could relate it to others
        self.numerator_relations = numerator_relations
        self._use_relations = numerator_relations is True
        self._G = None
        self._secinfo = {}       # sector -> (key, {line: canonical position})
        self._rep = {}           # key -> (representative sector, {canonical position: line}, automorphisms)
        self._canon_cache = {}
        self._sizes_done = set()
        self._bycheap = {}       # quick fingerprint -> sectors
        self._maps = {}          # (sector, frozen line map) -> momentum map or None
        self._relmaps = {}       # sector -> line maps that need momentum-shift relations

    # ------------------------------------------------------------------ helpers
    def sector(self, a):
        return tuple(1 if x > 0 else 0 for x in a)

    def is_zero(self, a):
        s = self.sector(a)
        if s not in self._zero:
            self._zero[s] = self.fam.is_zero_sector(s)
        return self._zero[s]

    def set_top(self, integrals):
        """Make the union of the sectors of these integrals the top sector (and forget every
        choice of representative made so far)."""
        self.top = tuple(1 if any(b[i] > 0 for b in integrals) else 0 for i in range(self.t))
        self._canon_cache, self._rep, self._relmaps = {}, {}, {}

    def _outside(self, a):
        """How many lines outside the top sector carry a propagator in a."""
        if self.top is None:
            return 0
        return sum(1 for i in range(self.t) if a[i] > 0 and not self.top[i])

    def _key(self, a):
        return (self._outside(a), weight(a))

    def canon(self, a):
        """The representative of a: for an integral without numerators the smallest image
        (by weight) inside its representative sector (sector symmetries); otherwise the
        smallest image under the symmetry group of the family."""
        a = tuple(a)
        if a in self._canon_cache:
            return self._canon_cache[a]
        if self.sector_symmetries and min(a) >= 0 and max(a) > 0:
            b = self._sector_canon(a)
        else:
            b = min((tuple(a[s[i]] for i in range(self.t)) for s in self.syms), key=self._key)
        self._canon_cache[a] = b
        return b

    # ------------------------------------------------------------------ sector symmetries
    def _monomials(self):
        """The monomials of G = U + F once, as an exponent array and coefficient classes."""
        if self._G is None:
            import numpy as np
            U, F = self.fam.UF()
            self._G = U + F
            items = sorted(((tuple(e) if hasattr(e, '__iter__') else (e,)), str(c)) for e, c in self._G.dict().items())
            names = sorted({c for _, c in items})
            self._cid = {c: n for n, c in enumerate(names)}
            self._E = np.array([e for e, _ in items], dtype=np.int64).reshape(-1, self.t)
            self._C = np.array([self._cid[c] for _, c in items], dtype=np.int64)
        return self._E, self._C

    def _restrict(self, sec):
        """Rows of the monomials that survive when the x_j of the missing lines are zero."""
        import numpy as np
        E, C = self._monomials()
        off = [j for j in range(self.t) if not sec[j]]
        keep = np.all(E[:, off] == 0, axis=1) if off else np.ones(len(E), dtype=bool)
        return E[keep], C[keep]

    def _cheap(self, sec):
        """A quick fingerprint that equivalent sectors share: the monomials (coefficient class,
        sorted exponents) and the exponent profile of every variable, each as a sorted list."""
        E, C = self._restrict(sec)
        lines = [i for i in range(self.t) if sec[i]]
        mons = tuple(sorted((int(c), tuple(sorted(int(x) for x in e if x))) for e, c in zip(E, C)))
        vars_ = tuple(sorted(tuple(sorted((int(E[r, i]), int(C[r])) for r in range(len(E)) if E[r, i])) for i in lines))
        return (sum(sec), mons, vars_)

    def _sector_graph(self, sec):
        """Coloured graph of G restricted to the sector: one vertex per present x_i, one per
        monomial (coloured by its coefficient), edges labelled by the exponents."""
        from sage.all import Graph
        E, C = self._restrict(sec)
        gr = Graph(multiedges=False, loops=False)
        xs = [('x', i) for i in range(self.t) if sec[i]]
        gr.add_vertices(xs)
        colours = {}
        for j in range(len(E)):
            v = ('m', j)
            gr.add_vertex(v)
            colours.setdefault(int(C[j]), []).append(v)
            for i in range(self.t):
                if E[j, i]:
                    gr.add_edge(v, ('x', i), int(E[j, i]))
        cols = sorted(colours)
        return gr, [xs] + [colours[c] for c in cols], tuple((c, len(colours[c])) for c in cols)

    def _info(self, sec):
        if sec not in self._secinfo:
            gr, part, cols = self._sector_graph(sec)
            C, cert = gr.canonical_label(partition=part, certificate=True, edge_labels=True)
            key = (sum(sec), cols, tuple(sorted(C.edges(labels=True, sort=True))))
            pos = {i: cert[('x', i)] for i in range(self.t) if sec[i]}
            self._secinfo[sec] = (key, pos)
        return self._secinfo[sec]

    def _representative(self, key, k, sec0=None):
        """The sector of smallest weight with this canonical form, its positions -> lines and
        the automorphisms of its polynomial (permutations of its lines).  Only sectors with
        the same quick fingerprint are compared (the fingerprints of all sectors with k lines
        are found once, the first time that size is met)."""
        if k not in self._sizes_done:
            for lines in combinations(range(self.t), k):
                sec = tuple(1 if i in lines else 0 for i in range(self.t))
                self._bycheap.setdefault(self._cheap(sec), []).append(sec)
            self._sizes_done.add(k)
        if key not in self._rep:
            cands = self._bycheap[self._cheap(sec0)] if sec0 is not None else \
                [sec for group in self._bycheap.values() for sec in group if sum(sec) == k]
            best = min((sec for sec in cands if self._info(sec)[0] == key), key=self._key)
            gr, part, _ = self._sector_graph(best)
            A = gr.automorphism_group(partition=part, edge_labels=True)
            perm_of = lambda g: {i: g(('x', i))[1] for i in range(self.t) if best[i]}
            auts = [perm_of(g) for g in A]
            gens = [perm_of(g) for g in A.gens()]
            line_at = {p: i for i, p in self._info(best)[1].items()}
            self._rep[key] = (best, line_at, auts, gens)
        return self._rep[key]

    # ------------------------------------------------------------------ numerators: momentum shifts
    def _global_images(self, sec):
        return {tuple(sec[g[i]] for i in range(self.t)) for g in self.syms}

    def _global_auts(self, sec):
        """The permutations of the lines of sec that come from symmetries of the whole family."""
        out = []
        for g in self.syms:
            if tuple(sec[g[i]] for i in range(self.t)) == sec:
                out.append({g[i]: i for i in range(self.t) if sec[i]})
        return out

    def _external_maps(self):
        """Signed permutations M of the external momenta (p_c -> sum_e M[c][e] p_e) that keep
        every scalar product p_a . p_b: the integrals depend on the momenta only through these."""
        if not hasattr(self, '_extmaps'):
            from sage.all import matrix, QQ
            from itertools import permutations, product as prod_
            exts = self.fam.kin.externals
            E = len(exts)
            sp = lambda a, b: self.fam.kin.ext_sp(exts[a], exts[b])
            out = []
            for pm in permutations(range(E)):
                for sg in prod_((1, -1), repeat=E):
                    M = matrix(QQ, E, E)
                    for c in range(E):
                        M[c, pm[c]] = sg[c]
                    if all(sg[a] * sg[b] * sp(pm[a], pm[b]) == sp(a, b) for a in range(E) for b in range(a, E)):
                        out.append(M)
            out.sort(key=lambda M: (M != 1, M != -1))        # identity first, then p -> -p
            self._extmaps = out
        return self._extmaps

    def _momentum_map(self, sec, perm):
        """A loop-momentum map l -> A l + B p (det A = +-1), with the external momenta mapped by
        a signed permutation that keeps their scalar products, which sends the propagator of
        every line i of sec into the propagator of line perm[i].  Returns the transformed
        momenta of all lines, or None."""
        from sage.all import matrix, vector, QQ
        from itertools import product as prod_
        key = (sec, tuple(sorted(perm.items())))
        if key in self._maps:
            return self._maps[key]
        fam = self.fam
        loops, exts = fam.loops, fam.kin.externals
        L = len(loops)
        C = [vector(QQ, [k.get(l, 0) for l in loops]) for k, _ in fam.props]
        P = [vector(QQ, [k.get(e, 0) for e in exts]) for k, _ in fam.props]
        lines = [i for i in range(self.t) if sec[i]]
        if any(fam.props[i][1] != fam.props[perm[i]][1] for i in lines):
            self._maps[key] = None
            return None
        basis = []
        for i in lines:
            if matrix(QQ, [C[j] for j in basis + [i]]).rank() == len(basis) + 1:
                basis.append(i)
            if len(basis) == L:
                break
        result = None
        if len(basis) == L:
            CI = matrix(QQ, [C[i] for i in basis]).inverse()
            for M in (self._external_maps() if exts else [None]):
                for signs in prod_((1, -1), repeat=L):
                    A = CI * matrix(QQ, [sg * C[perm[i]] for sg, i in zip(signs, basis)])
                    if abs(A.det()) != 1:
                        continue
                    B = CI * matrix(QQ, [sg * P[perm[i]] - P[i] * M for sg, i in zip(signs, basis)]) if exts else None
                    ok = True
                    for i in lines:
                        nl = C[i] * A
                        npv = (C[i] * B + P[i] * M) if exts else None
                        if not any(nl == sg * C[perm[i]] and (not exts or npv == sg * P[perm[i]]) for sg in (1, -1)):
                            ok = False
                            break
                    if ok:
                        new = []
                        for j in range(self.t):
                            q = {l: c for l, c in zip(loops, C[j] * A) if c != 0}
                            if exts:
                                q.update({e: c for e, c in zip(exts, C[j] * B + P[j] * M) if c != 0})
                            new.append(q)
                        result = new
                        break
                if result is not None:
                    break
        self._maps[key] = result
        return result

    def _linear_form(self, q, m2):
        """(q^2 -+ m^2) written as sum_k c_k D_k + c_0 (coefficients in the fraction field)."""
        fam = self.fam
        sign = 1 if fam.kin.euclidean else -1
        co, ext = fam.dot(q, q)
        lin = [self.K(0)] * self.t
        const = self.K(ext + sign * m2)
        for n_, cn in enumerate(co):
            if cn == 0:
                continue
            dens, c0 = fam.sp_in_dens(n_)
            const += cn * c0
            for k_, ck in enumerate(dens):
                lin[k_] += cn * ck
        return lin, const

    def _mapped(self, a, sec, perm, newmom):
        """I(a) written in the sector perm(sec) through the momentum map: {index: coefficient},
        coefficients in the polynomial ring of the kinematics.  The arithmetic runs in one flat
        polynomial ring (y_1..y_t and the kinematic variables); the numerators written through
        the propagators, and their powers, are cached per map."""
        from sage.all import PolynomialRing, QQ
        R = self.fam.kin.R
        if not hasattr(self, '_YR'):
            self._kn = R.ngens()
            self._YR = PolynomialRing(QQ, ['y%d' % i for i in range(self.t)] + ['kk%d' % i for i in range(self._kn)])
            self._kin_in_Y = R.hom(self._YR.gens()[self.t:], self._YR)
            self._powcache = {}
        Y = self._YR
        y = Y.gens()[:self.t]
        mkey = (sec, tuple(sorted(perm.items())))
        base = [0] * self.t
        for i in range(self.t):
            if sec[i]:
                base[perm[i]] = a[i]
        num = Y(1)
        for j in range(self.t):
            if a[j] < 0:
                ck = (mkey, j, -a[j])
                if ck not in self._powcache:
                    fk = (mkey, j, 1)
                    if fk not in self._powcache:
                        lin, const = self._linear_form(newmom[j], self.fam.props[j][1])
                        conv = lambda c: self._kin_in_Y(R(c.numerator())) / QQ(c.denominator())
                        self._powcache[fk] = sum(conv(c) * y[k] for k, c in enumerate(lin) if c != 0) + conv(const)
                    self._powcache[ck] = self._powcache[fk] ** (-a[j])
                num *= self._powcache[ck]
        groups = {}
        for e, c in num.dict().items():
            b = tuple(base[k] - e[k] for k in range(self.t))
            groups.setdefault(b, {})[tuple(e[self.t:])] = c
        out = {}
        uni = R.ngens() == 1
        for b, dct in groups.items():
            c = R({k[0]: v for k, v in dct.items()}) if uni else R(dct)
            if c != 0:
                out[b] = c
        return out

    def _relation_maps(self, sec):
        """The line maps that need a momentum shift for integrals with numerators in sec:
        into its representative sector, or the automorphisms of the representative that the
        symmetries of the whole family do not give."""
        if sec in self._relmaps:
            return self._relmaps[sec]
        key, pos = self._info(sec)
        best, line_at, auts, gens = self._representative(key, sum(sec), sec)
        maps = []
        if sec != best:
            if best not in self._global_images(sec):
                # into the representative, then (if no shift realises that) through its automorphisms
                first = {i: line_at[pos[i]] for i in range(self.t) if sec[i]}
                for aut in auts:
                    cand = {i: aut[j] for i, j in first.items()}
                    if self._momentum_map(sec, cand) is not None:
                        maps.append(cand)
                        break
        else:
            known = self._global_auts(sec)
            maps += [g for g in gens if g not in known and any(g[i] != i for i in g)
                     and self._momentum_map(sec, g) is not None]
        self._relmaps[sec] = maps
        return maps

    def relations_needed(self, masters):
        """True if some master has numerators in a sector where the momentum-shift relations
        could relate it to other integrals (then the reduction is repeated with them)."""
        if not self.sector_symmetries:
            return False
        for m in masters:
            m = tuple(m)
            if min(m) < 0 and max(m) > 0 and not self.is_zero(m) and self._relation_maps(self.sector(m)):
                return True
        return False

    def symmetry_relations(self, a, polynomial=False):
        """Extra equations for an integral with numerators: I(a) = (the same integral after a
        momentum shift into its representative sector, or after an automorphism of that
        sector that the symmetries of the whole family do not give).  Each is a dict
        {index: coefficient} meaning sum c I = 0, the coefficients in the fraction field (or,
        with polynomial=True, in the polynomial ring of the kinematics).  Empty unless the
        relations are switched on (numerator_relations=True, or "auto" once they are needed)."""
        a = tuple(a)
        if not (self.sector_symmetries and self._use_relations) or min(a) >= 0 or max(a) <= 0:
            return []
        if self.is_zero(a):
            return []
        sec = self.sector(a)
        rels = []
        R = self.fam.kin.R
        for perm in self._relation_maps(sec):
            newmom = self._momentum_map(sec, perm)
            rel = {b: -c for b, c in self._mapped(a, sec, perm, newmom).items()}
            rel[a] = rel.get(a, R(0)) + 1
            rel = {b: c for b, c in rel.items() if c != 0}
            if rel:
                rels.append(rel if polynomial else {b: self.K(c) for b, c in rel.items()})
        return rels

    def _sector_canon(self, a):
        sec = self.sector(a)
        key, pos = self._info(sec)
        best, line_at, auts, _ = self._representative(key, sum(sec), sec)
        a0 = [0] * self.t
        for i in range(self.t):
            if sec[i]:
                a0[line_at[pos[i]]] = a[i]
        imgs = []
        for perm in auts:
            img = [0] * self.t
            for i, j in perm.items():
                img[j] = a0[i]
            imgs.append(tuple(img))
        return min(imgs, key=self._key)

    def clean(self, expr):
        out = {}
        for a, c in expr.items():
            if c == 0 or self.is_zero(a):
                continue
            b = self.canon(a)
            out[b] = out.get(b, self.K(0)) + c
            if out[b] == 0:
                del out[b]
        return out

    def substitute(self, expr):
        out = {}
        for a, c in expr.items():
            if a in self.rules:
                for b, cb in self.rules[a].items():
                    out[b] = out.get(b, self.K(0)) + c * cb
            else:
                out[a] = out.get(a, self.K(0)) + c
        return {a: c for a, c in out.items() if c != 0}

    # ------------------------------------------------------------------ elimination
    def add_identity(self, ident):
        c = self.clean(ident)
        self.seen.update(c.keys())
        e = self.substitute(c)
        if not e:
            return False
        top = max(e, key=weight)
        ct = e.pop(top)
        rule = {a: -c / ct for a, c in e.items()}
        # back substitution into existing rules
        for a in list(self.rules):
            r = self.rules[a]
            if top in r:
                c = r.pop(top)
                for b, cb in rule.items():
                    r[b] = r.get(b, self.K(0)) + c * cb
                    if r[b] == 0:
                        del r[b]
        self.rules[top] = rule
        return True

    def seeds(self, rmax, smax):
        """Seeds sector by sector, simplest sectors first."""
        sectors = [s for s in iproduct((0, 1), repeat=self.t) if sum(s) > 0]
        sectors.sort(key=lambda s: (sum(s), s))
        out = []
        for s in sectors:
            if self.top is not None and any(s[i] and not self.top[i] for i in range(self.t)):
                continue                               # outside the top sector
            if self.fam.is_zero_sector(s):
                continue
            pos = [i for i in range(self.t) if s[i]]
            neg = [i for i in range(self.t) if not s[i]]
            for r in range(rmax + 1):
                for dots in _compositions(r, len(pos)):
                    for sn in range(smax + 1):
                        for nums in _compositions(sn, len(neg)):
                            a = [0] * self.t
                            for i, dd in zip(pos, dots):
                                a[i] = 1 + dd
                            for i, nn in zip(neg, nums):
                                a[i] = -nn
                            out.append(tuple(a))
        return out

    def run(self, rmax, smax=0, verbose=False):
        self.rmax, self.smax = rmax, smax
        self._run(rmax, smax, verbose)
        if self.numerator_relations == "auto" and not self._use_relations and self.relations_needed(self.masters()):
            if verbose:
                print("masters with numerators in symmetric sectors: again with momentum-shift relations")
            self.rules, self.seen, self._use_relations = {}, set(), True
            self._run(rmax, smax, verbose)
        return self

    def _run(self, rmax, smax, verbose):
        n_used = 0
        for a in self.seeds(rmax, smax):
            for ident in self.fam.ibp(a) + self.symmetry_relations(a):
                n_used += self.add_identity(ident)
        if verbose:
            print("independent identities used:", n_used)

    def reduce(self, a):
        a = self.canon(tuple(a))
        if self.is_zero(a):
            return {}
        if a in self.rules:
            return dict(self.rules[a])
        return {a: self.K(1)}

    def masters(self, targets=None):
        """Master integrals.  With targets: exactly the integrals the targets reduce to.
        Without: the unreduced integrals inside the seeded range (at most rmax dots and
        smax numerator powers); integrals just outside it are left unreduced by
        construction and are not masters."""
        ms = set()
        if targets is not None:
            for t in targets:
                ms.update(self.reduce(t).keys())
            return sorted(ms, key=weight)
        for r in self.rules.values():
            ms.update(r.keys())
        ms.update(a for a in self.seen if a not in self.rules)
        inside = [a for a in ms
                  if sum(x - 1 for x in a if x > 0) <= self.rmax and sum(-x for x in a if x < 0) <= self.smax]
        return sorted(inside, key=weight)


def _compositions(n, k):
    """All k-tuples of non-negative integers summing to n."""
    if k == 0:
        if n == 0:
            yield ()
        return
    if k == 1:
        yield (n,)
        return
    for first in range(n + 1):
        for rest in _compositions(n - first, k - 1):
            yield (first,) + rest


def _closure(gens, t):
    """The group generated by the given permutations (always including the identity)."""
    ident = tuple(range(t))
    group = {ident}
    frontier = [ident]
    while frontier:
        new = []
        for g in frontier:
            for h in gens:
                comp = tuple(g[h[i]] for i in range(t))
                if comp not in group:
                    group.add(comp); new.append(comp)
        frontier = new
    return sorted(group)
