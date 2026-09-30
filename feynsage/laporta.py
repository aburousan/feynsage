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
"""
from itertools import product as iproduct


def weight(a):
    pos = [x for x in a if x > 0]
    neg = [-x for x in a if x < 0]
    return (len(pos), sum(x - 1 for x in pos), sum(neg), tuple(x > 0 for x in a), tuple(a))


class Reducer:

    def __init__(self, family, symmetries=()):
        self.fam = family
        self.t = family.t
        self.syms = _closure([tuple(s) for s in symmetries], self.t)
        self._zero = {}
        self.rules = {}
        self.seen = set()
        self.K = family.kin.K

    # ------------------------------------------------------------------ helpers
    def sector(self, a):
        return tuple(1 if x > 0 else 0 for x in a)

    def is_zero(self, a):
        s = self.sector(a)
        if s not in self._zero:
            self._zero[s] = self.fam.is_zero_sector(s)
        return self._zero[s]

    def canon(self, a):
        """Smallest image of a under the symmetry group (by weight)."""
        imgs = [tuple(a[s[i]] for i in range(self.t)) for s in self.syms]
        return min(imgs, key=weight)

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
        n_used = 0
        for a in self.seeds(rmax, smax):
            for ident in self.fam.ibp(a):
                n_used += self.add_identity(ident)
        if verbose:
            print("independent identities used:", n_used)
        return self

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
