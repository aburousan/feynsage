r"""
Fast IBP reduction by finite-field sampling and rational reconstruction
(the strategy of FIRE and Kira/FireFly, in a small form).

1. The IBP system is generated once, with polynomial coefficients in the variables
   (d, and at most one kinematic invariant; set the other scales to numbers).
2. At a numerical sample point, modulo a prime p, the system is a sparse matrix over
   GF(p).  Sage row-reduces it in compiled code.  Columns are ordered from the most
   complicated integral to the simplest, so each target comes out as a combination
   of master integrals.  One row reduction gives *all* coefficients at that point.
3. Every coefficient is a rational function.  It is rebuilt by Thiele interpolation
   in d, and, if there is a second variable, each of its (normalised) coefficients is
   rebuilt again by Thiele interpolation in that variable.  Reconstruction stops only
   when two new, unused points agree with the interpolant.
4. The rational numbers are recovered from residues modulo several primes (Chinese
   remainder theorem + rational reconstruction); a result is accepted when one more
   prime no longer changes it.
"""
from sage.all import (GF, ZZ, QQ, matrix, previous_prime, PolynomialRing, crt,
                      lcm, prod)

import ctypes
import heapq
import os
import subprocess
import sys

from .laporta import weight


# ---------------------------------------------------------------------- the system
class IBPSystem:
    """All identities of a Reducer's seeds, with polynomial coefficients."""

    def __init__(self, reducer, rmax, smax=0):
        import numpy as np
        self.red = reducer
        fam = reducer.fam
        self.R = fam.kin.R
        self.nvars = self.R.ngens()
        if self.nvars > 2:
            raise ValueError("ff reduction handles d and at most one invariant: set the other scales to numbers")
        t = fam.t
        templates = ibp_templates(fam)
        # symmetry relations of the seeds with numerators (sector symmetries, momentum shifts)
        sym_rows = []
        for sd in reducer.seeds(rmax, smax):
            for rel in reducer.symmetry_relations(sd, polynomial=True):
                den = lcm([QQ(x).denominator() for c in rel.values() for x in c.coefficients()])
                sym_rows.append({b: _as_dict(self.R(c * den)) for b, c in rel.items()})
        # every coefficient is an integer combination of these monomials in (d, invariant)
        mons = sorted({e for tpl in templates for _, c0, cs in tpl for dct in [c0] + cs for e in dct}
                      | {e for row in sym_rows for dct in row.values() for e in dct})
        midx = {e: k for k, e in enumerate(mons)}
        nmon = len(mons)
        seeds = np.array(reducer.seeds(rmax, smax), dtype=np.int64).reshape(-1, t)
        N = len(seeds)
        # an index vector is packed into one integer, 6 bits per entry
        B, OFF = 6, 32
        if t * B > 62 or seeds.max(initial=0) + 2 >= OFF or -seeds.min(initial=0) + 2 >= OFF:
            raise ValueError("index vectors too long or too large for the packed representation")
        weights = np.array([1 << (B * i) for i in range(t)], dtype=np.int64)
        cid = {}                                   # packed index -> provisional column, -1 if zero
        canon_of = []                              # provisional column -> canonical index tuple
        cmap = {}
        R_, C_, E_ = [], [], []
        for ti, tpl in enumerate(templates):
            for shift, c0, cs in tpl:
                c0v = np.zeros(nmon, dtype=np.int64)
                for e, v in c0.items():
                    c0v[midx[e]] = v
                csm = np.zeros((t, nmon), dtype=np.int64)
                for a, ca in enumerate(cs):
                    for e, v in ca.items():
                        csm[a, midx[e]] = v
                coef = c0v[None, :] + seeds @ csm
                keep = np.any(coef != 0, axis=1)
                if not keep.any():
                    continue
                idx = seeds[keep] + np.array(shift, dtype=np.int64)
                keys = (idx + OFF) @ weights
                uk, inv = np.unique(keys, return_inverse=True)
                look = np.empty(len(uk), dtype=np.int64)
                for u, key in enumerate(uk.tolist()):
                    c = cid.get(key)
                    if c is None:
                        b = tuple(((key >> (B * i)) & 63) - OFF for i in range(t))
                        cb = None if reducer.is_zero(b) else reducer.canon(b)
                        if cb is None:
                            c = -1
                        else:
                            c = cmap.get(cb)
                            if c is None:
                                c = cmap[cb] = len(canon_of)
                                canon_of.append(cb)
                        cid[key] = c
                    look[u] = c
                cols = look[inv]
                nz = cols >= 0
                rows = (np.nonzero(keep)[0] * len(templates) + ti)[nz]
                R_.append(rows); C_.append(cols[nz]); E_.append(coef[keep][nz])
        base_row = N * len(templates)
        for k, row in enumerate(sym_rows):
            for b, dct in row.items():
                if not dct:
                    continue
                if max(b) + 2 >= OFF or -min(b) + 2 >= OFF:
                    raise ValueError("index vectors too large for the packed representation")
                key = sum((x + OFF) << (B * i) for i, x in enumerate(b))
                c = cid.get(key)
                if c is None:
                    cb = None if reducer.is_zero(b) else reducer.canon(b)
                    if cb is None:
                        c = -1
                    else:
                        c = cmap.get(cb)
                        if c is None:
                            c = cmap[cb] = len(canon_of)
                            canon_of.append(cb)
                    cid[key] = c
                if c < 0:
                    continue
                vec = np.zeros((1, nmon), dtype=np.int64)
                for e, v in dct.items():
                    vec[0, midx[e]] = v
                R_.append(np.array([base_row + k], dtype=np.int64)); C_.append(np.array([c], dtype=np.int64)); E_.append(vec)
        R_ = np.concatenate(R_); C_ = np.concatenate(C_); E_ = np.concatenate(E_)
        # final column order: most complicated integral first
        order = sorted(range(len(canon_of)), key=lambda c: weight(canon_of[c]), reverse=True)
        perm = np.empty(len(canon_of), dtype=np.int64)
        perm[np.array(order, dtype=np.int64)] = np.arange(len(canon_of))
        C_ = perm[C_]
        self.cols = [canon_of[c] for c in order]
        self.colidx = {c: i for i, c in enumerate(self.cols)}
        ncols = len(self.cols)
        # add up entries with the same (row, column), drop the ones that cancel
        key2 = R_ * ncols + C_
        o = np.argsort(key2, kind='stable')
        key2, E_ = key2[o], E_[o]
        starts = np.nonzero(np.r_[True, key2[1:] != key2[:-1]])[0]
        E_ = np.add.reduceat(E_, starts, axis=0)
        key2 = key2[starts]
        nz = np.any(E_ != 0, axis=1)
        key2, E_ = key2[nz], E_[nz]
        R_, C_ = key2 // ncols, key2 % ncols
        # rows: renumber, then simplest equations first (largest leading column), shorter first
        ur, R_ = np.unique(R_, return_inverse=True)
        nrows = len(ur)
        lead = np.full(nrows, ncols, dtype=np.int64); np.minimum.at(lead, R_, C_)
        length = np.bincount(R_, minlength=nrows)
        rorder = np.lexsort((length, -lead))
        rpos = np.empty(nrows, dtype=np.int64); rpos[rorder] = np.arange(nrows)
        R_ = rpos[R_]
        o = np.lexsort((C_, R_))
        self.ecols = np.ascontiguousarray(C_[o], dtype=np.int32)
        self.E = np.ascontiguousarray(E_[o], dtype=np.int64)
        self.rowptr = np.zeros(nrows + 1, dtype=np.int64)
        np.cumsum(np.bincount(R_, minlength=nrows), out=self.rowptr[1:])
        self.mons = mons
        self.nrows = nrows
        self._rows_int = None
        self._trimmed = False

    def _c_trim(self, p, pt, target_cols):
        """Indices of the rows the targets need, from one elimination in C."""
        import numpy as np
        mv = np.array([_mono_mod(e, pt, p) for e in self.mons], dtype=np.uint64)
        tg = np.array(target_cols, dtype=np.int32)
        out = np.zeros(self.nrows, dtype=np.int32)
        cnt = _clib().elim_trim(self.nrows, len(self.cols), self.rowptr.ctypes.data, self.ecols.ctypes.data,
                                self.E.ctypes.data, len(self.mons), mv.ctypes.data, p,
                                len(tg), tg.ctypes.data, out.ctypes.data)
        return out[:cnt].tolist()

    def _py_rows(self, which):
        """Rows as [(column, [(coefficient, exponent tuple), ...]), ...] for the Python path."""
        out = []
        for n in which:
            a, b = int(self.rowptr[n]), int(self.rowptr[n + 1])
            row = []
            for e in range(a, b):
                terms = [(int(v), self.mons[m]) for m, v in enumerate(self.E[e].tolist()) if v]
                row.append((int(self.ecols[e]), terms))
            out.append(row)
        return out

    def sample(self, targets, p, point):
        """{(target, master): value mod p} at one point, or None if the point is degenerate.

        Sparse Gaussian elimination on Python integers.  Column 0 is the most
        complicated integral, so each row is solved for its smallest column index
        (as in Laporta).  Forward elimination only; back substitution is done just for
        the rows the targets need."""
        if self._trimmed and not set(targets) <= self._trim_targets:
            raise ValueError("the system was trimmed for other targets")
        pt = [int(v) % p for v in point]
        if not self._trimmed:
            if p < 2**63 and _clib() is not None:
                tcols = [self.colidx.get(self.red.canon(tuple(t)), -1) for t in targets]
                self._rows_int = self._py_rows(self._c_trim(p, pt, tcols))
                self._trim_targets = set(targets)
                self._trimmed = True
            elif self._rows_int is None:
                self._rows_int = self._py_rows(range(self.nrows))
        piv = {}                                  # column -> row dict, normalised to 1 at column
        src = {}                                  # column -> (equation number, pivots used on it)
        for n, row in enumerate(self._rows_int):
            hit = set()
            r = {}
            for j, terms in row:
                v = _eval_mod(terms, pt, p)
                if v:
                    r[j] = v
            heap = list(r)
            heapq.heapify(heap)
            while heap:
                j = heapq.heappop(heap)
                if j not in r:                    # cancelled earlier, or a duplicate entry
                    continue
                if j in piv:
                    hit.add(j)
                    c = r.pop(j)
                    for k, v in piv[j].items():
                        if k != j:
                            if k in r:
                                w = (r[k] - c * v) % p
                                if w:
                                    r[k] = w
                                else:
                                    del r[k]
                            else:
                                r[k] = (-c * v) % p
                                heapq.heappush(heap, k)
                    continue
                inv = pow(r[j], -1, p)
                piv[j] = {k: v * inv % p for k, v in r.items()}
                src[j] = (n, hit)
                break
        if not self._trimmed:
            # keep only the equations the targets need (Kira's trick): the pivot rows of the
            # target columns, and recursively every pivot used to build or back-substitute
            # them.  A random point gives the generic structure with probability close to 1.
            need, stack = set(), []
            for t in targets:
                j = self.colidx.get(self.red.canon(tuple(t)))
                if j in piv:
                    stack.append(j)
            while stack:
                j = stack.pop()
                if j in need:
                    continue
                need.add(j)
                stack.extend(k for k in src[j][1] | set(piv[j]) if k in piv and k not in need)
            self._rows_int = [self._rows_int[n] for n in sorted(src[j][0] for j in need)]
            self._trim_targets = set(targets)
            self._trimmed = True
        memo = {}

        def solve(j):                             # integral j in terms of non-pivot columns
            if j in memo:
                return memo[j]
            out = {}
            for k, v in piv[j].items():
                if k == j:
                    continue
                if k in piv:
                    for m, w in solve(k).items():
                        out[m] = (out.get(m, 0) - v * w) % p
                else:
                    out[k] = (out.get(k, 0) - v) % p
            out = {m: w for m, w in out.items() if w}
            memo[j] = out
            return out

        sys.setrecursionlimit(max(10000, 4 * len(self.cols)))
        F = GF(p)
        out = {}
        for t in targets:
            tc = self.red.canon(tuple(t))
            if self.red.is_zero(tc):
                continue
            j = self.colidx.get(tc)
            if j is None:
                raise ValueError("target %s is outside the seeded range: raise rmax or smax" % (t,))
            if j not in piv:                      # the target is itself a master integral
                out[(t, tc)] = F(1)
                continue
            for k, v in solve(j).items():
                out[(t, self.cols[k])] = F(v)
        return out


def ibp_templates(fam):
    r"""
    The L(L+E) IBP identities with symbolic indices a_1..a_t, derived once.
    Each identity is a list of terms (shift, c0, [c_1..c_t]): the integral a + shift
    with coefficient c0 + sum_i a_i c_i.  The c's are integer polynomials in the
    variables of fam.kin.R, stored as {exponent tuple: int}; each identity is
    multiplied by one common factor so that all of them are integral.
    """
    if not fam.complete:
        raise ValueError("family is not complete: add irreducible numerators")
    K, R, t = fam.kin.K, fam.kin.R, fam.t
    d = K(R.gen(0))
    out = []
    for lr in fam.loops:
        for vname in fam.loops + fam.kin.externals:
            v = {vname: QQ(1)}
            terms = {}                         # shift -> [c0, c_1..c_t] in K

            def put(shift, which, c):
                row = terms.setdefault(tuple(shift), [K(0)] * (t + 1))
                row[which] += c
            if vname == lr:
                put([0] * t, 0, d)
            for i, (k, m2) in enumerate(fam.props):
                cir = k.get(lr, 0)
                if cir == 0:
                    continue
                co, ext = fam.dot(v, k)
                base = [0] * t; base[i] += 1
                pref = -2 * cir                # times a_i
                put(base, i + 1, pref * ext)
                for j, cj in enumerate(co):
                    if cj == 0:
                        continue
                    dens, const = fam.sp_in_dens(j)
                    put(base, i + 1, pref * cj * const)
                    for kk, ck in enumerate(dens):
                        if ck == 0:
                            continue
                        b = list(base); b[kk] -= 1
                        put(b, i + 1, pref * cj * ck)
            allc = [c for row in terms.values() for c in row if c != 0]
            if not allc:
                continue
            den = lcm([c.denominator() for c in allc])
            den = den * lcm([QQ(x).denominator() for c in allc for x in R(c * den).coefficients()])
            tpl = []
            for shift, row in terms.items():
                ints = [_as_dict(R(c * den)) for c in row]
                if any(ints):
                    tpl.append((shift, ints[0], ints[1:]))
            out.append(tpl)
    return out


def _as_dict(poly):
    """{exponent tuple: int} of an integer polynomial."""
    return {(tuple(e) if hasattr(e, '__iter__') else (e,)): int(c) for e, c in poly.dict().items()}


def _mono_mod(e, pt, p):
    t = 1
    for x, k in zip(pt, e):
        t = t * pow(x, k, p) % p
    return t


def _eval_mod(terms, pt, p):
    s = 0
    for c, e in terms:
        t = c % p
        for x, k in zip(pt, e):
            if k:
                t = t * pow(x, k, p) % p
        s += t
    return s % p


# ---------------------------------------------------------------------- compiled first probe
_LIB = None


def _clib():
    """Compile csrc/elim.c once (with the system C compiler) and load it; None if that fails."""
    global _LIB
    if _LIB is not None:
        return _LIB or None
    here = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'csrc')
    src = os.path.join(here, 'elim.c')
    name = 'libelim' + ('.dylib' if sys.platform == 'darwin' else '.so')
    lib = os.path.join(here, name)
    if not os.access(here, os.W_OK) and not os.path.exists(lib):   # read-only install: build in the cache
        cache = os.path.join(os.path.expanduser('~'), '.cache', 'feynsage')
        os.makedirs(cache, exist_ok=True)
        lib = os.path.join(cache, name)
    try:
        if not os.path.exists(lib) or os.path.getmtime(lib) < os.path.getmtime(src):
            subprocess.run([os.environ.get('CC', 'cc'), '-std=c99', '-O3', '-shared', '-fPIC', '-o', lib, src],
                           check=True, capture_output=True)
        L = ctypes.CDLL(lib)
        L.elim_trim.restype = ctypes.c_int64
        L.elim_trim.argtypes = [ctypes.c_int32, ctypes.c_int32, ctypes.c_void_p, ctypes.c_void_p,
                                ctypes.c_void_p, ctypes.c_int32, ctypes.c_void_p, ctypes.c_uint64,
                                ctypes.c_int32, ctypes.c_void_p, ctypes.c_void_p]
        _LIB = L
    except Exception:
        _LIB = False
    return _LIB or None



# ---------------------------------------------------------------------- Thiele interpolation over GF(p)
class Thiele:
    """Incremental Thiele (continued-fraction) interpolation of a rational function."""

    def __init__(self, F):
        self.F, self.xs, self.a = F, [], []

    def value(self, x):
        v = self.a[-1]
        for k in range(len(self.a) - 2, -1, -1):
            if v == 0:
                return None
            v = self.a[k] + (x - self.xs[k]) / v
        return v

    def add(self, x, y):
        r = y
        for xk, ak in zip(self.xs, self.a):
            if r == ak:
                return False
            r = (x - xk) / (r - ak)
        self.xs.append(x); self.a.append(r)
        return True

    def rational(self, Rx):
        """(numerator, denominator) in Rx, normalised so the lowest nonzero coefficient of
        the denominator is 1 (a unique form, the same for every prime)."""
        x = Rx.gen()
        num, den = Rx(self.a[-1]), Rx(1)
        for k in range(len(self.a) - 2, -1, -1):
            num, den = self.a[k] * num + (x - self.xs[k]) * den, num
        g = num.gcd(den)
        num, den = num // g, den // g
        c = den[min(den.exponents())]
        return num / c, den / c


class _Recon:
    """Reconstruct many univariate rational functions from shared samples."""

    def __init__(self, F, agree=2):
        self.F, self.agree = F, agree
        self.T, self.ok = {}, {}

    def feed(self, x, vals):
        for key, y in vals.items():
            T = self.T.setdefault(key, Thiele(self.F))
            if not T.a:
                T.add(x, y); self.ok[key] = 0
                continue
            if self.ok.get(key, 0) >= self.agree:
                continue
            if T.value(x) == y:
                self.ok[key] += 1
            else:
                self.ok[key] = 0
                T.add(x, y)

    def done(self, keys):
        return all(self.ok.get(k, 0) >= self.agree for k in keys)


def _univariate(mapper, F, keys_hint=None, maxpts=400, batch=1):
    """mapper([points]) -> [{key: value}].  Rebuild all keys as rational functions of x.
    Points are asked for in batches (for parallel sampling) and fed in order."""
    rec = _Recon(F)
    keys = set(keys_hint or [])
    n = 0
    while n < maxpts:
        xs = [F.random_element() for _ in range(batch)]
        for x, vals in zip(xs, mapper([[x] for x in xs])):
            if vals is None:
                continue
            keys.update(vals)
            full = {k: vals.get(k, F(0)) for k in keys}
            rec.feed(x, full)
            n += 1
            if rec.done(keys):
                Rx = PolynomialRing(F, 'x')
                return {k: rec.T[k].rational(Rx) for k in keys}
    raise RuntimeError("reconstruction did not converge")


_POOL_SYS = None


def _pool_sample(args):
    targets, p, point = args
    return _POOL_SYS.sample(targets, p, point)


# ---------------------------------------------------------------------- driver
def _primes(start=2**62):
    p = start
    while True:
        p = previous_prime(p)
        yield p


def _auto_relations(once, reducer, *args, **kw):
    """Run a reduction; with numerator_relations="auto", repeat it with the momentum-shift
    relations if a master with numerators sits in a sector where they could matter."""
    out = once(reducer, *args, **kw)
    if getattr(reducer, 'numerator_relations', False) == "auto" and not reducer._use_relations:
        masters = {m for row in out.values() for m in row}
        if reducer.relations_needed(masters):
            if kw.get('verbose'):
                print("masters with numerators in symmetric sectors: again with momentum-shift relations")
            reducer._use_relations = True
            out = once(reducer, *args, **kw)
    return out


def reduce_ff(reducer, targets, rmax, smax=0, point=None, verbose=False, nproc=1):
    r"""
    Reduce the targets with finite fields (see _reduce_ff_once); with the reducer's
    numerator_relations="auto" the reduction is repeated with the momentum-shift relations
    when a master with numerators sits in a sector where they could matter.
    """
    return _auto_relations(_reduce_ff_once, reducer, targets, rmax, smax=smax, point=point, verbose=verbose, nproc=nproc)


def _reduce_ff_once(reducer, targets, rmax, smax=0, point=None, verbose=False, nproc=1):
    r"""
    Reduce the targets with finite fields.  The reducer's kinematics may have the
    variable d and at most one invariant.  Returns {target: {master: coefficient}}
    with coefficients in the fraction field of the family's polynomial ring.
    nproc > 1 evaluates the sample points in parallel worker processes (fork); the
    first probe, which trims the system, runs before the workers start.
    """
    import time
    t0 = time.time()
    for t in targets:
        if sum(x - 1 for x in t if x > 0) > rmax or sum(-x for x in t if x < 0) > smax:
            raise ValueError("target %s needs more dots or numerator powers than rmax=%d, smax=%d"
                             % (tuple(t), rmax, smax))
    sysm = IBPSystem(reducer, rmax, smax)
    R = sysm.R
    K = R.fraction_field()
    nv = sysm.nvars
    if verbose:
        print("system: %d equations, %d integrals (%.1fs)" % (sysm.nrows, len(sysm.cols), time.time() - t0))

    residues = []          # per prime: {key: (num, den) over GF(p)[vars]}
    results_prev = None
    # univariate: 62-bit primes (plain Python ints); two variables: Singular gcd needs p < 2^29
    pool = None
    if nproc and nproc > 1:
        import multiprocessing as mproc
        global _POOL_SYS
        p0 = next(_primes(2**62 if nv == 1 else 2**29))
        sysm.sample(targets, p0, [GF(p0).random_element() for _ in range(nv)])   # trim first
        _POOL_SYS = sysm
        pool = mproc.get_context('fork').Pool(nproc)
    for p in _primes(2**62 if nv == 1 else 2**29):
        F = GF(p)
        Rp = PolynomialRing(F, [str(g) for g in R.gens()])
        if pool is not None:
            mapper = lambda pts, p=p: pool.map(_pool_sample, [(targets, p, pt) for pt in pts])
        else:
            mapper = lambda pts, p=p: [sysm.sample(targets, p, pt) for pt in pts]
        if nv == 1:
            fx = _univariate(mapper, F, batch=max(1, nproc or 1))
            recon = {k: (Rp(n.change_ring(F)(Rp.gen(0))), Rp(d.change_ring(F)(Rp.gen(0)))) for k, (n, d) in fx.items()}
        else:
            recon = _bivariate(sysm, targets, p, F, Rp, mapper)
        residues.append((p, recon))
        res = _lift(residues, R, K)
        if res is not None and res == results_prev:
            break
        results_prev = res
    out = {tuple(t): {} for t in targets}
    for (t, m), val in res.items():
        if val != 0:
            out.setdefault(t, {})[m] = val
    if pool is not None:
        pool.close()
        pool.join()
    if verbose:
        print("primes used: %d, total %.1fs" % (len(residues), time.time() - t0))
    return out


def _bivariate(sysm, targets, p, F, Rp, mapper=None):
    """Reconstruct f(d, y): in d at fixed y, then every normalised coefficient in y."""
    x0, y0 = Rp.gens()
    # 1. degrees at one random y
    y1 = F.random_element()
    mapper = mapper or (lambda pts: [sysm.sample(targets, p, pt) for pt in pts])
    base = _univariate(lambda pts: mapper([[x[0], y1] for x in pts]), F)
    shape = {k: (n.degree(), d.degree()) for k, (n, d) in base.items()}
    npts = max(a + b for a, b in shape.values()) + 4
    # 2. for each new y, rebuild in d with a fixed number of points, feed coefficients to Thiele in y
    coef_rec = _Recon(F)
    coef_keys = set()
    Rx = PolynomialRing(F, 'x')
    for _ in range(400):
        y = F.random_element()
        xs = [F.random_element() for _ in range(npts + 2)]
        samples = mapper([[x, y] for x in xs])
        if any(s is None for s in samples):
            continue
        cvals = {}
        good = True
        for k in shape:
            T = Thiele(F)
            pts = [(x, s.get(k, F(0))) for x, s in zip(xs, samples)]
            for x, v in pts[:-2]:
                T.add(x, v)
            if any(T.value(x) != v for x, v in pts[-2:]):
                good = False; break
            n, d = T.rational(Rx)
            for j in range(shape[k][0] + 1):
                cvals[(k, 'n', j)] = n[j]
            for j in range(shape[k][1] + 1):
                cvals[(k, 'd', j)] = d[j]
        if not good:
            continue
        coef_keys.update(cvals)
        coef_rec.feed(y, cvals)
        if coef_rec.done(coef_keys):
            break
    Ry = PolynomialRing(F, 'y')
    out = {}
    for k, (dn, dd) in shape.items():
        parts = {}
        for kind, deg in (('n', dn), ('d', dd)):
            parts[kind] = [coef_rec.T[(k, kind, j)].rational(Ry) for j in range(deg + 1)]
        # common denominator in y
        dens = [b for (a, b) in parts['n'] + parts['d']]
        L = lcm(dens)
        num = sum((a * (L // b))(x0.parent().gen(1)) * x0**j for j, (a, b) in enumerate(parts['n']))
        den = sum((a * (L // b))(x0.parent().gen(1)) * x0**j for j, (a, b) in enumerate(parts['d']))
        out[k] = (Rp(num), Rp(den))
    return out


def _lift(residues, R, K):
    """Combine residues from several primes into rational functions over QQ."""
    keys = residues[0][1].keys()
    out = {}
    for k in keys:
        forms = []
        for p, rec in residues:
            n, d = rec[k]
            g = n.gcd(d)
            n, d = n // g, d // g
            lc = d.lc()                      # normalise: leading coefficient of the denominator = 1
            forms.append((p, n / lc, d / lc))
        mons_n = set().union(*[set(n.dict()) for _, n, _ in forms])
        mons_d = set().union(*[set(d.dict()) for _, _, d in forms])
        M = prod(p for p, _, _ in forms)

        def lift(mon, which):
            res = [int(f[which].dict().get(mon, 0)) for f in forms]
            ps = [f[0] for f in forms]
            c = crt(res, ps)
            try:
                return QQ(ZZ(c).rational_reconstruction(M))
            except (ArithmeticError, ValueError):
                return None
        gens = R.gens()
        num, den = R(0), R(0)
        for mon in mons_n:
            c = lift(mon, 1)
            if c is None:
                return None
            num += c * prod(g**e for g, e in zip(gens, mon if hasattr(mon, '__iter__') else (mon,)))
        for mon in mons_d:
            c = lift(mon, 2)
            if c is None:
                return None
            den += c * prod(g**e for g, e in zip(gens, mon if hasattr(mon, '__iter__') else (mon,)))
        out[k] = K(num) / K(den)
    return out


def reduce_exact_trimmed(reducer, targets, rmax, smax=0, verbose=False):
    r"""
    Exact Laporta reduction on the trimmed system (see _reduce_exact_trimmed_once), with the
    same automatic momentum-shift relations as reduce_ff.
    """
    return _auto_relations(_reduce_exact_trimmed_once, reducer, targets, rmax, smax=smax, verbose=verbose)


def _reduce_exact_trimmed_once(reducer, targets, rmax, smax=0, verbose=False):
    r"""
    Exact Laporta reduction (rational functions in d and the invariants, no sampling) on the
    trimmed system: one elimination modulo a large prime at a random point finds the
    equations the targets need (as in reduce_ff); only those are then eliminated exactly.
    Same answer as reducer.run(rmax, smax) for these targets, in a fraction of the time.
    Returns {target: {master: coefficient}} like reduce_ff.
    """
    import time
    from .laporta import Reducer
    t0 = time.time()
    sysm = IBPSystem(reducer, rmax, smax)
    nv = sysm.nvars
    p = next(_primes(2**62 if nv == 1 else 2**29))
    sysm.sample([tuple(t) for t in targets], p, [GF(p).random_element() for _ in range(nv)])
    rows = sysm._rows_int
    if verbose:
        print("system: %d equations, %d kept for the targets (%.1fs)" % (sysm.nrows, len(rows), time.time() - t0))
    R = sysm.R
    gens = R.gens()
    exact = Reducer.__new__(Reducer)
    # the same family, symmetries and caches (so both pick the same representatives), no rules yet
    exact.__dict__.update(reducer.__dict__)
    exact.__dict__.update({'_zero': dict(reducer._zero), 'rules': {}, 'seen': set()})
    for row in rows:
        ident = {}
        for col, terms in row:
            poly = sum((c * prod((g ** e for g, e in zip(gens, exps)), R(1)) for c, exps in terms), R(0))
            if poly:
                key = sysm.cols[col]
                ident[key] = ident.get(key, exact.K(0)) + exact.K(poly)
        exact.add_identity(ident)
    out = {tuple(t): exact.reduce(t) for t in targets}
    if verbose:
        print("exact elimination of the trimmed system: %.1fs" % (time.time() - t0))
    return out
