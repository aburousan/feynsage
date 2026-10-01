r"""
Integral families: propagators, scalar products, graph polynomials, zero sectors and
integration-by-parts identities.

A family is fixed by
- the loop momenta ``loops`` (names),
- the kinematics (external basis, invariants, scalar-product rules),
- a list of propagators, each a pair (momentum, mass^2).  The momentum is a dict
  over loop and external names; mass^2 is an expression in the invariants.

Euclidean propagators are D = k^2 + m^2, Minkowski ones D = k^2 - m^2.
The integral with indices a = (a_1,...,a_t) is  Int prod_r [dl_r]  prod_i D_i^(-a_i);
negative indices are numerators.

The family must be *complete*: the t propagators must be linearly independent
functions of the N_sp = L(L+1)/2 + L E scalar products that involve loop momenta,
with t = N_sp.  Irreducible numerators are added as extra propagators that are
only ever used with index <= 0 (this is the convention of LiteRed and Kira).
"""
from itertools import combinations
from sage.all import QQ, PolynomialRing, matrix, vector, prod

from .momenta import add


class IntegralFamily:

    def __init__(self, name, loops, kinematics, propagators):
        self.name = name
        self.loops = list(loops)
        self.kin = kinematics
        self.props = [(dict(k), kinematics.R(m2) if not isinstance(m2, str) else kinematics.R(m2))
                      for k, m2 in propagators]
        self.t = len(self.props)
        L, E = len(self.loops), len(self.kin.externals)
        # scalar products with at least one loop momentum
        self.sp_pairs = [(self.loops[i], self.loops[j]) for i in range(L) for j in range(i, L)]
        self.sp_pairs += [(l, p) for l in self.loops for p in self.kin.externals]
        self.nsp = len(self.sp_pairs)
        # U and F need nothing more; IBP needs the propagators to be a basis of the scalar products
        self.complete = (self.nsp == self.t)
        if self.complete:
            self._build_sp_map()

    # ------------------------------------------------------------------ scalar products
    def _sp_index(self, a, b):
        for n, (x, y) in enumerate(self.sp_pairs):
            if (x, y) == (a, b) or (x, y) == (b, a):
                return n
        return None

    def dot(self, u, v):
        """u.v for two momenta (dicts): returns (vector of coefficients of the loop scalar
        products, external part in the kinematics ring)."""
        coeffs = [self.kin.R(0)] * self.nsp
        ext = self.kin.R(0)
        for a, ca in u.items():
            for b, cb in v.items():
                c = ca * cb
                n = self._sp_index(a, b)
                if n is not None:
                    coeffs[n] += c
                else:
                    ext += c * self.kin.ext_sp(a, b)
        return coeffs, ext

    def _build_sp_map(self):
        """Write every propagator as  D_i = sum_j A_ij sp_j + c_i  and invert."""
        sign = 1 if self.kin.euclidean else -1
        A, c = [], []
        for k, m2 in self.props:
            co, ext = self.dot(k, k)
            A.append(co)
            c.append(ext + sign * m2)
        self.A = matrix(self.kin.K, A)
        self.c = vector(self.kin.K, c)
        if self.A.det() == 0:
            raise ValueError("propagators do not form a basis of the scalar products")
        self.Ainv = self.A.inverse()
        # sp_j = sum_k Ainv[j,k] (D_k - c_k)

    def sp_in_dens(self, j):
        """Scalar product number j as (coefficients of D_k, constant)."""
        row = self.Ainv.row(j)
        return list(row), -(row * self.c)

    # ------------------------------------------------------------------ graph polynomials
    def UF(self, sector=None, names='x'):
        r"""
        Symanzik polynomials of the sector (tuple of 0/1, default: all lines), by
        completing the square in  sum_i x_i D_i = l.M.l - 2 Q.l + J:
        U = det M,  F = U (J - Q.M^{-1}.Q)   (Euclidean)   or   F = -U (J - Q.M^{-1}.Q)
        (Minkowski), so that in both cases F = F_0 + U sum x_i m_i^2 with the sign
        conventions of the parametric formula.
        """
        if sector is None:
            sector = (1,) * self.t
        R = PolynomialRing(self.kin.R, ['%s%d' % (names, i + 1) for i in range(self.t)])
        xs = R.gens()
        L = len(self.loops)
        M = [[R(0)] * L for _ in range(L)]
        Q = [dict() for _ in range(L)]       # Q_r as momentum dict with coefficients in R
        J = R(0)
        sign = 1 if self.kin.euclidean else -1
        for i, (k, m2) in enumerate(self.props):
            if not sector[i]:
                continue
            x = xs[i]
            c = [k.get(l, 0) for l in self.loops]
            q = {a: v for a, v in k.items() if a not in self.loops}
            for r in range(L):
                for s in range(L):
                    M[r][s] += x * c[r] * c[s]
                for a, v in q.items():
                    Q[r][a] = Q[r].get(a, R(0)) - x * c[r] * v
            J += x * (self.kin.square_external(q) + sign * m2)
        Mm = matrix(R, M)
        U = Mm.det()
        adj = Mm.adjugate()
        QMQ = R(0)
        for r in range(L):
            for s in range(L):
                for a, va in Q[r].items():
                    for b, vb in Q[s].items():
                        QMQ += adj[r, s] * va * vb * self.kin.ext_sp(a, b)
        F = U * J - QMQ           # = U (J - Q M^{-1} Q), a polynomial
        if not self.kin.euclidean:
            F = -F
        return U, F

    def is_zero_sector(self, sector):
        r"""
        Lee's criterion (R. N. Lee, 2013): the sector integrals vanish (scaleless) iff
        there is a vector k with  sum_i k_i x_i d(U+F)/dx_i = U+F,  i.e. iff the linear
        system  <k, e> = 1  has a solution for every exponent vector e of U+F.
        """
        if sum(sector) == 0:
            return True
        U, F = self.UF(sector)
        G = U + F
        if G == 0:
            return True
        idx = [i for i in range(self.t) if sector[i]]
        exps = [e if hasattr(e, '__iter__') else (e,) for e in G.exponents()]   # one line: ints
        rows = [[e[i] for i in idx] for e in exps]
        Mx = matrix(QQ, rows)
        rhs = vector(QQ, [1] * len(rows))
        try:
            Mx.solve_right(rhs)
            return True
        except ValueError:
            return False

    # ------------------------------------------------------------------ IBP identities
    def ibp(self, a):
        r"""
        All L (L + E) IBP identities for the index vector a: the integral of
        d/dl_r . ( v  prod D^{-a} ) with v in {loops, externals} vanishes.
        Each identity is a dict {index tuple: coefficient in the fraction field}.
        """
        if not self.complete:
            raise ValueError("family has %d propagators but %d scalar products: add irreducible numerators"
                             % (self.t, self.nsp))
        K = self.kin.K
        d = K(self.kin.R.gen(0))
        out = []
        vecs = self.loops + self.kin.externals
        for r, lr in enumerate(self.loops):
            for vname in vecs:
                v = {vname: QQ(1)}
                ident = {}

                def put(idx, c):
                    idx = tuple(idx)
                    ident[idx] = ident.get(idx, K(0)) + c
                    if ident[idx] == 0:
                        del ident[idx]
                if vname == lr:
                    put(a, d)
                for i, (k, m2) in enumerate(self.props):
                    ai = a[i]
                    cir = k.get(lr, 0)
                    if ai == 0 or cir == 0:
                        continue
                    # v . dD_i/dl_r = 2 c_ir (v . k_i)
                    co, ext = self.dot(v, k)
                    base = list(a); base[i] += 1
                    pref = -ai * 2 * cir
                    put(base, pref * ext)
                    for j, cj in enumerate(co):
                        if cj == 0:
                            continue
                        dens, const = self.sp_in_dens(j)
                        put(base, pref * cj * const)
                        for kk, ck in enumerate(dens):
                            if ck == 0:
                                continue
                            b = list(base); b[kk] -= 1
                            put(b, pref * cj * ck)
                if ident:
                    out.append(ident)
        return out
