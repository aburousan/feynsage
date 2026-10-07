r"""
Derivatives of Feynman integrals with respect to the kinematic invariants, and differential equations
for master integrals.

    derivative(fam, "J(1,1)", "pp")         dJ(1,1)/dpp as a combination of integrals of the family
    diff_reduce(fam, "J(1,1)", "pp")        the same, reduced to master integrals (IBP)
    differential_equation(fam, ["J(0,1)", "J(1,1)"], "pp")    the matrix A in  d m/dpp = A m

How.  Let s_kl = p_k.p_l be the scalar products of the independent external momenta.  The operator
p_i . d/dp_j acts on the integrand like the IBP operators (it differentiates the propagators), and on a
function of the s_kl it gives  sum_{k<=l} (delta_jk s_il + delta_jl s_ik) df/ds_kl.  These E^2 equations
give every df/ds_kl as a combination of integrals.  An invariant x enters through the given values of
the s_kl, so  df/dx = sum df/ds_kl  ds_kl/dx  (the chain rule on the surface fixed by the kinematics).
"""
from sage.all import QQ, SR, matrix, vector

from .easy import _target, ibp_reduce


def _put(dct, idx, c):
    idx = tuple(idx)
    dct[idx] = dct.get(idx, 0) + c
    if dct[idx] == 0:
        del dct[idx]


def ext_operator(fam, a, pi, pj):
    """p_i . d/dp_j applied to the integral with powers a: {powers: coefficient}."""
    out = {}
    v = {pi: QQ(1)}
    for i, (k, m2) in enumerate(fam.props):
        ai, eij = a[i], k.get(pj, 0)
        if ai == 0 or eij == 0:
            continue
        co, ext = fam.dot(v, k)                    # p_i . k_i
        base = list(a); base[i] += 1
        pref = -ai * 2 * eij
        _put(out, base, fam.kin.K(pref * ext))
        for j, cj in enumerate(co):
            if cj == 0:
                continue
            dens, const = fam.sp_in_dens(j)
            _put(out, base, fam.kin.K(pref * cj * const))
            for kk, ck in enumerate(dens):
                if ck == 0:
                    continue
                b = list(base); b[kk] -= 1
                _put(out, b, fam.kin.K(pref * cj * ck))
    return out


def _sp_derivatives(fam, a):
    """{(k, l): {powers: coefficient}} = d F(a) / d(p_k.p_l) for k <= l (independent externals)."""
    E = fam.kin.externals
    K = fam.kin.K
    pairs = [(E[k], E[l]) for k in range(len(E)) for l in range(k, len(E))]
    rows, rhs = [], []
    for pi in E:
        for pj in E:
            row = []
            for (pk, pl) in pairs:
                c = K(0)
                if pj == pk:
                    c += fam.kin.ext_sp(pi, pl)
                if pj == pl:
                    c += fam.kin.ext_sp(pi, pk)
                row.append(c)
            rows.append(row)
            rhs.append(ext_operator(fam, a, pi, pj))
    M = matrix(K, rows)
    piv = M.transpose().pivots()                   # independent rows
    if len(piv) < len(pairs):
        raise ValueError("the derivatives with respect to the scalar products are not fixed by these kinematics")
    Msq = M.matrix_from_rows(piv)
    Minv = Msq.inverse()
    out = {}
    for n, pr in enumerate(pairs):
        acc = {}
        for r_idx, r in enumerate(piv):
            c = Minv[n, r_idx]
            if c == 0:
                continue
            for idx, v in rhs[r].items():
                _put(acc, idx, c * v)
        out[pr] = acc
    return out


def derivative(fam, target, x):
    """d F(a)/dx for an invariant x (a name such as "s" or "pp"): {powers: coefficient}."""
    a = _target(target)
    R, K = fam.kin.R, fam.kin.K
    xv = R(str(x))
    out = {}
    for (pk, pl), dct in _sp_derivatives(fam, a).items():
        ds = fam.kin.ext_sp(pk, pl).derivative(xv)
        if ds == 0:
            continue
        for idx, v in dct.items():
            _put(out, idx, K(ds) * v)
    # masses: d/dm2 of a propagator k^2 -+ m2 raises its power
    sign = 1 if fam.kin.euclidean else -1
    for i, (k, m2) in enumerate(fam.props):
        dm = R(m2).derivative(xv) if m2 != 0 else 0
        if dm == 0 or a[i] == 0:
            continue
        b = list(a); b[i] += 1
        _put(out, b, K(-a[i] * sign * dm))
    return out


def _to_masters(fam, dct, method="auto"):
    targets = sorted(dct)
    red = ibp_reduce(fam, targets, method=method)
    out = {}
    for t, c in dct.items():
        for m, cm in red.table[t].items():
            _put(out, m, c * fam.kin.K(cm))
    return out, red


def diff_reduce(fam, target, x, method="auto"):
    """d F(a)/dx reduced to master integrals, as a Sage expression in the family's integrals."""
    out, _ = _to_masters(fam, derivative(fam, target, x), method)
    return _as_sr(fam, out)


def _as_sr(fam, dct):
    tot = SR(0)
    for m, c in sorted(dct.items()):
        tot += SR(str(c)).factor() * fam.integral(*m)
    return tot


def differential_equation(fam, masters, x, method="auto"):
    r"""
    The matrix A (rational functions of d and the invariants) with  d m/dx = A m  for the master
    integrals m = masters (given as "J(1,1)" or tuples, in the order wanted).  The masters can be any
    basis of the integrals the reduction leaves; the change of basis is done here.
    """
    ms = [_target(m) for m in masters]
    K = fam.kin.K
    xv = fam.kin.R(str(x))
    # reduce the derivatives of the chosen masters and the chosen masters themselves
    ders = [derivative(fam, m, x) for m in ms]
    allt = sorted(set().union(*ders) | set(ms))
    red = ibp_reduce(fam, allt, method=method)
    basis = sorted({b for t in allt for b in red.table[t]})
    def row_of(dct):
        acc = {}
        for t, c in dct.items():
            for b, cb in red.table[t].items():
                _put(acc, b, K(c) * K(cb))
        return [acc.get(b, K(0)) for b in basis]
    Tm = matrix(K, [row_of({m: K(1)}) for m in ms])          # chosen masters in the reduction's basis
    if Tm.nrows() != len(basis) or Tm.det() == 0:
        raise ValueError("the integrals %s are not a basis of the masters %s" % (masters, basis))
    Dm = matrix(K, [row_of(dct) for dct in ders])            # d m/dx in the reduction's basis
    A = Dm * Tm.inverse()
    return A.apply_map(lambda c: SR(str(c)).factor())
