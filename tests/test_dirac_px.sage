# Fermion lines, loop integrals with Dirac matrices, projectors and tensor projections, against
# Package-X 2.1.1 (its documented FermionLineExpand examples and LoopRefine numbers).
#     sage tests/test_dirac_px.sage
import os, sys, time
sys.path.insert(0, os.path.abspath('.'))
from feynsage import *
from feynsage.dirac import (line_expand, loop_line, loop_matrix, projector, form_factor, spur, trace,
                            transverse, longitudinal, chisholm, _vertex_basis)
from feynsage.pv import eps, mu, DiscB

t_start = time.time()
fails = []


def same(name, flag):
    print("%-58s %s" % (name, "ok" if flag else "FAIL"))
    if not flag:
        fails.append(name)


def close(name, a, b, tol=1e-12):
    a, b = CC(a), CC(b)
    err = abs(a - b) / max(abs(b), 1e-300)
    same("%s (rel %.0e)" % (name, err), err < tol)


def is_(r, expected):
    """r == expected, both {(structure, chiral): coefficient}."""
    keys = set(r) | set(expected)
    return all(bool((SR(r.get(k, 0)) - SR(expected.get(k, 0))).expand() == 0) for k in keys)


m, M, P = var('m M P')
mu_, nu_ = ('i', 'mu'), ('i', 'nu')
# documented FermionLineExpand examples
same("<u(p2)| g^mu (p2/ + m) g^nu |v(p1)> = 2 p2^mu <g^nu>",
     is_(line_expand(("u", "p2", "m"), ["mu", "p2 + m", "nu"], ("v", "p1", "m")), {((nu_,), '1'): 2 * SR.var('p2__mu')}))
same("<u(p2)| p1/ |u(p1)> = m <1>", is_(line_expand(("u", "p2", "m"), ["p1"], ("u", "p1", "m")), {((), '1'): m}))
same("<u(p2)| g^mu p2/ g^nu |u(p1)> (sigma, g, metric)",
     is_(line_expand(("u", "p2", "m"), ["mu", "p2", "nu"], ("u", "p1", "m")),
         {(('sigma', mu_, nu_), '1'): I * m, ((nu_,), '1'): 2 * SR.var('p2__mu'), ((), '1'): -SR.var('g__mu__nu') * m}))
same("<v(P)| g^mu P/ g^nu |u(p)> (v spinor signs)",
     is_(line_expand(("v", "P", "M"), ["mu", "P", "nu"], ("u", "p", "m")),
         {(('sigma', mu_, nu_), '1'): -I * M, ((nu_,), '1'): 2 * SR.var('P__mu'), ((), '1'): M * SR.var('g__mu__nu')}))
same("PL inside: <u(p2)| g^mu PL (p2/ + m) g^nu PL |v(p1)>",
     is_(line_expand(("u", "p2", "m"), ["mu", "PL", "p2 + m", "nu", "PL"], ("v", "p1", "m")),
         {(('sigma', mu_, nu_), 'L'): I * m, ((nu_,), 'L'): 2 * SR.var('p2__mu'), ((), 'L'): -SR.var('g__mu__nu') * m}))
same("Gordon: (p1 + p2)^mu <u u> = 2m <g^mu> - i <sigma^{mu q}>",
     is_(line_expand(("u", "p2", "m"), ["p1__mu + p2__mu"], ("u", "p1", "m")),
         {((mu_,), '1'): 2 * m, (('sigma', mu_, ('p', 'p2')), '1'): -I, (('sigma', mu_, ('p', 'p1')), '1'): I}))
same("g^mu g^nu g_mu = (2 - d) g^nu",
     is_(line_expand(("u", "p2", "m"), ["mu", "nu", "mu"], ("u", "p1", "m"), gordon=False), {((nu_,), '1'): 2 - SR.var('d')}))

# QED vertex (Package-X tutorial): structures and form factors
kin = {"p1^2": "m^2", "p2^2": "m^2", "p1.p2": "m^2 - q2/2"}
vtx = ["rho", "p2 - l + m", "mu", "p1 - l + m", "rho"]
props = (["l", "0"], ["l - p2", "m"], ["l - p1", "m"])
r = loop_line(("u", "p2", "m"), vtx, ("u", "p1", "m"), *props, kin=kin).coefficients()
num = lambda e: SR(e).subs(m=1, q2=-QQ(3)/2, mu=1).series(eps, 1).truncate().expand()
g = num(r[((mu_,), '1')])
close("vertex g^mu, 1/eps (Package-X)", g.coefficient(eps, -1).n(), 3.8241218773715148471)
close("vertex g^mu, finite (Package-X)", g.coefficient(eps, 0).n(), 4.4497712184915136478)
close("vertex sigma^{mu p2} (Package-X)", num(r[(('sigma', mu_, ('p', 'p2')), '1')]).n(), CC(0, 0.80689196496328995632))
for name, val in (("F2", 1.6137839299265799126), ("F3", 0)):
    v = num(form_factor(vtx, projector(name, 'mu', 'p1', m, 'p2', m), *props, kin=kin).coefficients().get('1', 0))
    close("vertex %s by projector (Package-X)" % name, v.n(), val) if val else same("vertex F3 = 0 (Package-X)", bool(v == 0))

# self-energy and vacuum polarisation
se = loop_matrix(["mu", "p - l + m", "mu"], ["l", "0"], ["l - p", "m"], kin={"p^2": "s"}).coefficients()
v = SR(se[(((('p', 'p'),)), '1')]).subs(s=-2, m=1, mu=1).series(eps, 1).truncate().expand()
close("self-energy p/ finite (Package-X)", v.coefficient(eps, 0).n(), 0.32395921650108226855)
v = SR(se[((), '1')]).subs(s=-2, m=1, mu=1).series(eps, 1).truncate().expand()
close("self-energy 1 finite (Package-X)", v.coefficient(eps, 0).n(), -0.59167373200865814837)
s = var('s')
tr = form.dirac_trace(['nu', 'l + m', 'mu', 'l - q + m'], vectors=['l', 'q'], dim='D')
vp = loop(str(form.to_loop(tr)).replace('D', 'd'), ["l", "m"], ["l - q", "m"], kin={"q^2": "s"})
pxT = 4*(12*m**2 + 5*s)/9 + 4*(2*m**2 + s)*DiscB(s, m, m)/3 + 4*s*(1/eps + log(mu**2/m**2))/3
close("vacuum polarisation, transverse (Package-X)", (transverse(vp, 'q', 's') - pxT).subs(s=-3, m=2, mu=5).n() + 1, 1)
same("vacuum polarisation, longitudinal = 0 (Ward identity)", bool(longitudinal(vp, 'q', 's').expand() == 0))

# projectors reproduce their own structures (unequal masses), Chisholm on explicit matrices
q2 = var('q2')
sp = {('p1', 'p1'): m**2, ('p2', 'p2'): M**2, ('p1', 'p2'): (m**2 + M**2 - q2)/2}
B = _vertex_basis('mu', 'p1', m, 'p2', M, False)
pt = {m: QQ(3)/2, M: QQ(5)/7, q2: -QQ(11)/3}
row = [spur(B[j], projector('F2', 'mu', 'p1', m, 'p2', M), ['p1', 'p2'], sp).subs({SR.var('fs_qsq'): q2}).subs(pt)
       for j in ('F1', 'F2', 'F3')]
same("Tr[B_j P_F2] = (0, 1, 0), unequal masses", [x.simplify_rational() for x in row] == [0, 1, 0])
c = chisholm({((mu_, nu_, ('i', 'rho')), '1'): SR(1)})
same("Chisholm: g^mu g^nu g^rho -> metric terms + i eps g5",
     c.get(((('i', 'rho'),), '1')) == SR.var('g__mu__nu') and c.get((('eps', mu_, nu_, ('i', 'rho')), '5')) == I)

print("time %.1f s" % (time.time() - t_start))
print("ALL DIRAC CHECKS PASS" if not fails else "FAILED: %s" % fails)
