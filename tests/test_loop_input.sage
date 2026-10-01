# Input handling of pv.loop(): the order of the propagators and Einstein contractions
# (both found in an independent review) must not change the answer.
import sys; sys.path.insert(0, '.')
from feynsage.pv import *
ok = True
def same(a, b):
    return bool((SR(a) - SR(b)).expand() == 0)
K = {"p^2": "s"}
# propagator order
ok &= same(loop("l^2", ["l", "m"], ["l+p", "m"], kin=K)["1"], loop("l^2", ["l+p", "m"], ["l", "m"], kin=K)["1"])
ok &= same(loop("l.p", ["l", "m"], ["l+p", "m"], kin=K)["1"], loop("l.p", ["l+p", "m"], ["l", "m"], kin=K)["1"])
# Einstein summation
ok &= same(loop("l^mu l^mu", ["l", "m"], ["l+p", "m"], kin=K)["1"], loop("l^2", ["l", "m"], ["l+p", "m"], kin=K)["1"])
ok &= same(loop("l^mu p^mu", ["l", "m"], ["l+p", "m"], kin=K)["1"], loop("l.p", ["l", "m"], ["l+p", "m"], kin=K)["1"])
ok &= same(loop("g(mu,nu)*p^mu*p^nu", ["l", "m"], kin=K)["1"], var('s') * A0(var('m')))
ok &= same(loop("g(mu,mu)", ["l", "m"])["1"], (4 * A0(var('m')) - 2 * var('m')^2))
ok &= same(loop("g(mu,nu)*l^nu", ["l", "m1"], ["l+p", "m2"], kin=K)["p^mu"], loop("l^mu", ["l", "m1"], ["l+p", "m2"], kin=K)["p^mu"])
print("loop() input handling:", "OK" if ok else "FAILED")

# a squared scalar product written by FORM: p.pp^2 means (p.pp)^2, and pp.l = l.pp
r1 = loop("l.pp*p.pp^2 + pp.l", ["l", "1"], ["l + p", "1"], kin={"p^2": 1, "pp^2": 1, "p.pp": "c"})
r2 = loop("(c^2 + 1)*l.pp", ["l", "1"], ["l + p", "1"], kin={"p^2": 1, "pp^2": 1, "p.pp": "c"})
assert bool((SR(r1["1"]) - SR(r2["1"])).expand() == 0), "p.pp^2 parsing"
print("FORM-style squares of scalar products: OK")
