import sys, time; sys.path.insert(0, '.')
from feynsage import *
# the only scale is q^2: set it to 1 (as Kira's symbol_to_replace_by_one); the q^2 dependence follows from dimensions
kin = Kinematics(['q'], [], {('q','q'): 1}, euclidean=False)
R = kin.R; d = R.gen(0); q2 = 1; K = kin.K
props = [(mom(l1=1),0), (mom(l1=1,q=1),0), (mom(l1=1,l2=1),0), (mom(l1=1,l2=1,q=1),0), (mom(l2=1),0)]
fam = IntegralFamily('kite', ['l1','l2'], kin, props)
# symmetries of this family (as in the course notebook): (1,2,3,4,5) -> (3,4,1,2,5) and (2,1,4,3,5)
red = Reducer(fam, symmetries=[(2,3,0,1,4), (1,0,3,2,4)])   # generators; the group is closed automatically
t0 = time.time(); red.run(rmax=4, smax=2); print("time %.1fs" % (time.time()-t0))
print("masters:", red.masters())
# reference values (LiteRed / course notebook, same family)
ref_kite = {(0,1,1,0,1): 2*(3*d-10)*(3*d-8)/((d-4)**2*q2**2), (1,1,1,1,0): -2*(d-3)/((d-4)*q2)}
got = red.reduce((1,1,1,1,1))
print("kite F(1,1,1,1,1):", {k: K(v).factor() for k, v in got.items()})
ok = all(K(got.get(red.canon(m),0)) - K(c) == 0 for m, c in ref_kite.items())
ref2 = {(0,1,1,0,1): -9*(d-7)*(3*d-16)*(3*d-14)*(3*d-10)*(3*d-8)*(d**2-11*d+20)/((d-10)*(d-8)**2*(d-4)*q2**6),
        (1,1,1,1,0): -4*(d-9)*(d-6)*(d-5)*(d-3)/((d-8)*q2**5)}
got2 = red.reduce((2,2,1,2,2))
ok2 = all(K(got2.get(red.canon(m),0)) - K(c) == 0 for m, c in ref2.items())
print("F(1,1,1,1,1) agrees:", ok, "   F(2,2,1,2,2) agrees:", ok2)
