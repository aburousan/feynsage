# finite-field reduction against the exact reducer and the Kira / LiteRed references
import sys, time; sys.path.insert(0, '.')
from feynsage import *

# --- kite, univariate (q^2 = 1)
kin = Kinematics(['q'], [], {('q','q'): 1}, euclidean=False)
d = kin.R.gen(0); K = kin.K
props = [(mom(l1=1),0), (mom(l1=1,q=1),0), (mom(l1=1,l2=1),0), (mom(l1=1,l2=1,q=1),0), (mom(l2=1),0)]
fam = IntegralFamily('kite', ['l1','l2'], kin, props)
red = Reducer(fam, symmetries=[(2,3,0,1,4), (1,0,3,2,4)])
t0 = time.time()
res = reduce_ff(red, [(1,1,1,1,1), (2,2,1,2,2)], rmax=4, smax=2, verbose=True)
print("kite ff time %.1fs" % (time.time()-t0))
ref1 = {(0,1,1,0,1): 2*(3*d-10)*(3*d-8)/(d-4)**2, (1,1,1,1,0): -2*(d-3)/(d-4)}
ref2 = {(0,1,1,0,1): -9*(d-7)*(3*d-16)*(3*d-14)*(3*d-10)*(3*d-8)*(d**2-11*d+20)/((d-10)*(d-8)**2*(d-4)),
        (1,1,1,1,0): -4*(d-9)*(d-6)*(d-5)*(d-3)/(d-8)}
for t, ref in [((1,1,1,1,1), ref1), ((2,2,1,2,2), ref2)]:
    got = res[t]
    print(t, sorted(got), all(K(got.get(red.canon(m), 0)) - K(c) == 0 for m, c in ref.items()))

# --- bubble, bivariate (d, s)
kin = Kinematics(['p'], ['s'], {('p','p'): 's'}, euclidean=False)
d, s = kin.R.gens(); K = kin.K
fam = IntegralFamily('bub', ['l'], kin, [(mom(l=1), 1), (mom(l=1, p=-1), 1)])
red = Reducer(fam, symmetries=[(1, 0)])
t0 = time.time()
res = reduce_ff(red, [(2,1), (2,2), (3,1), (0,2)], rmax=3, verbose=True)
print("bubble ff time %.1fs" % (time.time()-t0))
kira = {
 (2,1): {(1,1): (-d+3)/(s-4), (1,0): (d-2)/(2*s-8)},
 (2,2): {(1,1): ((d**2-9*d+18)*s+4*d-12)/(s**3-8*s**2+16*s), (1,0): ((d-2)*s-2*d**2+10*d-12)/(s**3-8*s**2+16*s)},
 (0,2): {(1,0): (d-2)/2},
 (3,1): {(1,1): ((d**2-7*d+12)*s-4*d+12)/(2*s**3-16*s**2+32*s),
         (1,0): ((d**2-6*d+8)*s**2+(-8*d**2+48*d-64)*s+8*d**2-40*d+48)/(8*s**3-64*s**2+128*s)},
}
ok = True
for a, ref in kira.items():
    got = res[a]
    for m, c in ref.items():
        ok = ok and K(got.get(red.canon(m), 0)) - K(c) == 0
print("bubble ff agrees with Kira:", ok)
