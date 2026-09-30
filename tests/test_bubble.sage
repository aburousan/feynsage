import sys; sys.path.insert(0, '.')
from feynsage import *
kin = Kinematics(['p'], ['s'], {('p','p'): 's'}, euclidean=False)
R = kin.R; d, s = R.gens()
fam = IntegralFamily('bub', ['l'], kin, [(mom(l=1), 1), (mom(l=1, p=-1), 1)])
red = Reducer(fam, symmetries=[(1, 0)]).run(rmax=3)
K = kin.K
print("masters:", red.masters())
kira = {  # Kira 3.1 output (kira_targets.m), bub[1,0] = bub[0,1] by symmetry
 (2,1): {(1,1): (-d+3)/(s-4), (1,0): (d-2)/(2*s-8)},
 (2,2): {(1,1): ((d**2-9*d+18)*s+4*d-12)/(s**3-8*s**2+16*s), (1,0): ((d-2)*s-2*d**2+10*d-12)/(s**3-8*s**2+16*s)},
 (0,2): {(1,0): (d-2)/2},
 (3,1): {(1,1): ((d**2-7*d+12)*s-4*d+12)/(2*s**3-16*s**2+32*s),
         (1,0): ((d**2-6*d+8)*s**2+(-8*d**2+48*d-64)*s+8*d**2-40*d+48)/(8*s**3-64*s**2+128*s)},
}
ok = True
for a, ref in kira.items():
    got = red.reduce(a)
    for m, c in ref.items():
        diff = K(got.get(red.canon(m), 0)) - K(c)
        print(a, m, "difference:", diff)
        ok = ok and diff == 0
print("ALL AGREE WITH KIRA" if ok else "MISMATCH")
