# Tests of topologies.py, models.py and diagrams.py (diagram generation like FeynArts + FeynCalc).
# Reference numbers: FeynArts 3.12 / FeynCalc 10.2 (run locally, 2026-10-08) and an explicit
# helicity sum with Dirac matrices, spinors and polarization vectors (no traces).
#     sage tests/test_diagrams.sage
import time
from feynsage import *
from feynsage.topologies import topologies
from feynsage.diagrams import insert_fields
from feynsage.models import SM, EL, SW, CW, MW, MZ, MH, ME, MM, MU, GS
from feynsage.qft import tensors as T
t_start = time.time()
FAIL = []
def check(name, ok):
    print(("ok   " if ok else "FAIL ") + name)
    if not ok:
        FAIL.append(name)

# ---------------------------------------------------------------- topologies vs FeynArts CreateTopologies
for (L, a, b), ex, n in [((0,1,2),(),1), ((0,2,2),(),4), ((0,2,3),(),25), ((0,2,4),(),220), ((1,1,1),(),3),
                         ((1,1,2),(),14), ((1,2,2),(),99), ((1,2,2),('tadpoles',),74),
                         ((1,2,2),('tadpoles','wf'),42), ((1,1,2),('tadpoles','wf'),4), ((2,1,1),(),29)]:
    check("topologies(%d, %d, %d, exclude=%s) = %d" % (L, a, b, ex, n), len(topologies(L, a, b, exclude=ex)) == n)

# ---------------------------------------------------------------- insertions vs FeynArts InsertFields (Particles)
sm = SM()
for (L,a,b), fi, fo, ex, n in [((0,2,2),['e','e~'],['mu','mu~'],(),4), ((0,2,2),['e','e~'],['e','e~'],(),8),
        ((0,2,2),['e','e~'],['W+','W-'],(),4), ((0,2,2),['e','e~'],['Z','H'],(),4), ((0,2,2),['u','u~'],['g','g'],(),3),
        ((0,2,2),['g','g'],['g','g'],(),4), ((0,2,2),['W+','W-'],['W+','W-'],(),7), ((0,2,3),['e','e~'],['mu','mu~','A'],(),16),
        ((0,2,2),['u','d~'],['W+','A'],(),4), ((1,1,1),['A'],['A'],(),17), ((1,1,1),['e'],['e'],(),63),
        ((1,1,2),['A'],['e','e~'],('tadpoles','wf'),8), ((1,2,2),['e','e~'],['mu','mu~'],('tadpoles','wf'),388)]:
    d = insert_fields(topologies(L, a, b, exclude=ex), fi, fo, sm)
    check("insert_fields %s -> %s, %d loop(s): %d diagrams" % (fi, fo, L, n), len(d) == n)
check("e e~ -> mu mu~ per topology [0, 4, 0, 0] as FeynArts",
      insert_fields(topologies(0, 2, 2), ['e-','e+'], ['mu-','mu+'], sm).counts() == [0, 4, 0, 0])

# ---------------------------------------------------------------- amplitudes
p1, p2, p3, p4 = momenta("p1 p2 p3 p4")
s, t, u = var('s t u')
def kin(ma, mb):
    return {dot(p1,p1): ma^2, dot(p2,p2): ma^2, dot(p3,p3): mb^2, dot(p4,p4): mb^2,
            dot(p1,p2): (s - 2*ma^2)/2, dot(p3,p4): (s - 2*mb^2)/2, dot(p1,p3): (ma^2 + mb^2 - t)/2,
            dot(p2,p4): (ma^2 + mb^2 - t)/2, dot(p1,p4): (ma^2 + mb^2 - u)/2, dot(p2,p3): (ma^2 + mb^2 - u)/2}
def point(rs, c, ma, mb):
    sv = rs^2; pin = sqrt(sv/4 - ma^2); pout = sqrt(sv/4 - mb^2)
    tv = ma^2 + mb^2 - sv/2 + 2*pin*pout*c
    return {s: sv, t: tv, u: 2*ma^2 + 2*mb^2 - sv - tv}
me, mmu, mw = QQ(511)/10^6, QQ(105658)/10^6, QQ(80379)/1000
num = {EL: sqrt(4*pi*1000/137036), ME: me, MM: mmu, MZ: QQ(911876)/10^4, MW: mw, MH: QQ(12525)/100,
       SW: sqrt(QQ(23122)/10^5), CW: sqrt(1 - QQ(23122)/10^5)}
def close(a, b, digits):
    return abs(RealField(200)(a)/RealField(200)(b) - 1) < 10^(-digits)

# e- e+ -> mu- mu+ in the SM (the FeynArts/FeynCalc notebook): gamma, Z, H, G0
d = insert_fields(topologies(0, 2, 2), ['e-','e+'], ['mu-','mu+'], sm)
m2 = (d.squared([p1, p2], [p3, p4]) / 4).subs(kin(ME, MM))
for (rs, c), ref in {(70, 1/2): '0.00472095317462194382258', (70, -1/2): '0.02159285030659615502670',
                     (30, 3/10): '0.00874837012487834168541', (91, 1/5): '64.1636708177280719319537'}.items():
    check("e e~ -> mu mu~ |M|^2 at sqrt s = %s, cos = %s = FeynCalc (20 digits)" % (rs, c),
          close(m2.subs(num).subs(point(rs, QQ(c), me, mmu)).n(digits=40), ref, 20))
m2u = (d.squared([p1, p2], [p3, p4], gauge='unitary') / 4).subs(kin(ME, MM))
check("e e~ -> mu mu~: unitary gauge = Feynman gauge when M_W = M_Z c_W",
      (m2u - m2).subs({MW: MZ*CW}).subs({SW: sqrt(1 - CW^2)}).subs(u=2*ME^2 + 2*MM^2 - s - t).simplify_rational() == 0)

# Bhabha in the SM: relative sign of the s- and t-channel diagrams
d = insert_fields(topologies(0, 2, 2), ['e-','e+'], ['e-','e+'], sm)
b = (d.squared([p1, p2], [p3, p4]) / 4).subs(kin(ME, ME)).subs(num).subs(point(70, QQ(3)/10, me, me))
check("Bhabha |M|^2 = FeynCalc (20 digits)", close(b.n(digits=40), '0.17182883675669547391741592', 20))

# e- e+ -> W+ W-: explicit helicity sum, and the gauge cancellation at high energy
d = insert_fields(topologies(0, 2, 2), ['e-','e+'], ['W+','W-'], sm)
e3, e4 = T._declare('eps3', 'vector'), T._declare('eps4', 'vector')
w = (polarization_sum(d.squared([p1, p2], [p3, p4]), [e3, e4]) / 4).subs(kin(ME, MW)).subs(num)
check("e e~ -> W+ W- = explicit helicity sum (20 digits)",
      close(w.subs(point(250, QQ(3)/10, me, mw)).n(digits=40), '0.037111744563644591303', 20))
hi = [w.subs(point(rs, QQ(3)/10, me, mw)).n() for rs in (1000, 10000, 100000)]
check("e e~ -> W+ W-: |M|^2 stays bounded (gauge cancellation)", max(hi)/min(hi) < 1.2)

# u u~ -> g g: Xianyu's solution of Peskin problem 17.3, eq. (17.23)
d = insert_fields(topologies(0, 2, 2), ['u','u~'], ['g','g'], sm, exclude_fields=['A','Z','H','G0','G+','W+'])
eg3, eg4 = T._declare('eps3', 'vector'), T._declare('eps4', 'vector')
r = d.squared([p1, p2], [p3, p4]).subs({MU: 0})
r = polarization_sum(polarization_sum(r, eg3, gauge="physical", n=p4), eg4, gauge="physical", n=p3)
r = (color_factor(r, N=3).subs(kin(0, 0)).subs(s=-t-u) / 36).simplify_rational()
check("u u~ -> g g = (32/27) g^4 [t/u + u/t - 9 (t^2 + u^2)/(4 s^2)]",
      (r - 32*GS^4/27*(t/u + u/t - 9*(t^2 + u^2)/(4*(t + u)^2))).simplify_rational() == 0)

# g g -> g g: Xianyu's eq. (17.34), (9/2) g^4 (3 - tu/s^2 - us/t^2 - st/u^2), at s = 100, t = -30, u = -70;
# physical polarization sums, the result must not depend on their reference vectors
d = insert_fields(topologies(0, 2, 2), ['g','g'], ['g','g'], sm)
pt = {s: 100, t: -30, u: -70}
numk = {k: v.subs(pt) for k, v in kin(0, 0).items()}
rg = color_factor(d.squared([p1, p2], [p3, p4]).subs(numk).subs({GS: 1}), N=3)
ge = {i: T._declare('eps%d' % i, 'vector') for i in (1, 2, 3, 4)}
def gsum(refs):
    x = rg
    for i in (1, 2, 3, 4):
        x = polarization_sum(x, ge[i], gauge="physical", n=refs[i]).subs(numk)
    return x/256
g1, g2 = gsum({1: p2, 2: p1, 3: p4, 4: p3}), gsum({1: p3, 2: p4, 3: p1, 4: p2})
check("g g -> g g = (9/2)(3 - tu/s^2 - us/t^2 - st/u^2) [Xianyu (17.34)]", (g1 - (9/2*(3 - t*u/s^2 - u*s/t^2 - s*t/u^2)).subs(pt)).simplify_rational() == 0)
check("g g -> g g does not depend on the reference vectors", (g2 - g1).simplify_rational() == 0)

print("time %.0f s" % (time.time() - t_start))
print("ALL DIAGRAM CHECKS PASS" if not FAIL else "FAILED: %s" % FAIL)
