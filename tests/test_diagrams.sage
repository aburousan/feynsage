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

# ---------------------------------------------------------------- one loop: H -> gamma gamma, g - 2 (P&S 21.1)
from feynsage.pv import finite_part, uv_part
from feynsage.models import MM, MT
from feynsage import u as u_spinor
k1, k2 = momenta("k1 k2")
daa = insert_fields(topologies(1, 1, 2, exclude=("tadpoles", "wf")), ["H"], ["A", "A"], sm,
                    exclude_fields=['e', 'mu', 'ta', 'ne', 'nm', 'nt', 'u', 'c', 'd', 's'])
check("H -> gamma gamma: 30 diagrams with W, Goldstones, ghosts, b and t (FeynArts: 30)", len(daa) == 30)
def transverse(n, ext):
    mu, nu = ext[1], ext[2]
    return n*(metric(mu, nu) - (comp(k2, mu)*comp(k1, nu) + comp(k1, mu)*comp(k2, nu))/dot(k1, k2))
Maa = sum(daa.loop_amplitude([k1 + k2], [k1, k2], transverse, kin={"k1^2": 0, "k2^2": 0, "k1.k2": "MH^2/2"}))
check("H -> gamma gamma: the 1/eps pole cancels", uv_part(Maa).simplify_full() == 0)
p1, p2 = momenta("p1 p2"); q2 = var('q2')
def F2proj(n, ext):
    return dirac_trace((slash(p1) + MM)*vertex_projector(p1, p2, MM, ext[1], "F2")*(slash(p2) + MM)*n, dim='d', gamma5_scheme='NDR-even')
vk = {"p1^2": "MM^2", "p2^2": "MM^2", "p1.p2": "MM^2 - q2/2"}
ALL = ['A','Z','W+','H','G0','G+','u+','u-','uA','uZ','g','ug','e','mu','ta','ne','nm','nt','u','c','t','d','s','b']
tv = topologies(1, 2, 1, exclude=('tadpoles', 'wf'))
at0 = lambda F: (F.subs(q2=10^-10) + F.subs(q2=-10^-10))/2
dq = insert_fields(tv, ['mu', 'A'], ['mu'], sm, exclude_fields=[f for f in ALL if f not in ('A', 'mu')])
F2q = at0(finite_part(sum(dq.loop_amplitude([p1, p2 - p1], [p2], F2proj, kin=vk)), 1)/EL).subs({EL: 1, MM: 1})
check("QED vertex: F2(0) = alpha/(2 pi) (Schwinger)", close(F2q.n(digits=60), (1/(8*pi^2)).n(digits=60), 15))
unit = EL^2/(4*sqrt(2)*SW^2*MW^2)*MM^2/(8*pi^2*sqrt(2))
heavy = {EL: 1, MM: 1, MW: 1000, SW: sqrt(QQ(23)/100), CW: sqrt(QQ(77)/100), MZ: 1000/sqrt(QQ(77)/100)}
dW = insert_fields(tv, ['mu', 'A'], ['mu'], sm, exclude_fields=[f for f in ALL if f not in ('W+', 'G+', 'nm')])
aW = [at0(finite_part(sum(dW.loop_amplitude([p1, p2 - p1], [p2], F2proj, kin=vk, xi=xi)), 1)/EL/unit).subs(heavy).n(digits=60)
      for xi in (1, 3)]
check("muon g - 2, W loops: a = (10/3) G_F m^2/(8 pi^2 sqrt 2) for m_W = 1000 m [P&S 21.1(a)]", abs(aW[0] - 10/3) < 1e-5)
check("muon g - 2, W loops: the same in R_xi gauge with xi = 3 [P&S 21.1(b)]", abs(aW[1] - aW[0]) < 1e-30)
dZ = insert_fields(tv, ['mu', 'A'], ['mu'], sm, exclude_fields=[f for f in ALL if f not in ('Z', 'G0', 'mu')])
aZ = at0(finite_part(sum(dZ.loop_amplitude([p1, p2 - p1], [p2], F2proj, kin=vk)), 1)/EL/unit).subs(heavy).n(digits=60)
w2 = QQ(23)/100
check("muon g - 2, Z loops: -(4/3 + 8/3 s^2 - 16/3 s^4) [P&S 21.1(c)]", abs(aZ + (4/3 + 8/3*w2 - 16/3*w2^2)) < 1e-4)

# ---------------------------------------------------------------- e- e+ -> W+_0 W-_0 by helicity (P&S 21.2, eq. 21.108)
from feynsage.diagrams import conjugate_amplitude
d = insert_fields(topologies(0, 2, 2), ['e-','e+'], ['W+','W-'], sm)
p1, p2, p3, p4 = momenta("p1 p2 p3 p4")
M = d.amplitude([p1, p2], [p3, p4])
P34 = (s - 2*MW^2)/2; NL = MW*sqrt(P34^2 - MW^4)
eL3, eL4 = (P34*p3 - MW^2*p4)/NL, (P34*p4 - MW^2*p3)/NL
lng = {T._declare(n, 'vector'): v for n, v in (('eps3', eL3), ('eps3_c', eL3), ('eps4', eL4), ('eps4_c', eL4))}
Mh = chiral(M, p1, "L")
LR = SR(spin_sum(Mh*conjugate_amplitude(Mh), rules=lng)).subs(kin(ME, MW)).subs({ME: 0})
J = chiral(vbar(p2)*slash(p3 - p4)*u_spinor(p1), p1, "L")
J2 = SR(spin_sum(J*conjugate_amplitude(J))).subs(kin(0, MW))
mw_ = QQ(804)/10; mz_ = mw_/sqrt(1 - w2); nw = {EL: 1, SW: sqrt(w2), CW: sqrt(1 - w2), MW: mw_, MZ: mz_}
rs, cc = 500, -QQ(1)/2
b_ = sqrt(1 - 4*mw_^2/rs^2); sv = rs^2
br = (1/(2*w2)*(-sv/(sv - mz_^2)*(mz_^2/(2*mw_^2) + 1) + 2/b_^2 - 8*mw_^2/(sv*b_^2*(1 + b_^2 + 2*b_*cc))) + mz_^2/mw_^2*((sv/2 + mw_^2)/(sv - mz_^2)))/sv
pt = point(rs, cc, 0, mw_)
check("e-_L e+_R -> W+_0 W-_0 = Peskin (21.108) (15 digits)", close(LR.subs(nw).subs(pt).n(digits=60), (J2.subs(nw).subs(pt)*br^2).n(digits=60), 15))

# ---------------------------------------------------------------- one-loop self-energies and vertices = FeynCalc
from feynsage.qft.color import f_color, color_factor, color_indices, color_delta
from feynsage.diagrams import loop_amplitude
from feynsage.pv import loop as pv_loop
check("massless tadpole squared is scaleless: int 1/(l^2)^2 = 0",
      all(SR(v).is_zero() for v in pv_loop("1", ["l", "0"], ["l", "0"]).coefficients().values()))
P2v, A22, A33, A23 = var('P2v A22 A33 A23')
pp = momenta("pp"); q2_, q3_ = momenta("q2_ q3_")
se_t = topologies(1, 1, 1, exclude=('tadpoles', 'wf'))
only = lambda keep: [f for f in ALL if f not in keep]
de = insert_fields(se_t, ['e'], ['e'], sm, exclude_fields=only(('A', 'e')))
ske = {"pp^2": "P2v"}
Ae = uv_part(sum(de.loop_amplitude([pp], [pp], lambda n, e: dirac_trace(slash(pp)*n, dim='d')/(4*P2v), ske, xi=3)))
Be = uv_part(sum(de.loop_amplitude([pp], [pp], lambda n, e: dirac_trace(n, dim='d')/4, ske, xi=3)))
check("electron self-energy pole, xi = 3: e^2 (xi p/ - (xi + 3) m)/(16 pi^2) [FeynCalc QED/OneLoop/El-El]",
      bool((Ae - 3*EL^2/(16*pi^2)).simplify_full() == 0) and bool((Be + 6*EL^2*ME/(16*pi^2)).simplify_full() == 0))
dgl = insert_fields(se_t, ['g'], ['g'], sm, exclude_fields=only(('g', 'ug')))
a1_, a2_, a3_ = color_indices("a1 a2 a3")
gl = [uv_part(sum(dgl.loop_amplitude([pp], [pp], lambda n, e: n*T.metric(e[0], e[1]), ske, xi=x))).subs({color_delta(a1_, a2_): 1})
      for x in (1, 3)]
check("gluon self-energy pole (gluon + ghost loops, four-gluon tadpole), xi = 1 and 3 [FeynCalc QCD/OneLoop/Gl-Gl]",
      all(bool((g_ - GS^2*(-3*P2v)*(9*x - 39)/(96*pi^2)).simplify_full() == 0) for g_, x in zip(gl, (1, 3))))
dgh = insert_fields(se_t, ['ug'], ['ug'], sm, exclude_fields=only(('g', 'ug')))
gh = uv_part(sum(dgh.loop_amplitude([pp], [pp], lambda n, e: n, ske))).subs({color_delta(a1_, a2_): 1})
check("ghost self-energy pole = -3 g^2 p^2/(32 pi^2) (open ghost line: no loop sign) [FeynCalc QCD/OneLoop/Gh-Gh]",
      bool((gh + 3*GS^2*P2v/(32*pi^2)).simplify_full() == 0))
vt = topologies(1, 1, 2, exclude=('tadpoles', 'wf'))
kv = {"q2_^2": "A22", "q3_^2": "A33", "q2_.q3_": "A23"}
dgv = insert_fields(vt, ['ug'], ['ug', 'g'], sm, exclude_fields=only(('g', 'ug')))
gv = [color_factor(uv_part(sum(dgv.loop_amplitude([q2_ + q3_], [q2_, q3_], lambda n, e: n*T.comp(q2_, e[2]), kv, xi=x)))
                   * f_color(a1_, a2_, a3_), N=3) for x in (0, 3)]
check("ghost-gluon vertex pole = -(9/4) i xi g^3 p_ghost/pi^2: zero in Landau gauge (Taylor) [Muta 2.5.142 via FeynCalc]",
      bool(SR(gv[0]).simplify_full() == 0) and bool((gv[1] + 27*I*GS^3*A22/(4*pi^2)).simplify_full() == 0))
dg3 = insert_fields(vt, ['g'], ['g', 'g'], sm, exclude_fields=only(('g', 'ug')))
ghl = [x for x in dg3 if any(f[0].startswith('ug') for f in x.fields)]
tri = color_factor(uv_part(sum(loop_amplitude(ghl, [q2_ + q3_], [q2_, q3_],
                   lambda n, e: n*T.metric(e[0], e[1])*T.comp(q2_, e[2]), kv)))*f_color(a1_, a2_, a3_), N=3)
check("three-gluon vertex, ghost triangles (three ghost vertices): pole = FeynCalc",
      bool((tri + I*9*GS^3*(2*A22 + A23)/(16*pi^2)).simplify_full() == 0))

print("time %.0f s" % (time.time() - t_start))
print("ALL DIAGRAM CHECKS PASS" if not FAIL else "FAILED: %s" % FAIL)
