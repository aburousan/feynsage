# feynsage against textbook formulas.  Every equation below was read on the printed page of
# M. E. Peskin and D. V. Schroeder, An Introduction to Quantum Field Theory (Addison-Wesley, 1995);
# the page numbers are those printed in the book.
#     sage tests/test_books.sage
import time
from feynsage import *
from feynsage.models import EL, ME, GS, MU
from feynsage.qft import tensors as T
from feynsage.qft.color import color_factor, color_indices, quark_colors, T_color, color_delta
from feynsage.diagrams import loop_amplitude
t_start = time.time()
FAIL = []
def check(name, ok):
    print(("ok   " if ok else "FAIL ") + name)
    if not ok:
        FAIL.append(name)
def close(a, b, digits=25):
    a, b = CC(SR(a).n(digits=digits + 15)), CC(SR(b).n(digits=digits + 15))
    return abs(a - b) <= 10^(-digits) * max(abs(b), 1e-300)

qed = QED(muons=False, taus=False)

# ---------------------------------------------------------------- tree level
# e+ e- -> mu+ mu- with the muon mass, (5.11) and (5.13) p. 136 (massless electron)
P = process("e- e+ -> mu- mu+", model=QED(taus=False), massless=["e"])
from feynsage.models import MM
E_, c_ = var('E_ c_')
assume(E_ > MM, MM > 0)
ours = P.squared().subs(P.to_angles()).subs({P.sqrt_s: 2*E_, P.cos_theta: c_})
check("e+e- -> mu+mu- (1/4) sum |M|^2 = Peskin (5.11), p. 136",
      bool((ours - EL^4*((1 + MM^2/E_^2) + (1 - MM^2/E_^2)*c_^2)).simplify_full() == 0))
al_ = var('al_')
sig = integrate(P.dsigma_dcos().subs({EL: sqrt(4*pi*al_)}).subs({P.sqrt_s: 2*E_, P.cos_theta: c_}), c_, -1, 1)
check("e+e- -> mu+mu- total cross section = Peskin (5.13), p. 136",
      bool((sig - 4*pi*al_^2/(3*(2*E_)^2)*sqrt(1 - MM^2/E_^2)*(1 + MM^2/(2*E_^2))).canonicalize_radical() == 0))
forget()

# Compton scattering, (5.87) p. 162: (1/4) sum |M|^2 in p.k and p.k'
P = process("e- gamma -> e- gamma", model=qed)
s, t, u = P.s, P.t, P.u
pk, pkp = (s - ME^2)/2, (ME^2 - u)/2
peskin = 2*EL^4*(pkp/pk + pk/pkp + 2*ME^2*(1/pk - 1/pkp) + ME^4*(1/pk - 1/pkp)^2)
check("Compton (1/4) sum |M|^2 = Peskin (5.87), p. 162",
      bool((P.squared() - peskin).subs(u=2*ME^2 - s - t).simplify_rational() == 0))

# e+ e- -> 2 gamma, (5.105) and (5.106) p. 168 (the overall minus sign of (5.105) is the crossing sign that
# the book says to remove; (5.106) has no 1/2 for the identical photons, dsigma_dcos() has it)
A = process("e- e+ -> gamma gamma", model=qed)
s, t, u = A.s, A.t, A.u
p1k1, p1k2 = (ME^2 - t)/2, (ME^2 - u)/2
peskin = 2*EL^4*(p1k2/p1k1 + p1k1/p1k2 + 2*ME^2*(1/p1k1 + 1/p1k2) - ME^4*(1/p1k1 + 1/p1k2)^2)
check("e+ e- -> 2 gamma (1/4) sum |M|^2 = Peskin (5.105) without the crossing sign, p. 168",
      bool((A.squared() - peskin).subs(u=2*ME^2 - s - t).simplify_rational() == 0))
E, c = var('E_', 'c_')
pp = sqrt(E^2 - ME^2)
al = EL^2/(4*pi)
peskin = 2*pi*al^2/(4*E^2)*(E/pp)*((E^2 + pp^2*c^2)/(ME^2 + pp^2*(1 - c^2)) + 2*ME^2/(ME^2 + pp^2*(1 - c^2))
                                    - 2*ME^4/(ME^2 + pp^2*(1 - c^2))^2)
assume(E > ME, ME > 0)
ours = 2*A.dsigma_dcos().subs({A.sqrt_s: 2*E, A.cos_theta: c})
check("e+ e- -> 2 gamma dsigma/dcos(theta) = Peskin (5.106), p. 168 (massive, exact)",
      bool((ours - peskin).canonicalize_radical().simplify_full() == 0))
forget()

# ---------------------------------------------------------------- one loop, QED
p = momenta("pp_")
P2 = var('P2_')
se = topologies(1, 1, 1, exclude=('tadpoles', 'wf'))
# photon self-energy: i Pi^{mu nu} = i M, Pi^{mu nu} = (q^2 g - q q) Pi_2;  g_{mu nu} Pi^{mu nu} = (d - 1) q^2 Pi_2
dA = insert_fields(se, ['A'], ['A'], qed)
gM = sum(dA.loop_amplitude([p], [p], lambda n, e: n*T.metric(e[0], e[1]), {"pp_^2": "P2_"}))
Pi2 = gM/((3 - 2*eps)*P2)
def Pi_hat(q2):                       # Pi_2(q^2) - Pi_2(0), Pi_2(0) from q^2 -> 0 on both sides
    f = lambda x: finite_part(Pi2.subs({P2: x}), 1)
    z = (f(10^-12) + f(-10^-12))/2
    return f(q2) - z
x = var('x_')
def mp_(r):                            # a Sage rational as an mpmath number (no Sage numbers inside mpmath)
    import mpmath
    r = QQ(r)
    return mpmath.mpf(int(r.numerator()))/int(r.denominator())
def peskin_791(q2):                   # (7.91) p. 252 at q^2 + i0 above threshold (q^2 > 4), m = 1, e = 1/3; plain
    import mpmath                     # Python ints (the r suffix): some mpmath builds reject Sage integers
    mpmath.mp.dps = 40r
    Q = mpmath.mpc(mp_(q2), mpmath.mpf('1e-30'))
    al_ = mpmath.mpf(1r)/9r/(4r*mpmath.pi)
    g = lambda xx: xx*(1r - xx)*mpmath.log(1r/(1r - xx*(1r - xx)*Q))
    r = mpmath.sqrt(mpmath.mpf(1r)/4r - 1r/mp_(q2))
    half = mpmath.mpf(1r)/2r
    v = -2r*al_/mpmath.pi*mpmath.quad(g, [0r, half - r, half, half + r, 1r])
    return CC(ComplexField(130)(str(v.real), str(v.imag)))
def integral(f, var_, at):
    """int_0^1 f dx by Maxima, Giac and SymPy: Maxima = Giac exactly; SymPy writes logs of negative numbers that
    Sage cannot cancel symbolically, so it is compared at a point to 40 digits"""
    r = {alg: integrate(f, x, 0, 1, algorithm=alg) for alg in ('maxima', 'giac', 'sympy')}
    assert bool((r['maxima'] - r['giac']).simplify_full() == 0), r
    assert abs(CC((r['maxima'] - r['sympy']).subs({var_: at}).n(digits=40))) < 1e-35, r
    return r['maxima']
# exact for spacelike q^2 = -Q^2: put Q^2 = (1 - y)^2/y (0 < y < 1), which makes every square root rational.
# Then 1 + x(1-x) Q^2 = (1 - x(1-y)) (y + x(1-y))/y, so Peskin's x integral has logs of positive numbers only;
# it is done by the computer (integral(): Maxima, Giac and SymPy must agree)
y_ = var('y_')
assume(y_ > 0, y_ < 1)
J = integral(x*(1 - x)*(log(1 - x*(1 - y_)) + log(y_ + x*(1 - y_)) - log(y_)), y_, 1/3)
peskin_sl = 2*(EL^2/(4*pi))/pi*J                          # -(2 alpha/pi) int x(1-x) log(1/(1 + x(1-x) Q^2))
Pi2y = disc_expand(finite_part(Pi2, 1)).subs({ME: 1}).subs({P2: -(1 - y_)^2/y_})
Pi2y = Pi2y.canonicalize_radical().log_expand().simplify_full()
d791 = Pi2y - Pi2y.limit(y_=1) - peskin_sl
check("photon self-energy Pi_2(q^2) - Pi_2(0) = Peskin (7.91), p. 252, exact for spacelike q^2 (m = 1)",
      bool(d791.simplify_full() == 0))
forget()
check("photon self-energy Pi_2(q^2) - Pi_2(0) = Peskin (7.91) above threshold, the +i0 branch (q^2 = 10, m = 1)",
      close(Pi_hat(10).subs({ME: 1, EL: 1/3}), peskin_791(10), 15))
check("photon self-energy pole: Pi_2 = -(e^2/(12 pi^2))/eps, Peskin (7.90), p. 252",
      bool((uv_part(Pi2) + EL^2/(12*pi^2)).simplify_full() == 0))

# electron self-energy, (10.41) p. 333 with the photon mass mu = 0 (dimensional regularization for both
# divergences); -i Sigma_2 = i M, Sigma_2 = A p/ + B m.  MS-bar normalization: drop -gamma + log(4 pi)
de = insert_fields(se, ['e'], ['e'], qed)
Am = sum(de.loop_amplitude([p], [p], lambda n, e: dirac_trace(slash(p)*n, dim='d')/(4*P2), {"pp_^2": "P2_"}))
Bm = sum(de.loop_amplitude([p], [p], lambda n, e: dirac_trace(n, dim='d')/4, {"pp_^2": "P2_"}))
# exact, for 0 < p^2 < m^2 = 1: Peskin's x integrals done by the computer.  Delta = (1 - x)(1 - x p^2), so
# log Delta = log(1 - x) + log(1 - x p^2), both logs of positive numbers
assume(P2 > 0, P2 < 1)
logD = log(1 - x) + log(1 - x*P2)
k_ = EL^2/(16*pi^2)
fa = k_*integral(2*x*logD + 2*x, P2, 3/10)                        # p/ part of (10.41)
fb = k_*integral(-4*logD - 2, P2, 3/10)                           # m part
ex = lambda F: disc_expand(SR(F).subs({ME: 1})).subs({LogM(1 - P2): log(1 - P2)})
check("electron self-energy Sigma_2(p) = Peskin (10.41), p. 333, exact for 0 < p^2 < m^2 (m = 1)",
      bool((ex(-finite_part(Am, 1)) - fa).simplify_full() == 0 and (ex(-finite_part(Bm, 1)) - fb).simplify_full() == 0))
forget()

# ---------------------------------------------------------------- one loop, QCD: the beta function
sm = SM()
ALL = ['A', 'Z', 'W+', 'H', 'G0', 'G+', 'u+', 'u-', 'uA', 'uZ', 'g', 'ug', 'e', 'mu', 'ta', 'ne', 'nm', 'nt',
       'u', 'c', 't', 'd', 's', 'b']
only = lambda keep: [f for f in ALL if f not in keep]
a1, a2 = color_indices("a1 a2"); i1, i2 = quark_colors("i1 i2")
kin = {"pp_^2": "P2_"}
# delta_3 = Pi(k^2) pole (gluon, ghost and one quark loop), Peskin (16.74) p. 527
dg = insert_fields(se, ['g'], ['g'], sm, exclude_fields=only(('g', 'ug', 'u')))
gM = sum(dg.loop_amplitude([p], [p], lambda n, e: n*T.metric(e[0], e[1]), kin)).subs({color_delta(a1, a2): 1})
c3 = (uv_part(gM)/(3*P2)).simplify_full()
check("delta_3 = g^2/(16 pi^2 eps) (5/3 C2(G) - 4/3 nf C(r)) with nf = 1 = Peskin (16.74), p. 527",
      bool((c3 - GS^2/(16*pi^2)*(5 - 2/3)).simplify_full() == 0))
# delta_2 = the p/ coefficient of Sigma_2 (Sigma = -M), quark self-energy
dq = insert_fields(se, ['u'], ['u'], sm, exclude_fields=only(('g', 'u')))
Aq = sum(dq.loop_amplitude([p], [p], lambda n, e: dirac_trace(slash(p)*n, dim='d')/(4*P2), kin))
c2 = (-uv_part(Aq)).subs({color_delta(i1, i2): 1}).simplify_full()
# delta_1 = minus the pole of the quark-gluon vertex in units of the tree vertex: project with Tr[gamma_mu X^mu]/(4 d)
q1, q2 = momenta("q1_ q2_")
dv = insert_fields(topologies(1, 2, 1, exclude=('tadpoles', 'wf')), ['u', 'g'], ['u'], sm, exclude_fields=only(('g', 'ug', 'u')))
a = color_indices("a2")[0] if isinstance(color_indices("a2"), (list, tuple)) else color_indices("a2")
proj = lambda n, e: dirac_trace(gamma(e[1])*n, dim='d')/16
kv = {"q1_^2": 0, "q2_^2": 0, "q1_.q2_": "-P2_/2"}
Mv = sum(dv.loop_amplitude([q1, q2 - q1], [q2], proj, kv))
# the tree vertex u(i1) g(a2) -> u(i3): M = -g_s T^a gamma^mu (i M = -i g_s T^a gamma^mu); contract the colour with T^a_{i1 i3}
i3 = quark_colors("i3")[0] if isinstance(quark_colors("i3"), (list, tuple)) else quark_colors("i3")
Tc = T_color(a, i1, i3)
num = color_factor(uv_part(Mv)*Tc, N=3)
den = color_factor(-GS*T_color(a, i3, i1)*Tc, N=3)
c1 = -(num/den).simplify_full()
beta = (-2*GS*(-c1 + c2 + c3/2)).simplify_full()
check("delta_1, delta_2, delta_3 from feynsage's diagrams give beta = -g^3/(16 pi^2) (11 - 2 nf/3), nf = 1: "
      "Peskin (16.73) p. 527 and (16.85) p. 531",
      bool((beta + GS^3/(16*pi^2)*(11 - 2/3)).simplify_full() == 0))

# ---------------------------------------------------------------- the same in a general covariant gauge
# J. C. Collins, Foundations of Perturbative QCD (Cambridge, 2011), Sec. 3.3, read on the printed pages 47-48:
# gluon propagator -g + (1 - xi) k k/k^2, the same convention as feynsage's R_xi; here xi is a symbol, given to
# the gluon only (SM(xi={"g": xi}) / loop_amplitude(..., xi={"g": xi})).  C_F = 4/3, C_A = 3, T_F = 1/2, n_f = 1.
xg = var('xi_g')
XG = {"g": xg}
CF_, CA_, k_ = QQ(4)/3, 3, GS^2/(16*pi^2)
Aq_x = sum(dq.loop_amplitude([p], [p], lambda n, e: dirac_trace(slash(p)*n, dim='d')/(4*P2), kin, xi=XG))
Bq_x = sum(dq.loop_amplitude([p], [p], lambda n, e: dirac_trace(n, dim='d')/4, kin, xi=XG))
pole = lambda e: SR(uv_part(e)).subs({color_delta(i1, i2): 1}).simplify_full()
# (3.22): the pole of the graph i M is i g^2 C_F/(16 pi^2 eps) [-3 m + xi (p/ - m)]
check("quark self-energy pole in a general gauge, i g^2 C_F/(16 pi^2 eps) [-3m + xi (p/ - m)]: Collins (3.22), p. 47",
      bool((pole(Aq_x) - k_*CF_*xg).simplify_full() == 0) and bool((pole(Bq_x) - k_*CF_*(-3 - xg)*MU).simplify_full() == 0))
gM_x = sum(dg.loop_amplitude([p], [p], lambda n, e: n*T.metric(e[0], e[1]), kin, xi=XG)).subs({color_delta(a1, a2): 1})
c3_x = (uv_part(gM_x)/(3*P2)).simplify_full()
check("gluon self-energy: Z3 - 1 = (alpha_s/(4 pi eps)) [(13/6 - xi/2) C_A - 4/3 T_F n_f]: Collins (3.25), p. 47",
      bool((c3_x - k_*((QQ(13)/6 - xg/2)*CA_ - QQ(2)/3)).simplify_full() == 0))
Mv_x = sum(dv.loop_amplitude([q1, q2 - q1], [q2], proj, kv, xi=XG))
c1_x = -(color_factor(uv_part(Mv_x)*Tc, N=3)/den).simplify_full()
check("quark-gluon vertex: Z1 - 1 = -(alpha_s/(4 pi eps)) [xi (C_F - C_A/2) + C_A (3/4 + 3 xi/4)]: Collins (3.29)-(3.30), p. 48",
      bool((c1_x + k_*(xg*(CF_ - CA_/2) + CA_*(QQ(3)/4 + 3*xg/4))).simplify_full() == 0))
# the ghost: as for the quark (Z2 - 1 = minus the p/ pole of M), Z~ - 1 = minus the p^2 pole of M
dgh = insert_fields(se, ['ug'], ['ug'], sm, exclude_fields=only(('g', 'ug')))
Mgh = sum(dgh.loop_amplitude([p], [p], lambda n, e: n, kin, xi=XG)).subs({color_delta(a1, a2): 1})
check("ghost self-energy: Z~ - 1 = (alpha_s/(4 pi eps)) C_A (3/4 - xi/4): Collins (3.26), p. 47",
      bool((-uv_part(Mgh)/P2 - k_*CA_*(QQ(3)/4 - xg/4)).simplify_full() == 0))
c2_x = -pole(Aq_x)
beta_x = (-2*GS*(-c1_x + c2_x + c3_x/2)).simplify_full()
check("and the beta function from them does not depend on xi: g0 of Collins (3.31), p. 48 (11/6 C_A - 2/3 T_F n_f)",
      bool((beta_x + GS^3/(16*pi^2)*(11 - QQ(2)/3)).simplify_full() == 0))

# ---------------------------------------------------------------- QED in a general covariant gauge: Z1 = Z2
# S. Pokorski, Gauge Field Theories, 2nd ed. (Cambridge, 2000), p. 193, (5.71): (Z1 - 1) = -(alpha/(2 pi)) 1/eps_UV in
# minimal subtraction, his n = 4 - eps (so his 1/eps is our 1/(2 eps)), and Z1 = Z2 (the Ward identity).  With a symbolic
# photon gauge parameter the vertex and self-energy poles both depend on xi; the Ward identity says they stay equal.
xa = var('xi_A')
de_x = insert_fields(se, ['e'], ['e'], qed)
A_x = sum(de_x.loop_amplitude([p], [p], lambda n, e: dirac_trace(slash(p)*n, dim='d')/(4*P2), {"pp_^2": "P2_"}, xi={"A": xa}))
dvq = insert_fields(topologies(1, 2, 1, exclude=('tadpoles', 'wf')), ['e', 'A'], ['e'], qed)
Mvq = sum(dvq.loop_amplitude([q1, q2 - q1], [q2], proj, kv, xi={"A": xa}))
d2_qed = (-uv_part(A_x)).simplify_full()
d1_qed = (-uv_part(Mvq)/EL).simplify_full()           # the tree vertex is M = e gamma^mu (Feynman rule i e gamma^mu)
check("QED Ward identity Z1 = Z2 in a general gauge (symbolic xi), = -e^2/(16 pi^2 eps) at xi = 1: Pokorski (5.71), p. 193",
      bool((d1_qed - d2_qed).simplify_full() == 0) and bool((d2_qed.subs({xa: 1}) + EL^2/(16*pi^2)).simplify_full() == 0))

# ---------------------------------------------------------------- electroweak processes
# A. Lahiri and P. B. Pal, A First Book of Quantum Field Theory, 2nd edition (preface July 2004), Ch. 15, read on
# the printed pages 343-351.  Couplings: g = e/sin(theta_W), alpha = e^2/(4 pi); Z couplings (15.56)-(15.57) and
# (15.72)-(15.73); Feynman-'t Hooft gauge; electron and muon masses zero where the book drops them.
from feynsage.models import SW, CW, MZ, MW, MH, MM
al_ew = EL^2/(4*pi)
simp_ew = lambda e: SR(e).subs({MW: MZ*CW}).subs({CW: sqrt(1 - SW^2)}).simplify_full()
for name, anti, t3, Q in (("nu_e", "nu_e~", QQ(1)/2, 0), ("e-", "e+", -QQ(1)/2, -1)):
    af, bf = t3 - 2*Q*SW^2, t3
    G = process("Z -> %s %s" % (name, anti), massless=["e"]).width()
    check("Gamma(Z -> %s %s) = alpha M_Z (a_f^2 + b_f^2)/(12 sin^2 cos^2): Lahiri & Pal (15.69), p. 345" % (name, anti),
          bool(simp_ew(G - al_ew*MZ/(12*SW^2*CW^2)*(af^2 + bf^2)) == 0))
assume(MH > 2*MM, MM > 0)
G = process("H -> mu- mu+").width()
check("Gamma(H -> mu- mu+) = alpha m^2 M_H/(8 sin^2 cos^2 M_Z^2) (1 - 4 m^2/M_H^2)^(3/2): Lahiri & Pal (15.70), p. 346",
      bool(simp_ew(G - al_ew*MM^2*MH/(8*SW^2*CW^2*MZ^2)*(1 - 4*MM^2/MH^2)^(3/2)).canonicalize_radical() == 0))
forget()
P = process("e- e+ -> mu- mu+", massless=["e", "mu"])
s_ew, c_ew = P.sqrt_s^2, P.cos_theta
gV, gA = -QQ(1)/2 + 2*SW^2, -QQ(1)/2
fZ = s_ew/(4*SW^2*CW^2*(s_ew - MZ^2))
a1_ew = 2*fZ*gV^2 + fZ^2*(gV^2 + gA^2)^2
a2_ew = 4*fZ*gA^2 + 8*fZ^2*gV^2*gA^2
check("e+e- -> mu+mu- with gamma and Z: |M|^2 = e^4 [(1 + a1)(1 + cos^2) + a2 cos]: Lahiri & Pal (15.81)-(15.82), p. 348",
      bool(simp_ew(SR(P.squared()).subs(P.to_angles()) - EL^4*((1 + a1_ew)*(1 + c_ew^2) + a2_ew*c_ew)) == 0))
ds = P.dsigma_dcos()
AFB = (integrate(ds, c_ew, 0, 1) - integrate(ds, c_ew, -1, 0))/integrate(ds, c_ew, -1, 1)
check("forward-backward asymmetry A = (3/8) a2/(1 + a1): Lahiri & Pal (15.84), p. 348",
      bool(simp_ew(AFB - QQ(3)/8*a2_ew/(1 + a1_ew)) == 0))
for text, gv, ga, ref_ in (("nu_mu e- -> nu_mu e-", -QQ(1)/2 + 2*SW^2, -QQ(1)/2, "(15.89), p. 350"),
                           ("nu_e e- -> nu_e e-", QQ(1)/2 + 2*SW^2, QQ(1)/2, "(15.89) with (15.93), pp. 350-351")):
    N = process(text)
    s_, t_, u_ = N.s, N.t, N.u
    GF = EL^2/(4*sqrt(2)*SW^2*MW^2)                              # (15.88)
    ref = 16*GF^2*((gv + ga)^2*((s_ - ME^2)/2)^2 + (gv - ga)^2*((ME^2 - u_)/2)^2 - (gv^2 - ga^2)*ME^2*(-t_/2))
    lim = lambda e: (SR(e).subs({u_: 2*ME^2 - s_ - t_}).subs({MW: MZ*CW})*MZ^4).limit(MZ=oo)
    check("%s at energies << M_Z: |M|^2 = 16 G_F^2 [...]: Lahiri & Pal %s" % (text, ref_),
          bool(simp_ew(lim(N.squared()) - lim(ref)) == 0))
v_ew = {EL: sqrt(4*pi/128), SW: sqrt(QQ(23)/100), CW: sqrt(QQ(77)/100), MZ: 91, MW: 80}
check("Gamma(Z -> nu nu~) 'about 165 MeV' and Gamma(W+ -> e+ nu_e) 'about 220 MeV': Lahiri & Pal Ex. 15.10-15.11, p. 345 "
      "(alpha = 1/128; the book does not give alpha)",
      abs(process("Z -> nu_e nu_e~").width().subs(v_ew).n()/0.165 - 1) < 0.1 and
      abs(process("W+ -> e+ nu_e", massless=["e"]).width().subs(v_ew).n()/0.220 - 1) < 0.1)

# ---------------------------------------------------------------- bound states, Peskin Sec. 5.3 and Problem 5.4
# the pair at rest: sum v ubar = Pi_0 = gamma5 (P/ + M)/(2 sqrt 2), Pi_1 = eps*/ (P/ + M)/(2 sqrt 2); for annihilation
# u vbar = gamma0 Pi^dagger gamma0.  M(-> B) = sqrt(2/M) psi(0) times the free amplitude with the spins traced against Pi
Pb, p1b, p2b, epc, k1b, epb = momenta("Pb p1b p2b epc k1b epb")
mub, nub, rhob, sib = lorentz_indices("mub nub rhob sib")
Mb, psi2 = var('Mb psi2')
rb = {dot(p1b, p1b): 0, dot(p2b, p2b): 0, dot(p1b, p2b): Mb^2/2, dot(Pb, Pb): Mb^2, dot(p1b, Pb): Mb^2/2, dot(p2b, Pb): Mb^2/2}
Tb = dirac_trace(gamma(mub)*slash(epc)*(slash(Pb) + Mb)/(2*sqrt(2)), rules=rb)
cb = Tb.coefficient(comp(epc, mub))
Lb = dirac_trace(slash(p2b)*gamma(mub)*slash(p1b)*gamma(nub), rules=rb)
Msq_b = (EL^2/Mb^2)^2*(2/Mb)*psi2*cb^2*contract(Lb*(-metric(mub, nub) + comp(Pb, mub)*comp(Pb, nub)/Mb^2), rules=rb)
al_ = EL^2/(4*pi)
check("sigma(e+e- -> B) = 64 pi^3 alpha^2 |psi(0)|^2/M^3 delta(E^2 - M^2): Peskin (5.54), p. 151",
      bool((2*pi/(2*Mb^2)*Msq_b/4 - 64*pi^3*al_^2*psi2/Mb^3).simplify_full() == 0))
check("Gamma(B -> e+e-) = 16 pi alpha^2 |psi(0)|^2/(3 M^2): Peskin (5.57), p. 151",
      bool((Msq_b/3/(2*Mb)/(8*pi) - 16*pi*al_^2*psi2/(3*Mb^2)).simplify_full() == 0))
rp = {dot(Pb, Pb): Mb^2, dot(k1b, k1b): 0, dot(Pb, k1b): Mb^2/2, dot(epb, Pb): 0}
G2b = gamma(nub)*(slash(Pb/2 - k1b) + Mb/2)*gamma(mub)/(-Mb^2/2) + gamma(mub)*(slash(k1b - Pb/2) + Mb/2)*gamma(nub)/(-Mb^2/2)
A1b = dirac_trace(G2b*(slash(Pb) + Mb)*slash(epb)/(2*sqrt(2)), rules=rp)
A0b = dirac_trace(G2b*(-(slash(Pb) + Mb)*gamma5()/(2*sqrt(2))), rules=rp)
S0b = contract(A0b*conjugate(A0b).subs({mub: rhob, nub: sib})*metric(mub, rhob)*metric(nub, sib), rules=rp)
me_ = var('me_')
Gps = (1/2/(2*Mb)/(8*pi)*EL^4*(2/Mb)*psi2*S0b).subs({psi2: (me_*al_)^3/(8*pi), Mb: 2*me_})
check("positronium: ortho (S = 1) does not decay to 2 gamma, para gives Gamma = alpha^5 m/2: Peskin Problem 5.4 (a), p. 171",
      A1b == 0 and bool((Gps - al_^5*me_/2).simplify_full() == 0))

print("time %.0f s" % (time.time() - t_start))
print("ALL BOOK CHECKS PASS" if not FAIL else "FAILED: %s" % FAIL)
