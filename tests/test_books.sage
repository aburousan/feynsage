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
ours = 2*A.dsigma_dcos()
pt = {A.sqrt_s: 2*E, A.cos_theta: c}
check("e+ e- -> 2 gamma dsigma/dcos(theta) = Peskin (5.106), p. 168 (massive, two points)",
      all(close(SR(ours).subs(pt).subs({E: e_, c: c0, ME: 1, EL: 1/3}), peskin.subs({E: e_, c: c0, ME: 1, EL: 1/3}))
          for e_, c0 in ((3, 1/3), (7/4, -5/7))))

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
def peskin_791(q2):                   # (7.91) p. 252, q^2 + i0 above threshold, at m = 1 and e = 1/3
    import mpmath
    mpmath.mp.dps = 40
    Q = mpmath.mpc(mp_(q2), mpmath.mpf(10)**-30)
    al_ = mpmath.mpf(1)/9/(4*mpmath.pi)
    g = lambda xx: xx*(1 - xx)*mpmath.log(1/(1 - xx*(1 - xx)*Q))
    pts = [0, 0.5, 1]
    if mp_(q2) > 4:
        r = mpmath.sqrt(mpmath.mpf(1)/4 - 1/mp_(q2))
        pts = [0, mpmath.mpf(1)/2 - r, mpmath.mpf(1)/2, mpmath.mpf(1)/2 + r, 1]
    v = -2*al_/mpmath.pi*mpmath.quad(g, pts)
    return CC(ComplexField(130)(str(v.real), str(v.imag)))
check("photon self-energy Pi_2(q^2) - Pi_2(0) = Peskin (7.91), p. 252, below and above threshold (m = 1)",
      all(close(Pi_hat(q2).subs({ME: 1, EL: 1/3}), peskin_791(q2), 15) for q2 in (-3, 1/2, 10)))
check("photon self-energy pole: Pi_2 = -(e^2/(12 pi^2))/eps, Peskin (7.90), p. 252",
      bool((uv_part(Pi2) + EL^2/(12*pi^2)).simplify_full() == 0))

# electron self-energy, (10.41) p. 333 with the photon mass mu = 0 (dimensional regularization for both
# divergences); -i Sigma_2 = i M, Sigma_2 = A p/ + B m.  MS-bar normalization: drop -gamma + log(4 pi)
de = insert_fields(se, ['e'], ['e'], qed)
Am = sum(de.loop_amplitude([p], [p], lambda n, e: dirac_trace(slash(p)*n, dim='d')/(4*P2), {"pp_^2": "P2_"}))
Bm = sum(de.loop_amplitude([p], [p], lambda n, e: dirac_trace(n, dim='d')/4, {"pp_^2": "P2_"}))
def peskin_1041(p2):
    # Sigma_2 = e^2/(16 pi^2) int dx Gamma(eps) e^(gamma eps) Delta^(-eps) ((4 - 2 eps) m - (2 - 2 eps) x p/),
    # Delta = (1 - x) m^2 - x (1 - x) p^2  (eps here = (4 - d)/2, Peskin's eps/2); finite parts at m = 1, e = 1/3:
    # Gamma(eps) e^(gamma eps) Delta^(-eps) (a - b eps) = a/eps - a log Delta - b + O(eps)
    import mpmath
    mpmath.mp.dps = 40
    P = mp_(p2)
    D = lambda xx: (1 - xx) - xx*(1 - xx)*P
    fa = mpmath.quad(lambda xx: 2*xx*mpmath.log(D(xx)) + 2*xx, [0, 1])       # p/: a = -2x, b = -2x
    fb = mpmath.quad(lambda xx: -4*mpmath.log(D(xx)) - 2, [0, 1])           # m:  a = 4,   b = 2
    k = mpmath.mpf(1)/9/(16*mpmath.pi**2)
    return RealField(130)(str(k*fa)), RealField(130)(str(k*fb))
check("electron self-energy Sigma_2(p) = Peskin (10.41), p. 333, with mu = 0 (p^2 = 0.3 and -2, m = 1)",
      all(close(-finite_part(Am, 1).subs({P2: p2, ME: 1, EL: 1/3}), peskin_1041(p2)[0], 15) and
          close(-finite_part(Bm, 1).subs({P2: p2, ME: 1, EL: 1/3}), peskin_1041(p2)[1], 15) for p2 in (3/10, -2)))

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

print("time %.0f s" % (time.time() - t_start))
print("ALL BOOK CHECKS PASS" if not FAIL else "FAILED: %s" % FAIL)
