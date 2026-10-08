# Tests of process(): the one-line front end (kinematics, averaging, polarization and colour sums,
# helicity states, widths, cross sections, parallel pairs).
#     sage tests/test_process.sage
import time, itertools
from feynsage import *
from feynsage.models import EL, SW, CW, MW, MZ, MH, MB, GS, QED
t_start = time.time()
FAIL = []
def check(name, ok):
    print(("ok   " if ok else "FAIL ") + name)
    if not ok:
        FAIL.append(name)

P = process("e- e+ -> mu- mu+", model=QED(), massless=["e", "mu"])
s, t, u = P.s, P.t, P.u
check("QED e+e- -> mu+mu-: |M|^2 = 2 e^4 (t^2 + u^2)/s^2", (P.squared() - 2*EL^4*(t^2 + u^2)/s^2).subs(u=-s-t).simplify_rational() == 0)
check("QED dsigma/dcos = pi alpha^2 (1 + cos^2)/(2 s) [Peskin 5.12]",
      (P.dsigma_dcos() - pi*(EL^2/(4*pi))^2*(1 + P.cos_theta^2)/(2*P.sqrt_s^2)).simplify_full() == 0)
check("e-_R e+_L -> mu-_R mu+_L = e^4 (1 + cos)^2 [Peskin 5.21, p. 143]",
      (P.squared(helicities={"e-": +1, "e+": -1, "mu-": +1, "mu+": -1}) - EL^4*(1 + P.cos_theta)^2).simplify_full() == 0)
check("e-_R e+_R -> mu mu = 0", P.squared(helicities={"e-": +1, "e+": +1}).simplify_full() == 0)
al = 1/137.035999084
check("sigma(10 GeV) = 4 pi alpha^2/(3 s) to 1e-12", abs(P.sigma(10, unit="pb", alpha=al)/float(4*pi*al^2/300*0.3893793721e9) - 1) < 1e-12)

W = process("e- e+ -> W+ W-", massless=["e"])
num = {EL: 1, MW: 804/10, MZ: (804/10)/sqrt(77/100), SW: sqrt(23/100), CW: sqrt(77/100)}
rs, c = 500, -QQ(1)/2
mw, mz, sw2 = 804/10, (804/10)/sqrt(77/100), 23/100
b = sqrt(1 - 4*mw^2/rs^2)
x29 = rs^2/(rs^2 - mz^2)*mz^2/(4*mw^2)*b*(b^2 - 3)*sqrt(1 - c^2)
v = W.squared(helicities={"e-": +1, "e+": -1, "W+": 0, "W-": 0}).subs(num).subs({W.sqrt_s: rs, W.cos_theta: c})
check("e-_R e+_L -> W+_0 W-_0 = Xianyu (21.29)", abs(v.n(digits=40)/x29^2 - 1) < 1e-25)
D = 1 + b^2 + 2*b*c
x31 = mz^2/(rs^2 - mz^2)*b*sqrt(1 - c^2)
v = W.squared(helicities={"e-": +1, "e+": -1, "W+": 1, "W-": 1}).subs(num).subs({W.sqrt_s: rs, W.cos_theta: c})
check("e-_R e+_L -> W+_+ W-_+ = Xianyu (21.31)", abs(v.n(digits=40)/x31^2 - 1) < 1e-25)
cP = -c
br = (1/(2*sw2)*(-rs^2/(rs^2 - mz^2)*(mz^2/(2*mw^2) + 1) + 2/b^2 - 8*mw^2/(rs^2*b^2*(1 + b^2 - 2*b*cP))) + mz^2/mw^2*((rs^2/2 + mw^2)/(rs^2 - mz^2)))/rs^2
E = rs/2; p = sqrt(E^2 - mw^2)
v = W.squared(helicities={"e-": -1, "e+": +1, "W+": 0, "W-": 0}).subs(num).subs({W.sqrt_s: rs, W.cos_theta: c})
check("e-_L e+_R -> W+_0 W-_0 = Peskin (21.108)", abs(v.n(digits=40)/(4*E*p*sqrt(1 - c^2)*br)^2 - 1) < 1e-25)
a1 = W.squared(nproc=1); a2 = W.squared(nproc=4)
pt = {W.s: 1000^2, W.t: -300000, W.u: 2*80^2 - 1000^2 + 300000}
check("e+e- -> W+W-: parallel pairs = one product", ((a1 - a2).subs(pt).subs(num)).n(digits=30) == 0 or abs((a1 - a2).subs(pt).subs(num).n(digits=30)) < 1e-25)

Q = process("u u~ -> g g", exclude_fields=['A', 'Z', 'H', 'G0', 'G+', 'W+'], massless=["u"])
s, t, u = Q.s, Q.t, Q.u
check("u u~ -> g g = (32/27) g^4 [t/u + u/t - 9(t^2+u^2)/(4s^2)] [Xianyu (17.23), Peskin problem 17.3]",
      (Q.squared() - 32*GS^4/27*(t/u + u/t - 9*(t^2 + u^2)/(4*s^2))).subs(s=-t-u).simplify_rational() == 0)
cfg = [{1: a, 2: b_, 3: c_, 4: d} for a, b_, c_, d in itertools.product((1, -1), repeat=4)]
tab = Q.helicity_table(cfg)
pt = {Q.sqrt_s: 7, Q.cos_theta: QQ(1)/3, GS: 1}
check("u u~ -> g g: sum over 16 helicity states = unpolarized |M|^2",
      abs((sum(tab).subs(pt)/Q.squared(average=False).subs(Q.to_angles()).subs(pt)).n(digits=30) - 1) < 1e-25)

H = process("H -> b b~")
ref = 3*EL^2*MH/(32*pi*SW^2)*MB^2/MW^2*(1 - 4*MB^2/MH^2)^(3/2)
check("Gamma(H -> b b~) = 3 alpha m_H m_b^2 beta^3/(8 sin^2 m_W^2)", ((H.width()/ref).subs({MH: 125, MB: 5, MW: 80, SW: sqrt(23/100), EL: 1}) - 1).n(digits=30) == 0)

from feynsage._parse import swap_i
x_ = var('x_')
check("conjugation: I*x -> -I*x (Sage's subs({I: -I}) leaves a lone I*monomial alone)",
      bool(swap_i(3*I*x_, -I) == -3*I*x_) and bool(swap_i((2 + I)*x_/(x_ + I), -I) == (2 - I)*x_/(x_ - I)))

G = process("g g -> u u~", exclude_fields=['A', 'Z', 'H', 'G0', 'G+', 'W+'])
s, t, u = G.s, G.t, G.u
from feynsage.models import MU
t1, t2, rho = (MU^2 - t)/s, (MU^2 - u)/s, 4*MU^2/s
ref = GS^4*(1/(6*t1*t2) - 3/8)*(t1^2 + t2^2 + rho - rho^2/(4*t1*t2))
check("g g -> Q Q~ with the quark mass = FeynCalc example QCD/Tree/GlGl-QQbar (massive result)",
      all(abs((G.squared()/ref).subs({s: a, t: b, u: 2*m^2 - a - b, MU: m}).n(digits=30) - 1) < 1e-25
          for a, b, m in ((100, -30, 3), (50, -7, 2), (1000, -300, 10))))

N = process("e- nu_mu -> nu_e mu-", massless=["e"])
s, t, u = N.s, N.t, N.u
from feynsage.models import MM
GF = EL^2/(4*sqrt(2)*MW^2*SW^2)
lim = 16*GF^2*s*(s - MM^2)            # M_W -> infinity, only the e- spin averaged [FeynCalc EW/Tree/ElNmu-MuNel]
pt = {s: 3, t: -1, u: MM^2 - 3 + 1, MM: 1/10, MW: 10^6, SW: 1/2, EL: 1/3}
check("e- nu_mu -> nu_e mu-: averaged over the e- spin only (one neutrino helicity)",
      abs((N.squared()/lim).subs(pt).n(digits=30) - 1) < 1e-9)

print("time %.0f s" % (time.time() - t_start))
print("ALL PROCESS CHECKS PASS" if not FAIL else "FAILED: %s" % FAIL)
