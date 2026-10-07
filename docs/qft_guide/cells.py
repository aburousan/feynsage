# The code cells of "Dirac algebra with feynsage, a beginner's guide", in order.  Each is
# (id, code, options); the cells share one namespace, like a notebook.  fig=<name> saves that
# matplotlib figure, text=True shows the value as plain text.
CELLS = [

# ------------------------------------------------------------------ first steps
("import", """
from feynsage import *
""", {}),

("symbols", """
p, q, k = momenta("p q k")
mu, nu, rho, tau = lorentz_indices("mu nu rho tau")
m = var('m')
""", {}),

("first", """
slash(p) * gamma(mu) * (slash(q) + m)
""", {}),

("order", """
gamma(mu) * gamma(nu) - gamma(nu) * gamma(mu)
""", {}),

("eulergamma", """
gamma(5), gamma(1/2)
""", {}),

# ------------------------------------------------------------------ traces
("tr2", """
dirac_trace(gamma(mu) * gamma(nu))
""", {}),

("tr4", """
dirac_trace(gamma(mu) * gamma(nu) * gamma(rho) * gamma(tau))
""", {}),

("tr3", """
dirac_trace(gamma(mu) * gamma(nu) * gamma(rho))
""", {}),

("trpq", """
dirac_trace(slash(p) * slash(q))
""", {}),

("trmix", """
dirac_trace(slash(p) * gamma(mu) * slash(q) * gamma(nu))
""", {}),

("trmass", """
dirac_trace((slash(p) + m) * gamma(mu) * (slash(q) + m) * gamma(nu))
""", {}),

("method", """
expr = slash(p) * slash(q) * slash(p) * slash(q)
expr.trace()
""", {}),

("batch", """
dirac_trace([gamma(mu) * gamma(mu), slash(p) * slash(p), 7])
""", {}),

# ------------------------------------------------------------------ scalar products and kinematics
("dots", """
dot(p + k, q)
""", {}),

("dotsym", """
dot(q, p) - dot(p, q)
""", {}),

("subs", """
s = var('s')
r = dirac_trace(slash(p) * slash(q) * slash(p) * slash(q))
r.subs({dot(p, p): m^2, dot(q, q): 0, dot(p, q): s/2})
""", {}),

("rules", """
dirac_trace(slash(p) * slash(q) * slash(p) * slash(q),
            rules={dot(p, p): m^2, dot(q, q): 0, dot(p, q): s/2})
""", {}),

("momrule", """
p1, p2, k1, k2 = momenta("p1 p2 k1 k2")
dirac_trace(slash(k2) * slash(p1), rules={k2: p1 + p2 - k1})
""", {}),

# ------------------------------------------------------------------ contractions
("contract1", """
contract(metric(mu, nu) * comp(p, mu) * comp(q, nu))
""", {}),

("contract2", """
T = dirac_trace(slash(p) * gamma(mu) * slash(q) * gamma(nu))
contract(T * comp(k, mu) * comp(k, nu))
""", {}),

("epseps", """
contract(epsilon(mu, nu, rho, tau) * epsilon(mu, nu, rho, tau))
""", {}),

("epseps2", """
contract(epsilon(tau, nu, rho, mu) * epsilon(tau, nu, rho, k))
""", {}),

("open1", """
contract(gamma(mu) * gamma(nu) * gamma(mu))
""", {}),

("open2", """
contract(gamma(mu) * gamma(nu) * gamma(rho) * gamma(mu))
""", {}),

("open3", """
contract(gamma(mu) * slash(p) * slash(q) * gamma(mu))
""", {}),

("ddim", """
d = var('d')
contract(gamma(mu) * gamma(nu) * gamma(mu), dim=d)
""", {}),

("ddimtr", """
dirac_trace(gamma(mu) * gamma(nu) * gamma(mu) * gamma(nu), dim=d)
""", {}),

# ------------------------------------------------------------------ gamma5
("g5tr", """
dirac_trace(gamma(mu) * gamma(nu) * gamma(rho) * gamma(tau) * gamma5())
""", {}),

("g5slash", """
dirac_trace(slash(p) * gamma(mu) * gamma5() * slash(q) * gamma(nu))
""", {}),

("projectors", """
PL()
""", {}),

("plpr", """
contract(PL() * PR()), contract(PL() * PL())
""", {}),

("ndr", """
dirac_trace(PL() * gamma(mu) * PR() * gamma(nu), dim=d)
""", {}),

("ndrerror", """
try:
    dirac_trace(gamma(mu) * gamma(nu) * gamma(rho) * gamma(tau) * gamma5(), dim=d)
except DiracError as err:
    print(err)
""", {}),

# ------------------------------------------------------------------ spinors
("spinor", """
J = ubar(p2, m) * gamma(mu) * u(p1, m)
J
""", {}),

("conj", """
conjugate(J)
""", {}),

("diraceq", """
contract(ubar(p2, m) * slash(p2) * gamma(mu) * slash(p1) * u(p1, m))
""", {}),

("diraceqv", """
contract(vbar(p2, m) * gamma(mu) * gamma5() * slash(p1) * v(p1, m))
""", {}),

("wrongorder", """
try:
    gamma(mu) * ubar(p, m)
except DiracError as err:
    print(err)
""", {}),

# ------------------------------------------------------------------ e+ e- -> mu+ mu-
("ee_M", """
M_mu = var('M_mu')
M = vbar(p2) * gamma(mu) * u(p1) * ubar(k1, M_mu) * gamma(mu) * v(k2, M_mu)
M
""", {}),

("ee_Mc", """
conjugate(M)
""", {}),

("ee_sum", """
Msq = spin_sum(M * conjugate(M)) / 4
Msq
""", {}),

("ee_frame", """
E, K, c, e = var('E K c e')        # beam energy, |k|, cos(theta), charge
assume(E > M_mu, M_mu > 0)
frame = {dot(p1, p2): 2*E^2,
         dot(k1, p1): E^2 - E*K*c, dot(k2, p2): E^2 - E*K*c,
         dot(k2, p1): E^2 + E*K*c, dot(k1, p2): E^2 + E*K*c}
msq = (e^4 / (4*E^2)^2 * Msq).subs(frame).subs(K=sqrt(E^2 - M_mu^2))
msq.full_simplify()
""", {}),

("ee_check511", """
peskin_511 = e^4 * ((1 + M_mu^2/E^2) + (1 - M_mu^2/E^2) * c^2)
(msq - peskin_511).full_simplify()
""", {}),

("ee_dsigma", """
a_em = var('a_em', latex_name=r'\\alpha')
Ecm = 2*E
dsigma = 1/(2*Ecm^2) * sqrt(E^2 - M_mu^2)/(16*pi^2*Ecm) * msq.subs(e=sqrt(4*pi*a_em))
sigma_tot = (2*pi * integrate(dsigma, c, -1, 1)).full_simplify()
sigma_tot
""", {}),

("ee_check513", """
peskin_513 = 4*pi*a_em^2/(3*Ecm^2) * sqrt(1 - M_mu^2/E^2) * (1 + M_mu^2/(2*E^2))
(sigma_tot - peskin_513).canonicalize_radical()       # E > 0, so sqrt(E^2) = E
""", {}),

("ee_plot", """
from feynsage import plotting
x = var('x')                    # x = E / M_mu
ratio = (sigma_tot / (4*pi*a_em^2/(3*Ecm^2))).subs(E=x*M_mu).full_simplify()
fig = plotting.quick_plot([ratio], (x, 1, 4), labels=[r"$\\sigma / \\sigma_{\\rm massless}$"])
""", {"fig": "fig"}),

# ------------------------------------------------------------------ Compton scattering
("c_pol", """
kp, pp = momenta("kp pp", latex_name=["k'", "p'"])
e_in = polarization("e_in", k, latex_name=r"\\epsilon")       # photon in
e_out = polarization("e_out", kp, latex_name=r"\\epsilon'")   # photon out
eo = conjugate_vector(e_out)                  # eps*(k')
eo
""", {}),

("c_M", """
Mc = (ubar(pp, m) * slash(eo) * (slash(p) + slash(k) + m) * slash(e_in) * u(p, m) / (2*dot(p, k))
      + ubar(pp, m) * slash(e_in) * (slash(p) - slash(kp) + m) * slash(eo) * u(p, m) / (-2*dot(p, kp)))
Mc
""", {}),

("c_sum", """
r = spin_sum(Mc * conjugate(Mc), rules={pp: p + k - kp})
r = polarization_sum(r, [e_in, e_out]) / 4
onshell = {dot(p, p): m^2, dot(k, k): 0, dot(kp, kp): 0, dot(k, kp): dot(p, k) - dot(p, kp)}
compton = r.subs(onshell)
len(compton.operands())
""", {"text": True}),

("c_check587", """
pk, pkp = dot(p, k), dot(p, kp)
peskin_587 = 2*(pkp/pk + pk/pkp + 2*m^2*(1/pk - 1/pkp) + m^4*(1/pk - 1/pkp)^2)
(compton - peskin_587).simplify_rational()
""", {}),

("c_ward", """
r2 = polarization_sum(spin_sum(Mc * conjugate(Mc), rules={pp: p + k - kp}), [e_in, e_out],
                      gauge="physical", n=p) / 4
(r2.subs(onshell) - compton).simplify_rational()
""", {}),

("c_kn", """
w = var('w', latex_name=r'\\omega'); th = var('th', latex_name=r'\\theta')
wp = w / (1 + w/m*(1 - cos(th)))                       # omega', Peskin (5.89)
lab = compton.subs({pk: m*w, pkp: m*wp}).subs(e=1)
dsig = 1/(2*w) * 1/(2*m) * 1/(8*pi) * wp^2/(w*m) * (4*pi*a_em)^2 * lab
peskin_591 = pi*a_em^2/m^2 * (wp/w)^2 * (wp/w + w/wp - sin(th)^2)
(dsig - peskin_591).full_simplify()
""", {}),

("c_thomson", """
low = limit(peskin_591, w=0)
thomson = pi*a_em^2/m^2 * (1 + cos(th)^2)
(low - thomson).trig_simplify(), integrate(low * sin(th), th, 0, pi)
""", {}),

("c_plot", """
th_ = var('theta')
curves = [1 + cos(th_)^2] + [(peskin_591 / (pi*a_em^2/m^2)).subs(w=W*m).subs(th=th_) for W in (1/2, 2, 10)]
fig = plotting.quick_plot(curves, (th_, 0, pi),
                          labels=[r"$\\omega = 0$", r"$\\omega = m/2$", r"$\\omega = 2m$", r"$\\omega = 10m$"])
""", {"fig": "fig"}),

# ------------------------------------------------------------------ chiral currents
("vma", """
L = dirac_trace(slash(p2) * gamma(mu) * PL() * slash(p1) * gamma(nu) * PL())
L
""", {}),

("vma2", """
Lk = dirac_trace(slash(k2) * gamma(mu) * PL() * slash(k1) * gamma(nu) * PL())
contract(L * Lk)
""", {}),

("vma3", """
Lv = dirac_trace(slash(p2) * gamma(mu) * slash(p1) * gamma(nu)) / 2
contract(Lv * Lk)
""", {}),

# ------------------------------------------------------------------ colour
("col_idx", """
a, b, cc, dd = color_indices("a b cc dd")
i, j = quark_colors("i j")
color_factor(T_color(a, i, j) * T_color(b, j, i))
""", {}),

("col_cf", """
color_factor(color_chain(a, a, i=i, j=j))
""", {}),

("col_ca", """
color_factor(f_color(a, cc, dd) * f_color(b, cc, dd))
""", {}),

("col_n3", """
color_factor(color_trace(a, b, a, b), N=3)
""", {}),

("col_delta", """
color_factor(color_delta(i, j) * color_delta(j, i))
""", {}),

("col_R", """
Ncol = color_factor(color_delta(i, j) * color_delta(j, i), N=3)
charges = {'u': 2/3, 'd': -1/3, 's': -1/3, 'c': 2/3, 'b': -1/3}
R = Ncol * sum(Q^2 for Q in charges.values())
R
""", {}),

("col_gg", """
color_factor(color_chain(a, b, i=i, j=j) * color_chain(b, a, i=j, j=i))
""", {}),

# ------------------------------------------------------------------ massive vector boson
("B_amp", """
q = momenta("q")                                  # the B
g, MB = var('g M_B')
eB = polarization("eB", q, mass=MB, latex_name=r"\\varepsilon")
Ma = g * vbar(p2) * slash(conjugate_vector(eB)) * u(p1)
Ma
""", {}),

("B_sum", """
ra = polarization_sum(spin_sum(Ma * conjugate(Ma)), eB)
ra
""", {}),

("B_onshell", """
ra = ra.subs({dot(p1, p1): 0, dot(p2, p2): 0, dot(q, p1): dot(p1, p2), dot(q, p2): dot(p1, p2),
              dot(q, q): 2*dot(p1, p2)})
ra / 4, ra / 3
""", {}),

("B_rates", """
s = var('s')
sigma_B = 1/(2*s) * (ra/4).subs({dot(p1, p2): s/2}) * 2*pi           # times delta(s - M_B^2)
Gamma_B = 1/(2*MB) * (ra/3).subs({dot(p1, p2): MB^2/2}) / (8*pi)
sigma_B, Gamma_B, sigma_B - 12*pi^2/MB * Gamma_B
""", {}),

("B_b_amp", """
t = var('t'); uu = SR.var('uu', latex_name='u')      # not var('u'): u is the spinor
ea = polarization("ea", p1, latex_name=r"\\epsilon")                  # the photon
eb = polarization("eb", p2, mass=MB, latex_name=r"\\varepsilon")     # the B
ac, bc = conjugate_vector(ea), conjugate_vector(eb)
Mb = e*g*(vbar(k2) * slash(bc) * (slash(k1) - slash(p1)) * slash(ac) * u(k1) / (-2*dot(k1, p1))
          + vbar(k2) * slash(ac) * (slash(k1) - slash(p2)) * slash(bc) * u(k1) / (MB^2 - 2*dot(k1, p2)))
Mb
""", {}),

("B_b_sum", """
mandelstam = {dot(k1, k1): 0, dot(k2, k2): 0, dot(p1, p1): 0, dot(p2, p2): MB^2,
              dot(k1, k2): s/2, dot(k1, p1): -t/2, dot(k2, p2): (MB^2 - t)/2,
              dot(k1, p2): (MB^2 - uu)/2, dot(k2, p1): -uu/2, dot(p1, p2): (s - MB^2)/2}
rb = polarization_sum(spin_sum(Mb * conjugate(Mb)), [ea, eb]) / 4
rb = rb.subs(mandelstam).subs(s=MB^2 - t - uu).simplify_rational()
rb
""", {}),

("B_check587", """
x587 = 2*e^2*g^2*(uu/t + t/uu + 2*s*MB^2/(t*uu))
(rb - x587.subs(s=MB^2 - t - uu)).simplify_rational()
""", {}),

("B_ward", """
Mk = e*g*(vbar(k2) * slash(p2) * (slash(k1) - slash(p1)) * slash(ac) * u(k1) / (-2*dot(k1, p1))
          + vbar(k2) * slash(ac) * (slash(k1) - slash(p2)) * slash(p2) * u(k1) / (MB^2 - 2*dot(k1, p2)))
polarization_sum(spin_sum(Mk * conjugate(Mk)), ea).subs(mandelstam).subs(s=MB^2 - t - uu).simplify_rational()
""", {}),

("B_angle", """
ph = var('ph')                                     # th = theta, from the Compton part
assume(s > MB^2, MB > 0)
dOmega = 1/(2*s) * (s - MB^2)/(2*sqrt(s)) / (16*pi^2*sqrt(s)) * x587
angles = {t: (MB^2 - s)*sin(th/2)^2, uu: (MB^2 - s)*cos(th/2)^2}
dcos = 2*pi * dOmega.subs(angles).subs(e=sqrt(4*pi*a_em))
x591 = a_em*g^2*(1 - MB^2/s)/(2*s*sin(th)^2) * (1 + cos(th)^2 + 4*s*MB^2/(s - MB^2)^2)
(dcos - x591).subs(th=2*ph).expand_trig().simplify_full()
""", {}),

("B_massless", """
(x591.subs({MB: 0, g: sqrt(4*pi*a_em)}) - 2*pi*a_em^2/s * (1 + cos(th)^2)/sin(th)^2).simplify_full()
""", {}),

("B_plot", """
cc = var('cc')                                     # cos(theta)
shape = (x591 / x591.subs(th=pi/2)).subs(th=arccos(cc))
curves = [shape.subs({MB: sqrt(f_*s)}).simplify_full() for f_ in (0, 1/2, 9/10)]
fig = plotting.quick_plot(curves, (cc, -0.95, 0.95),
                          labels=[r"$M_B \\to 0$", r"$M_B^2 = s/2$", r"$M_B^2 = 9s/10$"])
""", {"fig": "fig"}),

("B_forward", """
c = var('c'); assume(c > -1, c < 1)
A = limit(((1 - c) * x591.subs(th=arccos(c))).simplify_full(), c=1)
A
""", {}),

("B_ww", """
x = var('x')
f = a_em/(2*pi) * (1 + (1 - x)^2)/x                 # f(x) without its log(s/m_e^2)
ww = (f * pi*g^2/s).subs(x=1 - MB^2/s)             # pi g^2 delta(M_B^2 - (1 - x) s), x integral done
(A - ww).simplify_full()
""", {}),

# ------------------------------------------------------------------ equivalent photons (P&S 6.2)
("epa_traces", """
e1, e2 = momenta("e1 e2", latex_name=[r"\\epsilon_1", r"\\epsilon_2"])
CC = dirac_trace(slash(pp) * slash(e1) * PR() * slash(p) * slash(e1))     # |C|^2
DD = dirac_trace(slash(pp) * slash(e2) * PR() * slash(p) * slash(e2))     # |D|^2
CD = dirac_trace(slash(pp) * slash(e1) * PR() * slash(p) * slash(e2))     # C D^*
CD
""", {}),

("epa_frame", """
Ep = SR.var('Ep', latex_name="E'")
qt = momenta("qt", latex_name=r"\\tilde{q}")
assume(E > Ep, Ep > 0, th > 0, th < pi)
vec = {p: vector([E, 0, 0, E]), pp: vector([Ep, Ep*sin(th), 0, Ep*cos(th)]),
       e1: vector([0, Ep*cos(th) - E, 0, -Ep*sin(th)]) / sqrt(E^2 + Ep^2 - 2*E*Ep*cos(th)),
       e2: vector([0, 0, 1, 0]),
       qt: vector([E - Ep, Ep*sin(th), 0, Ep*cos(th) - E])}              # q~ = (q0, -q)
mdot = lambda a, b: a[0]*b[0] - a[1]*b[1] - a[2]*b[2] - a[3]*b[3]
frame = {dot(a, b): mdot(vec[a], vec[b]) for a in vec for b in vec}
frame[epsilon(e1, e2, p, pp)] = -matrix([vec[e1], vec[e2], vec[p], vec[pp]]).det()   # eps^{0123} = +1
w0 = SR.wild(0)
root = (E^2 - 2*E*Ep + Ep^2)^w0 == (E - Ep)^(2*w0)      # sqrt((E - E')^2) = E - E' because E > E'
small = lambda X: X.subs(frame).series(th, 3).truncate().simplify_full().subs(root).factor()
[small(X) for X in (CC, DD, CD)]
""", {}),

("epa_619", """
C2 = E*Ep*(E + Ep)^2/(E - Ep)^2 * th^2
D2 = E*Ep * th^2
CDs = I*E*Ep*(E + Ep)/(E - Ep) * th^2
[(small(X) - Y).simplify_full() for X, Y in ((CC, C2), (DD, D2), (CD, CDs))]
""", {}),

("epa_cross", """
CDL = dirac_trace(slash(pp) * slash(e1) * PL() * slash(p) * slash(e2))
(CD + CDL).subs(frame).simplify_full()
""", {}),

("epa_B", """
J = ubar(pp) * (slash(p) - slash(pp)) * u(p)              # q.J with q = p - p'
Jt = ubar(pp) * slash(qt) * u(p)                          # q~.J
Q = vec[p] - vec[pp]
q2_, qqt = mdot(Q, Q), mdot(Q, vec[qt])
B2 = (spin_sum(Jt * conjugate(Jt)) / 2).subs(frame) * q2_^2/(q2_^2 - qqt^2)^2
spin_sum(J * conjugate(J)).subs(frame).simplify_full(), B2.series(th, 5).truncate().simplify_full()
""", {}),

("epa_avg", """
Mh = SR.var('Mh', latex_name=r'|\\widehat{\\mathcal{M}}|^2')
avg = (C2 * Mh/2 + D2 * Mh/2).subs(Ep=(1 - x)*E)
avg.factor()
""", {}),

("epa_q2", """
lt = var('lt')                       # theta and m both small: scale them together
Epx = (1 - x)*E
q2m = -2*(E*Epx - sqrt(E^2 - m^2)*sqrt(Epx^2 - m^2)*cos(th)) + 2*m^2
assume(x > 0, x < 1)
q2_small = q2m.subs(th=lt*th, m=lt*m).series(lt, 3).truncate().subs(lt=1)
q2_small = q2_small.subs((E^2*x^2 - 2*E^2*x + E^2)^w0 == (E*(1 - x))^(2*w0)).subs((E^2)^w0 == E^(2*w0))
q2_small.simplify_full(), (q2_small - (-(1 - x)*E^2*th^2 - x^2*m^2/(1 - x))).simplify_full()
""", {}),

("epa_theta", """
a2, tmax = var('a2 tmax'); assume(a2 > 0, tmax > 0)
I_th = integrate(th^3/(th^2 + a2)^2, th, 0, tmax)
I_th, (I_th + log(a2)/2 - log(tmax)).subs(a2=0).simplify_log()
""", {}),

("epa_N", """
Mt, L = var('M_t L')                            # target mass; L = log(E^2/m^2)
q2_0 = -(1 - x)*E^2*(th^2 + a2)
per = (Epx*E * th/(8*pi^2) * e^2/q2_0^2 * avg).simplify_full()   # times dx dtheta and the target phase space
sigma_gamma = (Mh/2) / (2*x*E * 2*Mt)          # photon cross section, per unit of target phase space
N = (per / th^3 * (th^2 + a2)^2 * L/2) / (2*E * 2*Mt) / sigma_gamma
N.subs(e=sqrt(4*pi*a_em)).simplify_full()
""", {}),

# ------------------------------------------------------------------ g - 2 (P&S 6.3)
("g2_higgs", """
from feynsage.dirac import projector, form_factor
q2 = SR.var('q2', latex_name='q^2'); mh = SR.var('mh', latex_name='m_h'); lam = SR.var('lam', latex_name=r'\\lambda')
P = projector('F2', 'mu', 'p1', 'm', 'p2', 'm')
kin = {'p1^2': 'm^2', 'p2^2': 'm^2', 'p1.p2': 'm^2 - q2/2'}
higgs = form_factor(['p2 - l + m', 'mu', 'p1 - l + m'], P, ['l', 'mh'], ['l - p2', 'm'], ['l - p1', 'm'], kin=kin)
Lh = higgs.coefficients()['1']
uv_part(Lh).simplify_full()
""", {}),

("g2_F2h", """
F2h = -lam^2/(32*pi^2) * finite_part(Lh, 1).subs(q2=0)
from feynsage import tex
LatexExpr(r"\\begin{aligned} \\delta F_2(0) &= -\\frac{\\lambda^2}{32\\pi^2}\\Big[ "
          + tex.lines(finite_part(Lh, 1).subs(q2=0).expand(), per_line=3) + r" \\Big] \\end{aligned}")
""", {}),

("g2_634", """
x = var('x')
matrix(SR, [[R, F2h.subs(m=1, mh=sqrt(R), lam=1).n(digits=60).n(digits=20),
              (integrate((1 - x)^2*(1 + x)/((1 - x)^2 + x*R), x, 0, 1)/(4*pi)^2).n(digits=60).n(digits=20)]
             for R in (1/2, 10, 1000)])                 # R = (m_h/m)^2, feynsage, eq. (6.34)
""", {}),

("g2_heavy", """
a, w = var('a w'); assume(a > 1, w > 0)
Ih = integrate(((1 - x)^2*(1 + x)/((x + a)*(x + 1/a))).partial_fraction(x), x, 0, 1)
R_ = 2 + a + 1/a                                  # (1 - x)^2 + x R_ = (x + a)(x + 1/a)
(Ih - (log(R_) - 7/6)/R_).subs(a=1/w).series(w, 3).truncate().simplify_full()
""", {}),

("g2_numbers", """
electron = F2h.subs(lam=3/10^6, m=511/1000, mh=60000)          # masses in MeV
muon_60 = F2h.subs(lam=6/10^4, m=10566/100, mh=60000)
muon_light = F2h.subs(lam=6/10^4, m=10566/100, mh=10566/100000)  # m_h = m_mu/1000
[v.n(digits=60).n(digits=4) for v in (electron, muon_60, muon_light)]
""", {}),

("g2_plot", """
import numpy as np, matplotlib.pyplot as plt
ratios = np.logspace(-3, 6, 37)                          # m_h / m
unit = [float(F2h.subs(lam=1, m=1, mh=QQ(float(rv))).n(digits=60).real()) for rv in ratios]
fig, ax = plt.subplots(figsize=(5.2, 3.2))
ax.loglog(ratios, [3e-6**2*u_ for u_ in unit], label=r"electron, $\\lambda = 3\\times10^{-6}$")
ax.loglog(ratios, [6e-4**2*u_ for u_ in unit], label=r"muon, $\\lambda = 6\\times10^{-4}$")
ax.axhline(1e-10, ls="--", color="C0"); ax.axhline(3e-8, ls="--", color="C1")
ax.set_xlabel(r"$m_h/m$"); ax.set_ylabel(r"$\\delta F_2(0)$"); ax.legend(fontsize=8)
fig.tight_layout()
""", {"fig": "fig"}),

("g2_axion", """
ma = SR.var('ma', latex_name='m_a')
axion = form_factor(['g5', 'p2 - l + m', 'mu', 'p1 - l + m', 'g5'], P, ['l', 'ma'], ['l - p2', 'm'], ['l - p1', 'm'], kin=kin)
F2a = lam^2/(32*pi^2) * finite_part(axion.coefficients()['1'], 1).subs(q2=0)
matrix(SR, [[R, F2a.subs(m=1, ma=sqrt(R), lam=1).n(digits=60).n(digits=20),
              (-integrate((1 - x)^3/((1 - x)^2 + x*R), x, 0, 1)/(4*pi)^2).n(digits=60).n(digits=20)]
             for R in (1/100, 2, 1000)])                # R = (m_a/m)^2, feynsage, eq. (6.40)
""", {}),

("g2_axion_light", """
Ia = integrate(((1 - x)^3/((x + a)*(x + 1/a))).partial_fraction(x), x, 0, 1)
light = -lam^2/(4*pi)^2 * integrate((1 - x)^3/(1 - x)^2, x, 0, 1)       # m_a -> 0
heavy = (Ia - (log(R_) - 11/6)/R_).subs(a=1/w).series(w, 3).truncate().simplify_full()
light, heavy, sqrt(10^-10 / (light/lam^2).abs()).n(digits=3)
""", {}),

("g2_axion_plot", """
ratios = np.logspace(-3, 4, 29)                          # m_a / m
lam_max = [float(sqrt(1e-10 / abs(F2a.subs(lam=1, m=1, ma=QQ(float(rv))).n(digits=60).real()))) for rv in ratios]
fig, ax = plt.subplots(figsize=(5.2, 3.2))
ax.loglog(ratios, lam_max, color="C3")
ax.fill_between(ratios, lam_max, 1, color="C3", alpha=0.15)
ax.text(ratios[1], 0.1, "excluded: $|\\\\delta F_2(0)| > 10^{-10}$", fontsize=8)
ax.set_xlabel(r"$m_a/m_e$"); ax.set_ylabel(r"largest allowed $\\lambda$"); ax.set_ylim(1e-5, 1)
fig.tight_layout()
""", {"fig": "fig"}),

# ------------------------------------------------------------------ conventions
("conv_note", """
set_convention("note")
spacetime_dimension()
""", {}),

("conv_std", """
set_convention("standard")
spacetime_dimension()
""", {}),

("conv_eu", """
set_convention(euclidean=True)
t5 = dirac_trace(gamma(mu) * gamma(nu) * gamma(rho) * gamma(tau) * gamma5())
set_convention(euclidean=False)
t5
""", {}),

("wick", """
to_euclidean(dot(p, p) - m^2)
""", {}),

# ------------------------------------------------------------------ speed and debugging
("speed", """
import time
exprs = [slash(p) * gamma(mu) * slash(q) * gamma(nu) * n for n in range(1, 201)]
t0 = time.time(); one_by_one = [dirac_trace(x) for x in exprs]; t1 = time.time()
together = dirac_trace(exprs); t2 = time.time()
print("200 traces one by one: %.3f s,  in one FORM run: %.3f s" % (t1 - t0, t2 - t1))
""", {}),

("debug", """
dirac_trace(slash(p) * gamma(mu) * slash(q) * gamma(mu), debug=True)
""", {}),

("errors", """
for bad in (lambda: slash(p + m),
            lambda: dirac_trace(gamma(var('muu')) * gamma(nu)),
            lambda: contract(gamma(mu) * gamma(mu) * gamma(mu))):
    try:
        bad()
    except Exception as err:
        print(type(err).__name__ + ":", str(err).splitlines()[0])
""", {}),
]
