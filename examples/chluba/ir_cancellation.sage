# The infrared divergence of double Compton scattering (Chluba's thesis, sec. 4.4.5, cut off there
# at a lowest frequency) cancels against the one-loop virtual correction, here in dimensional
# regularisation.  The soft-photon factor depends only on the charged legs (the electron, P -> P'),
# so its virtual partner is the IR pole of the electron vertex at t = (P - P')^2 (Yennie-Frautschi-Suura).
#   real, one soft photon with energy < Delta E:  dsigma_C * ( -alpha I(t) / (2 pi eps) + finite )
#   virtual:                                       dsigma_C * 2 Re F1(t),  F1 = (alpha/4 pi) [A + B]
#   I(t) = 2 P.P' Int_0^1 dx / (m^2 - x(1-x) t) - 2     (angular integral of the eikonal factor)
# so the IR pole of A + B must be  I(t)/eps  (up to the t-independent UV pole), and
# F2 = -(alpha/4 pi) B -> alpha/(2 pi) as t -> 0 checks the normalisation.
import sys, json, time; sys.path.insert(0, '.')
from feynsage import form, loop, eps, mu
import mpmath as mp
O = 'examples/chluba/out/'
D_ = SR.var('D')
m = var('m')

def trace(factors):
    return form.dirac_trace(factors, vectors=['p', 'pp', 'l'], dim='D').subs(m=1)

# projections: Gamma^mu = A gamma^mu + B (P + P')^mu / (2 m); rows X = gamma_mu and X = (P + P')_mu
num_g = trace(['p + m', 'mu', 'pp + m', 'nu', 'pp + l + m', 'mu', 'p + l + m', 'nu'])
num_P = trace(['p + m', 'pp + m', 'nu', 'pp + l + m', 'p + pp', 'p + l + m', 'nu'])
dot_ = function('dot')
def tree(factors, tval):
    tab = {('p', 'p'): 1, ('pp', 'pp'): 1, ('p', 'pp'): 1 - tval / 2}
    e = form.dirac_trace(factors, vectors=['p', 'pp'], dim='D').subs(m=1)
    return e.substitute_function(dot_, lambda a, b: tab.get((str(a), str(b)), tab.get((str(b), str(a))))).expand()

res = {}
def form_factors(tval):
    kin = {"p^2": 1, "pp^2": 1, "p.pp": 1 - tval / 2}
    vs = []
    for num in (num_g, num_P):
        r = loop(form.to_loop(num).replace('D', 'd'), ["l", "0"], ["l + pp", "1"], ["l + p", "1"], kin=kin)
        vs.append(SR(r["1"]))
    PP2 = 4 - tval                                   # (P + P')^2
    Mx = matrix(SR, [[tree(['p + m', 'mu', 'pp + m', 'mu'], tval), tree(['p + m', 'p + pp', 'pp + m'], tval) / 2],
                     [tree(['p + m', 'pp + m', 'p + pp'], tval), tree(['p + m', 'pp + m'], tval) * PP2 / 2]])
    Mx = Mx.subs({D_: 4 - 2 * eps})
    A, B = Mx.solve_right(vector(SR, vs))
    F1 = (A + B).series(eps, 1).truncate().expand()
    F2 = (-B).series(eps, 1).truncate().expand()
    return F1, F2

def I_ang(tval):
    tv = mp.mpf(float(tval))
    return 2 * (1 - tv / 2) * mp.quad(lambda x: 1 / (1 - x * (1 - x) * tv), [0, 1]) - 2

t0 = time.time()
ts = [QQ(-1) / 2, -1, -3, -8]
poles = {}
for tv in ts:
    F1, F2 = form_factors(tv)
    poles[tv] = F1.subs(mu=1).coefficient(eps, -1)
    print("t = %5s:  1/eps pole of A + B = %s" % (tv, poles[tv].n(digits=15)))
lines = []
for tv in ts[1:]:
    got = (poles[tv] - poles[ts[0]]).n(digits=15)
    want = I_ang(tv) - I_ang(ts[0])
    lines.append((str(tv), float(got.real()), float(want)))
    print("t = %5s:  pole(t) - pole(-1/2) = %.12f   I(t) - I(-1/2) = %.12f   difference %.1e" % (tv, got.real(), want, abs(got.real() - float(want))))
F1s, F2s = form_factors(-QQ(1) / 10 ** 6)
f2 = F2s.subs(mu=1).coefficient(eps, 0).n(digits=12)
print("F2 at t = -1e-6, in units alpha/(4 pi): %s  (expected 2, i.e. F2(0) = alpha/2pi)" % f2)
res['ir_lines'] = lines; res['F2_small_t'] = str(f2); res['time'] = time.time() - t0
json.dump(res, open(O + 'chluba_ir.json', 'w'), indent=1)
print("time %.1f s" % (time.time() - t0))

# ---- Lightman's soft-photon emission law (thesis eq. 4.24) from the same eikonal integral.
# Soft photons per Compton scattering: dN = (alpha/pi) I(t) domega2/omega2.  Cold electrons and soft
# incident photons: -t = 2 omega0^2 (1 - cos theta) small, averaged over the Thomson distribution.
x_, t_, c_, w0 = var('x_ t_ c_ w0')
Iser = (2 * (1 - t_ / 2) * integrate((1 / (1 - x_ * (1 - x_) * t_)).series(t_, 3).truncate(), x_, 0, 1) - 2).series(t_, 2).truncate()
thomson = (1 + c_ ** 2)                                   # dsigma/dOmega up to constants
avg = integrate(Iser.subs(t_=-2 * w0 ** 2 * (1 - c_)) * thomson, c_, -1, 1) / integrate(thomson, c_, -1, 1)
avg = avg.series(w0, 3).truncate()
print("I(t) for small t:", Iser, "   Thomson average:", avg, "   so dN/domega2 = (alpha/pi) * (%s) / omega2" % avg)
res['lightman_I'] = str(Iser); res['lightman_avg'] = str(avg)
json.dump(res, open(O + 'chluba_ir.json', 'w'), indent=1)
