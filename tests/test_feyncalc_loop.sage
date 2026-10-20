# One-loop results of feynsage against FeynArts + FeynCalc + Package-X (FeynHelpers): poles and finite parts.
# The FeynCalc numbers were made with tests/feyncalc_ref/virtual.wl (TID, then PaXEvaluate with
# PaXImplicitPrefactor -> 1).  FeynCalc's loop measure differs from feynsage's by the factor i 16 pi^4 and
# Package-X's scheme by c = -gamma_E - log(pi) per 1/eps:  FeynCalc = i 16 pi^4 (feynsage + c eps feynsage_pole).
#     sage tests/test_feyncalc_loop.sage
import time
from feynsage import *
from feynsage.models import EL, ME, MM, ML, GS, MB, MH, MW, SW, CW, MZ, SM
from feynsage.qft import tensors as T
from feynsage.qft.color import f_color, color_factor, color_indices
t_start = time.time()
FAIL = []
def check(name, ok):
    print(("ok   " if ok else "FAIL ") + name)
    if not ok:
        FAIL.append(name)
kappa = I*16*pi^4; c = -euler_gamma - log(pi)
def agree(V, pt, fc1, fc0, tol=1e-12):
    p1 = CC(SR(pole_parts(V)[1]).subs(pt).n(100)); p0 = CC(SR(finite_part(V, 1)).subs(pt).n(100))
    a1 = CC((kappa*p1).n(100)); a0 = CC((kappa*(p0 + c*p1)).n(100))
    return abs(a1/CC(fc1) - 1) < tol and abs(a0/CC(fc0) - 1) < tol

# ---------------------------------------------------------------- e- e+ -> mu- mu+ in QED, masses kept
P = process("e- e+ -> mu- mu+", model="QED")
FC = {'vacuum polarization (e, mu, tau loops)': ([5, 6, 7], [
          (-200.06784743719362071335161533526722368074*I, 619.05138362145934536158964617658147586696 + 951.266865875816336575640410192588905442*I),
          (-177.25809504356488079426745835777567439023*I, 553.98050401877602481258224833759074933422 + 990.85661795033123037053240296270860676094*I)]),
      'vertex corrections (IR divergent)': ([3, 4], [
          (629.38160901957578959192769605763632168744 + 911.38453568446903190666831642430014512032*I, -3277.14468801132338390883578017958229471654 - 1463.21335657258597245357357213625884613671*I),
          (557.12894815317122040493877913079721653074 + 900.70532026623221497141598920074978826793*I, -3300.26846277094335983703224375024160145794 - 1983.70242408116354429062097010042040047247*I)]),
      'boxes (IR divergent)': ([1, 2], [
          (-129.61928091206705389096178552882983878278*I, 30.51927356226111452769876359371407771533 + 717.1567777126659797339532296254961889909*I),
          (-53.32409276212250152618043325151452109942*I, -6.10795077766278664634353498125407695041 + 334.31792874644840368860400646803560635221*I)])}
for name, (dg, vals) in FC.items():
    V = P.virtual(average=False, diagrams=dg)
    ok = all(agree(V, {P.s: sv, P.t: tv, EL: 1, ME: 1/2, MM: 1, ML: 3/2}, f1, f0)
             for (sv, tv), (f1, f0) in zip([(30, -9), (50, -20)], vals))
    check("e+e- -> mu+mu- one loop, sum M1 M0*, %s = FeynCalc at two points" % name, ok)

# ---------------------------------------------------------------- H -> b b~ with the gluon loop (colour, IR divergent)
H = process("H -> b b~")
idx = [k + 1 for k, d in enumerate(H.loop_diagrams()) if any(f[0] == 'g' for f in d.fields)]
V = H.virtual(average=False, diagrams=idx)
check("H -> b b~, QCD vertex correction, sum M1 M0* = FeynCalc",
      agree(V, {MH: 125, MB: 5, MW: 80, SW: 1/2, EL: 1/3, GS: 1},
            3342.88143025970956497460040541199262233378 + 8974.9528265304294908239211403194309927143*I,
            -37994.29414928721361255247613418727659428998 - 56108.82864210761150812928381528385710726108*I))

# ---------------------------------------------------------------- Z -> nu nu~ with every electroweak loop
Z = process("Z -> nu_e nu_e~", exclude_fields=['g', 'ug'])
V = Z.virtual(average=False)
check("Z -> nu_e nu_e~, all seven electroweak one-loop diagrams (W, Z, G+-, e), sum M1 M0* = FeynCalc",
      agree(V, {MZ: 91, MW: 80, CW: 80/91, SW: sqrt(1 - (80/91)^2), ME: 1/2, EL: 1/3, MH: 125},
            29974.97383356800849970537448388649175572252*I,
            291.38703292824653146033768208175878954099 - 358688.0337110493866668715744325876765375835*I))

# ---------------------------------------------------------------- the ghost-gluon vertex in R_xi gauge, xi = 3
sm = SM()
ALL = ['A','Z','W+','H','G0','G+','u+','u-','uA','uZ','g','ug','e','mu','ta','ne','nm','nt','u','c','t','d','s','b']
only = lambda keep: [f for f in ALL if f not in keep]
p2, p3 = momenta("p2 p3"); a1, a2, a3 = color_indices("a1 a2 a3")
d = insert_fields(topologies(1, 1, 2, exclude=('tadpoles', 'wf')), ['ug'], ['ug', 'g'], sm, exclude_fields=only(('g', 'ug')))
M = sum(d.loop_amplitude([p2 + p3], [p2, p3], lambda n, e: n*T.comp(p3, e[2]), {"p2^2": 2, "p3^2": 3, "p2.p3": "-1/2"}, xi=3))
M = color_factor(M*f_color(a1, a2, a3), N=3).subs({GS: 1})
check("ghost-gluon vertex at xi = 3 (gluon propagator with 1/k^4), projected on the gluon momentum = FeynCalc",
      agree(M, {}, -532.95863765882536541706251399331216130694,
            1236.42732437085308388880613800105358716556 - 1674.33894073619028947572101362347534092017*I))

# ---------------------------------------------------------------- exact, in d dimensions, before any eps expansion
# FeynCalc's TID result (virtual*.wl with SYMBOLIC=1, pave_<tag>.m.gz): coefficients times A0, B0, C0, D0, with D,
# s, t and the masses symbols, against feynsage's virtual(masters=True).  Both sets of masters are put in the
# canonical form of feynsage.pv.pave(); FeynCalc's TID writes the integral as i pi^2 * PaVe, feynsage's M has
# 1/(16 pi^2) * master, so FeynCalc = i 16 pi^4 feynsage term by term.  No scheme constant enters: nothing is
# expanded in eps.  The difference is reduced over Q with i^2 = -1 and s_W^2 + c_W^2 = 1 (after M_W = M_Z c_W).
import os, gzip, re
from feynsage.pv import pave, PaVeA0, PaVeB0, PaVeC0, PaVeD0
REF = os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), 'feyncalc_ref')
def read_pave(path):
    txt = gzip.open(path, 'rt').read()
    txt = txt.replace('Pi', 'pi').replace('[', '(').replace(']', ')').replace('^', '**')
    fA = lambda a: pave('A', [a]); fB = lambda *a: pave('B', list(a))
    fC = lambda *a: pave('C', list(a)); fD = lambda *a: pave('D', list(a))
    names = {n: SR.var(n) for n in ['s', 't', 'ME', 'MM', 'ML', 'EL', 'MH', 'MB', 'MW', 'MZ', 'SW', 'CW', 'GS']}
    names.update({'D': var('d'), 'A0': fA, 'B0': fB, 'C0': fC, 'D0': fD, 'I': I, 'pi': pi, 'Sqrt': sqrt})
    return SR(sage_eval(txt, locals=names))
def same_in_d(a, b):
    """a - b = 0 as a function: masters become symbols, then exact rational-function arithmetic over Q."""
    e = (SR(a) - SR(b)).subs({MW: MZ*CW})
    ms = sorted(set(e.find(PaVeA0(SR.wild(0)))) | set(e.find(PaVeB0(*[SR.wild(i) for i in range(3)]))) |
                set(e.find(PaVeC0(*[SR.wild(i) for i in range(6)]))) | set(e.find(PaVeD0(*[SR.wild(i) for i in range(10)]))), key=str)
    sym = {m: SR.var('fsM%d' % k) for k, m in enumerate(ms)}
    e = e.subs(sym).subs({I: SR.var('fsI'), pi: SR.var('fsPi')})
    gens = sorted(set(str(v) for v in e.variables()) | {'fsI', 'SW', 'CW'})
    R = PolynomialRing(QQ, gens)
    J = R.ideal([R('fsI^2 + 1'), R('SW^2 + CW^2 - 1')])
    F = R.fraction_field()(str(e).replace('^', '**'))
    return J.reduce(F.numerator()) == 0, len(ms)
kappa16 = I*16*pi^4
SYM = [('ee_mumu_vp', "e- e+ -> mu- mu+", "QED", [5, 6, 7], 'vacuum polarization'),
       ('ee_mumu_vertex', "e- e+ -> mu- mu+", "QED", [3, 4], 'vertex corrections'),
       ('ee_mumu_box', "e- e+ -> mu- mu+", "QED", [1, 2], 'boxes'),
       ('ee_aa', "e- e+ -> gamma gamma", "QED", None, 'all diagrams (e, mu, tau loops, IR-divergent boxes)'),
       ('hbb', "H -> b b~", None, None, 'QCD vertex'),
       ('z_nn', "Z -> nu_e nu_e~", None, None, 'all electroweak loops')]
for tag, text, model, dg, what in SYM:
    path = os.path.join(REF, 'pave_%s.m.gz' % tag)
    if not os.path.exists(path):
        print("MISS %s (no FeynCalc reduction)" % tag); continue
    P = process(text, model=model) if model else process(text, exclude_fields=['g', 'ug'] if tag == 'z_nn' else ())
    if tag == 'hbb':                                       # the gluon loop only, as above
        dg = [k + 1 for k, d_ in enumerate(P.loop_diagrams()) if any(f[0] == 'g' for f in d_.fields)]
    fs = P.virtual(average=False, diagrams=dg, masters=True, polarization_gauge='feynman')   # FeynCalc's -g
    ok, nm = same_in_d(read_pave(path), kappa16*fs)
    check("%s one loop, %s: FeynCalc's TID reduction = feynsage's, exactly in d (%d masters)" % (text, what, nm), ok)

print("time %.0f s" % (time.time() - t_start))
print("ALL FEYNCALC LOOP CHECKS PASS" if not FAIL else "FAILED: %s" % FAIL)
