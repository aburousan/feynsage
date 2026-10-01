# Tensor coefficients (PVC, PVD) with IR poles and at vanishing Gram determinants, against
# Package-X 2.1.1 LoopRefine (1/eps^2, 1/eps and eps^0 at mu = 1).
#  - IR: 8 triangle points x 8 coefficients, 7 box points x 6 coefficients.  At the box point
#    (0,-2,-3,-1; -4,-5; 0,0,0,1) Package-X's tensor values are wrong because its numerical
#    ScalarC0[-4,-3,-1,0,0,1] = +0.671 is wrong: all invariants are spacelike, so C0 must be negative,
#    and direct integration gives -0.58597680967236472265 = feynsage.  Independent checks (2026-10-01):
#    Mathematica NIntegrate of the Feynman-parameter integral -0.585976809672364722902; FeynCalc 10.2
#    FCFeynmanParametrize integrated -0.585976809672364725; Package-X itself is negative and smooth
#    nearby (m2 = 0.99, 1.01, s2 = -0.9, -1.1 agree with NIntegrate to 16 digits) and gives a third
#    value, -0.721, for the same integral relabelled as (-3,-1,-4; 1,0,0).  LoopTools 2.16 (both of its
#    C0 algorithms) fails here as well; see tests/looptools/README.md.  With that C0 put into
#    Package-X's own reduction its D1 becomes 0.0217878..., which is feynsage's value; we check that.
#  - Vanishing Gram (feynsage computes these from Feynman parameters): 11 points.
import sys, json, re, time; sys.path.insert(0, '.')
from feynsage import PVC, PVD, eps, mu, C0
import mpmath as mp
ok, t0 = True, time.time()

def coeffs(e):
    e = SR(e).subs(mu=1).expand()
    return [complex(e.coefficient(eps, j).n(prec=120)) for j in (-2, -1, 0)]

worst = {}
bad_px = [0, -2, -3, -1, -4, -5, 0, 0, 0, 1]
for r in json.load(open('tests/data/packagex_ir_tensor.json')):
    a = [QQ(x) for x in r['args']]
    if a == bad_px:
        continue
    got = coeffs((PVC if r['kind'] == 'C' else PVD)(*r['index'], *a))
    err = max(abs(g - complex(*w)) / max(1, abs(complex(*w))) for g, w in zip(got, r['packagex']))
    key = "%s%s" % (r['kind'], r['index'])
    worst[key] = max(worst.get(key, 0), err)
for k, v in worst.items():
    print("%-14s IR points, worst difference from Package-X %.1e" % (k, v))
ok &= max(worst.values()) < 1e-12

# the Package-X C0 bug and the corrected box value
c0 = complex(C0(-4, -3, -1, 0, 0, 1).n(prec=120))
with mp.workdps(25r):
    Y = {(0, 0): 0r, (1, 1): 0r, (2, 2): 2r, (0, 1): 4r, (1, 2): 4r, (0, 2): 2r}   # m^2_i + m^2_j - p_ij^2 (x2 on the diagonal)
    y = lambda i, j: Y[(min(i, j), max(i, j))]
    f = lambda x1, x2: 1r / (sum(y(i, j) * xx * yy for i, xx in enumerate([1r - x1 - x2, x1, x2])
                                for j, yy in enumerate([1r - x1 - x2, x1, x2])) / 2r)
    direct = -mp.quad(lambda x1: mp.quad(lambda x2: f(x1, x2), [0r, 1r - x1]), [0r, 1r])
print("C0(-4,-3,-1; 0,0,1): feynsage %.17f, direct integration %.17f, Package-X +0.67125 (wrong)" % (c0.real, float(direct)))
ok &= abs(c0 - complex(direct)) < 1e-15
d1 = coeffs(PVD(0, 1, 0, 0, *[QQ(x) for x in bad_px]))
px_fixed = -0.09714019792837501855 - 0.21428571428571428571 * complex(C0(-5, -3, -2, 0, 1, 0).n(prec=120))
print("D1 at that box point: feynsage %.15f, Package-X reduction with the correct C0 %.15f" % (d1[2].real, px_fixed.real))
ok &= abs(d1[2] - px_fixed) < 1e-15

w = 0
for r in json.load(open('tests/data/packagex_gram.json')):
    head = r['call']
    f = PVC if head.startswith('PVC') else PVD
    nums = [QQ(x) for x in re.findall(r"-?\d+(?:/\d+)?", head)]
    ni = 3 if f is PVC else 4
    got = coeffs(f(*[int(x) for x in nums[:ni]], *nums[ni:]))
    w = max(w, max(abs(g - complex(*x)) for g, x in zip(got, r['packagex'])))
print("vanishing Gram determinant: 11 points, worst difference from Package-X %.1e" % w)
ok &= w < 1e-14
m = var('m')
c1 = PVC(0, 1, 0, m ** 2, 0, m ** 2, 0, m, m)
c00 = PVC(1, 0, 0, m ** 2, 0, m ** 2, 0, m, m)
exact_ok = bool((c1 - 1 / (2 * m ** 2)).simplify_full() == 0) and \
    bool((c00 - (QQ(1) / 4 + (1 / eps + log(mu ** 2 / m ** 2)) / 4)).log_expand().expand() == 0)
print("g-2 at q^2 = 0, exact: C1 = %s, C00 = %s   (Package-X: 1/(2 m^2), 1/4 + (1/eps + log(mu^2/m^2))/4): %s"
      % (c1.simplify_full(), c00.log_expand().expand(), exact_ok))
ok &= exact_ok
# lines l, l+p, l+2p with p^2 = 0 (momenta proportional, Gram zero): the shift is (x1 + 2 x2) p,
# so the p^mu coefficient is Int (x1 + 2 x2)/m^2 = 1/(2 m^2)  (found by the independent audit)
from feynsage import loop
cp = loop("l^mu", ["l", "m"], ["l+p", "m"], ["l+2*p", "m"], kin={"p^2": 0})["p^mu"]
print("l, l+p, l+2p at p^2 = 0: p^mu coefficient %s (expected 1/(2 m^2))" % cp)
ok &= bool((cp - 1 / (2 * m ** 2)).simplify_full() == 0)
print("time %.1fs" % (time.time() - t0))
print("ALL IR-TENSOR AND GRAM CHECKS PASS" if ok else "SOME CHECK FAILED")
