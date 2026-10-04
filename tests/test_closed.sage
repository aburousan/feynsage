# Symbolic closed forms of C0 and D0 (Package-X's C0Expand, D0Expand) and products of fermion lines
# (FermionLineProduct).  Reference values: feynsage's numerical closed forms (themselves checked
# against Package-X and LoopTools), Package-X 2.1.1 numbers, and explicit 4x4 Dirac matrices.
#     sage tests/test_closed.sage
import os, sys, time, random
sys.path.insert(0, os.path.abspath('.'))
from feynsage import *
from feynsage.scalar import c0_value, d0_value
from feynsage.ir import classify_d0

t_start = time.time()
fails = []


def same(name, flag):
    print("%-62s %s" % (name, "ok" if flag else "FAIL"))
    if not flag:
        fails.append(name)


def close(name, a, b, tol=1e-15):
    a, b = complex(a), complex(b)
    err = abs(a - b) / max(abs(b), 1e-300)
    same("%s (rel %.0e)" % (name, err), err < tol)


s1, s2, s3, s4, s12, s23, m0, m1, m2, m3, M = var('s1 s2 s3 s4 s12 s23 m0 m1 m2 m3 M')
# 1. C0 with symbols, against the numerical closed form (and Package-X at two points)
C = c0_expand(s1, s12, s2, m0, m1, m2)
same("C0 general: one condition, lambda > 0", len(C.conditions) == 1)
close("C0(-1, 2, -3; 1, 2, 3)", C.n(s1=-1, s12=2, s2=-3, m0=1, m1=2, m2=3, prec=80), c0_value(-1, 2, -3, 1, 2, 3))
close("C0 on a threshold, s1 = (m0 + m1)^2", C.n(s1=9, s12=-2, s2=1, m0=1, m1=2, m2=2, prec=80), c0_value(9, -2, 1, 1, 2, 2))
close("C0 above all thresholds", C.n(s1=30, s12=200, s2=40, m0=1, m1=2, m2=3, prec=80), c0_value(30, 200, 40, 1, 2, 3))
close("C0(-2/3, 12, -2/3; 4, 9/2, 3) = Package-X C0Expand", C.n(s1=-QQ(2)/3, s12=12, s2=-QQ(2)/3, m0=4, m1=QQ(9)/2, m2=3, prec=80),
      -0.036437555935892087822)
for name, args in (("C0 with m0 = 0", (s1, s12, s2, 0, m1, m2)), ("C0 massless", (s1, s12, s2, 0, 0, 0)),
                   ("C0 with s12 = 0", (s1, 0, s2, m0, m1, m2)), ("C0 equal masses", (s1, s12, s2, M, M, M))):
    Cx = c0_expand(*args)
    pt = {s1: -QQ(7)/3, s12: QQ(5)/2, s2: -11, m0: QQ(3)/2, m1: 2, m2: QQ(5)/3, M: QQ(4)/3}
    num = [SR(a).subs(pt) for a in args]
    close(name, Cx.n(pt, prec=80), c0_value(*num))
same("C0 outside lambda > 0 is refused", C.holds(s1=1, s12=1, s2=1, m0=1, m1=1, m2=1) is False)

# 2. D0 with symbols
random.seed(int(7))
R = lambda a, b, dd: QQ(random.randint(a, b)) / random.randint(1, dd)
for name, ms, npts in (("D0 all masses (pair 1,3)", (m0, m1, m2, m3), 6), ("D0 m2 = 0", (m0, m1, 0, m3), 4),
                       ("D0 m0 = m1 = 0", (0, 0, m2, m3), 4), ("D0 three massless", (0, 0, 0, m3), 4),
                       ("D0 massless", (0, 0, 0, 0), 4)):
    Dx = d0_expand(s1, s2, s3, s4, s12, s23, *ms)
    n = worst = 0
    while n < npts:
        pt = {x: R(-40, 40, 4) for x in (s1, s2, s3, s4, s12, s23)}
        pt.update({x: R(1, 16, 3) for x in (m0, m1, m2, m3)})
        num = [SR(a).subs(pt) for a in (s1, s2, s3, s4, s12, s23) + ms]
        if 0 in num[:6] or Dx.holds(pt) is not True:
            continue
        try:
            if classify_d0(*num) is not None:
                continue
            ref = d0_value(*num)
        except (ZeroDivisionError, ValueError, ArithmeticError):
            continue
        a = complex(Dx.n(pt, prec=80))
        worst = max(worst, abs(a - ref) / max(abs(ref), 1e-300))
        n += 1
    same("%s at %d random points (worst %.0e)" % (name, n, worst), worst < 1e-15)
Dm = d0_expand(s1, s2, s3, s4, s12, s23, m0, m1, m2, m3)
close("D0(-5/3, 11/3, -24, -21; -1, -10; 3, 7, 2, 4) = Package-X D0Expand",
      Dm.n(s1=-QQ(5)/3, s2=QQ(11)/3, s3=-24, s4=-21, s12=-1, s23=-10, m0=3, m1=7, m2=2, m3=4, prec=80),
      0.000497520392345595613177707660252)
D1 = d0_expand(s1, s2, s3, s4, s12, s23, 0, m1, m2, m3)
close("D0(28, -17, -3, -19/3; -3, 20; 0, 9/2, 8, 2) = Package-X D0Expand",
      D1.n(s1=28, s2=-17, s3=-3, s4=-QQ(19)/3, s12=-3, s23=20, m1=QQ(9)/2, m2=8, m3=2, prec=80),
      complex(0.00208941422854927396499186918102, 0.0011875166947368114299763868379))
same("D0 inline() gives the same number", abs(complex(Dm.inline().n(s1=-1, s2=-2, s3=-3, s4=-4, s12=-5, s23=-6, m0=1, m1=2, m2=3, m3=4, prec=80))
                                              - complex(Dm.n(s1=-1, s2=-2, s3=-3, s4=-4, s12=-5, s23=-6, m0=1, m1=2, m2=3, m3=4, prec=80))) < 1e-15)

# regressions found by review: the labelling rule of Denner (4.43) and a leading Landau singularity
close("D0(10, 50, -5, 13; 20/3, -1; 1, 2, 3, 4) = Package-X (both r real, k02 k13 < 0)",
      d0_value(10, 50, -5, 13, QQ(20)/3, -1, 1, 2, 3, 4), complex(-0.0117571737205770615801951, 0.0145422137675868033794622))
for pr in ((0, 1), (2, 3), (1, 2)):
    close("  symbolic D0 there with pair %s" % (pr,),
          d0_expand(s1, s2, s3, s4, s12, s23, m0, m1, m2, m3, pair=pr).n(s1=10, s2=50, s3=-5, s4=13, s12=QQ(20)/3, s23=-1, m0=1, m1=2, m2=3, m3=4, prec=80),
          complex(-0.0117571737205770615801951, 0.0145422137675868033794622))
try:
    d0_value(QQ(2001)/1000, 4, 0, 4, 0, 4, 1, 1, 1, 1)
    same("D0 at a vanishing Cayley determinant is refused", False)
except ZeroDivisionError:
    same("D0 at a vanishing Cayley determinant is refused", True)

# 3. expand_c0d0 on a loop() result
r = loop("1", ["l", "m"], ["l + p1", "m"], ["l + p1 + p2", "m"], kin={"p1^2": "s1", "p2^2": "s2", "p1.p2": "x"}).coefficients()['1']
E = expand_c0d0(r)
pt = dict(s1=-1, s2=-2, x=3, m=1)
close("expand_c0d0(loop triangle) = the triangle", E.n(**pt, prec=80), r.subs(**pt).n(prec=80))

# 4. products of fermion lines (FermionLineProduct)
U3, U1, U4, U2 = ("u", "p3", "m"), ("u", "p1", "m"), ("u", "p4", "m"), ("u", "p2", "m")
d = var('d')
def coeffs(f1, f2, **k):
    return {str(key): SR(c) for key, c in line_product((U3, f1, U1), (U4, f2, U2), **k).items()}
V = "(((('i', 'fs1'),), '1'), ((('i', 'fs1'),), '1'))"
A = "(((('i', 'fs1'),), '5'), ((('i', 'fs1'),), '5'))"
c = coeffs(["mu", "nu", "rho"], ["mu", "nu", "rho"])
same("g^mu g^nu g^rho (x) g_mu g_nu g_rho = 10 V(x)V + 6 A(x)A at d = 4", c[V].subs(d=4) == 10 and c[A].subs(d=4) == 6)
c = coeffs(["mu", "nu", "rho"], ["mu", "rho", "nu"])
same("... (x) g_mu g_rho g_nu = (2 - d) V(x)V - (d-1)(d-2)(d-3) A(x)A (= Package-X)",
     bool((c[V] - (2 - d)).expand() == 0) and bool((c[A] + (d - 1) * (d - 2) * (d - 3)).expand() == 0))
c = coeffs(["mu", "nu", "rho", "sigma"], ["mu", "nu", "rho", "sigma"])
same("four gammas: 1(x)1 = 3d^2 - 2d (40 at d = 4; Package-X: 27d^2 - 18d, wrong)",
     bool((c["(((), '1'), ((), '1'))"] - (3 * d ** 2 - 2 * d)).expand() == 0))
c = coeffs(["mu", "nu", "rho", "PL"], ["rho", "nu", "mu", "PL"])
key = "(((('i', 'fs1'),), 'L'), ((('i', 'fs1'),), 'L'))"
same("g^mu g^nu g^rho PL (x) g_rho g_nu g_mu PL = 4 V_L(x)V_L at d = 4", c[key].subs(d=4) == 4)

print("time %.1f s" % (time.time() - t_start))
print("ALL CLOSED-FORM AND LINE-PRODUCT CHECKS PASS" if not fails else "FAILED: %s" % fails)
