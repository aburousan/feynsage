# The IR-finite box D0 in closed form (feynsage.scalar: Denner 1993 eq. 4.43 for nonzero masses,
# Denner-Dittmaier 2010 eqs. 3.76-3.84 for vanishing masses) against
#   1. Package-X 2.1.1 at 150 points (120 with all masses nonzero, 30 with one to four zero masses);
#   2. itself under relabellings of the four lines (D0 is symmetric under all 24);
#   3. direct integration over Feynman parameters at the two points where Package-X is only
#      accurate to 1e-9 and 2e-8 (one external p^2 = 0; reference values from scipy nquad, 1e-14).
# Then d0_closed (the exact expression) is checked against the number at a few points.
import sys, json, time, itertools; sys.path.insert(0, '.')
from feynsage.scalar import d0_value, d0_closed
ok, t0 = True, time.time()
rows = json.load(open('tests/data/packagex_d0.json'))
direct = {10: 0.00607769685971634, 113: 0.0147906701148701}
worst, n = 0, 0
for i, r in enumerate(rows):
    if r['packagex'] is None:
        continue
    a = [QQ(x) for x in r['args']]
    got = d0_value(*a)
    want = complex(*r['packagex'])
    err = abs(got - want) / abs(want)
    if i in direct:
        print("point %d: Package-X differs by %.1e; direct integration differs by %.1e" % (i, err, abs(got - direct[i]) / abs(got)))
        ok &= abs(got - direct[i]) / abs(got) < 1e-12
        continue
    worst = max(worst, err)
    n += 1
print("D0 against Package-X: worst relative difference %.1e over %d points" % (worst, n))
ok &= worst < 1e-9

def relabel(a, pm):
    P = {(0, 1): a[0], (1, 2): a[1], (2, 3): a[2], (0, 3): a[3], (0, 2): a[4], (1, 3): a[5]}
    pp = lambda i, j: P[(min(i, j), max(i, j))]
    m = a[6:]
    return [pp(pm[0], pm[1]), pp(pm[1], pm[2]), pp(pm[2], pm[3]), pp(pm[0], pm[3]), pp(pm[0], pm[2]), pp(pm[1], pm[3])] + [m[q] for q in pm]
spread = 0
for i in range(0, len(rows), 10):
    a = [QQ(x) for x in rows[i]['args']]
    v0 = d0_value(*a)
    for pm in itertools.permutations(range(4)):
        spread = max(spread, abs(d0_value(*relabel(a, pm)) - v0) / abs(v0))
print("D0 under the 24 relabellings (15 points): worst spread %.1e" % spread)
ok &= spread < 1e-25
wc = 0
for i in (0, 5, 60, 121, 130, 144):
    a = [QQ(x) for x in rows[i]['args']]
    e = d0_closed(*a)                       # raises if the exact expression and the number disagree
    wc = max(wc, abs(complex(e.n(prec=200)) - d0_value(*a)))
print("d0_closed: exact expression = number to %.1e at 6 points" % wc)
print("time %.1fs" % (time.time() - t0))
print("ALL D0 CHECKS PASS" if ok else "SOME D0 CHECK FAILED")
