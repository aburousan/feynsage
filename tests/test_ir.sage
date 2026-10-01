# IR-divergent C0 and D0 (feynsage.ir, Ellis-Zanderighi) against Package-X 2.1.1 LoopRefine:
# the 1/eps^2, 1/eps and eps^0 coefficients at mu = 1, for all 6 triangles and 16 boxes, each also
# at a random relabelling.  Package-X cannot evaluate the finite parts of boxes 15 and 16
# (ScalarD0IR16 stays unevaluated), nor the threshold points; those are checked against a
# photon-mass regularised D0 instead: D0(lambda) = A log(lambda^2) + B + O(lambda) and the
# replacement log(lambda^2) -> 1/eps + log(mu^2) gives c_-1 = A, c_0 = B (soft singularity only).
import sys, json, time; sys.path.insert(0, '.')
import mpmath as mp
from feynsage import ir
from feynsage import scalar as S
ok, t0 = True, time.time()
rows = json.load(open('tests/data/packagex_ir.json'))
worst, nref = {}, 0
for r in rows:
    a = [QQ(x) for x in r['args']]
    key = "%s%d" % (r['kind'], r['ez'])
    cls = (ir.classify_c0 if r['kind'] == 'C' else ir.classify_d0)(*a)
    if cls is None or cls[0] != r['ez']:
        print("classification wrong:", key, a, cls); ok = False; continue
    if r['packagex'] is None:
        continue
    got = ir.ir_value(r['kind'], a)
    want = [complex(*w) for w in r['packagex']]
    err = max(abs(g - w) / max(1, abs(w)) for g, w in zip(got, want))
    worst[key] = max(worst.get(key, 0), err)
    nref += 1
for k in sorted(worst, key=lambda k: (k[0], int(k[1:]))):
    print("%-4s worst difference from Package-X %.1e" % (k, worst[k]))
ok &= max(worst.values()) < 1e-12
print("%d points compared with Package-X" % nref)

# triangle 5 (Package-X needs Analytic -> True there): C0(0, m^2, m^2; 0, 0, m) = 1/m^2 - (1/eps + log(mu^2/m^2))/(2 m^2)
c = ir.ir_value('C', [0, 4, 4, 0, 0, 2])
ok &= abs(c[1] + 1/8) < 1e-15 and abs(c[2] - (0.25 + float(log(4)) / 8)) < 1e-15
print("triangle 5 against Package-X (Analytic -> True): %s" % (abs(c[2] - (0.25 + float(log(4)) / 8)) < 1e-15))

# boxes 14-16 against a photon mass
def photon(a, lams=(mp.mpf('1e-9'), mp.mpf('1e-11'))):
    vals = []
    for lam in lams:
        b = [x if not (i >= 6 and x == 0) else lam for i, x in enumerate(a)]
        vals.append(complex(S._d0_terms(*[S.V(x, None) for x in b], exact=False).n))
    A = (vals[0] - vals[1]) / complex(2 * mp.log(lams[0]) - 2 * mp.log(lams[1]))
    return A, vals[1] - A * complex(2 * mp.log(lams[1]))
wp = 0
S.DPS, S.ED, S.EDM = 120r, mp.mpf('1e-75'), mp.mpf('1e-95')   # room for a tiny mass
with mp.workdps(120r):
    for r in rows:
        if r['kind'] != 'D' or r['ez'] not in (14, 15, 16):
            continue
        a = [QQ(x) for x in r['args']]
        c = ir.ir_value('D', a)
        A, B = photon([mp.mpf(int(x.numerator())) / int(x.denominator()) for x in a])
        wp = max(wp, abs(c[1] - A), abs(c[2] - B))
print("boxes 14-16 against photon-mass regularisation (lambda = 1e-11): worst %.1e (O(lambda) expected)" % wp)
ok &= wp < 1e-8

# boxes 3 and 5 on s12 s23 = p2 p4 (EZ 4.21, 4.25): continuity from neighbouring points
wc = 0
for a in [(0, 2, 0, -3, -2, 3, 0, 0, 0, 0), (0, 2, 5, -3, -2, 3, 0, 0, 0, 0)]:
    a = [QQ(x) for x in a]
    at = ir.ir_value('D', a)
    b = list(a); b[5] += QQ(1) / 10 ** 9
    near = ir.ir_value('D', b)
    wc = max(wc, abs(at[1] - near[1]), abs(at[2] - near[2]))
print("boxes 3, 5 at the vanishing overall denominator: continuity %.1e (O(1e-9) expected)" % wc)
ok &= wc < 1e-8
print("time %.1fs" % (time.time() - t0))
print("ALL IR CHECKS PASS" if ok else "SOME IR CHECK FAILED")
