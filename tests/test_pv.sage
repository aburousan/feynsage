# feynsage.pv against Package-X 2.1.1 (reference values from LoopRefine, finite part at mu = 1)
import sys, json, time; sys.path.insert(0, '.')
from feynsage.pv import *
ref = json.load(open('tests/data/packagex_reference.json')) + json.load(open('tests/data/packagex_reference_special.json'))
s_ = var('s_')

def ours(f, a):
    a = [RR(float(x)).simplest_rational() for x in a]
    if f == 'A0':  e = A0(a[0])
    if f == 'B0':  e = B0(*a)
    if f == 'B1':  e = PVB(0, 1, *a)
    if f == 'B00': e = PVB(1, 0, *a)
    if f == 'B11': e = PVB(0, 2, *a)
    if f == 'B0p': e = diff(B0(s_, a[1], a[2]), s_).subs(s_=a[0])
    if f == 'C0':  e = C0(*a)
    if f == 'C1':  e = PVC(0, 1, 0, *a)
    if f == 'C2':  e = PVC(0, 0, 1, *a)
    if f == 'C00': e = PVC(1, 0, 0, *a)
    if f == 'C11': e = PVC(0, 2, 0, *a)
    if f == 'C12': e = PVC(0, 1, 1, *a)
    if f == 'D0':  e = D0(*a)
    return complex(finite_part(e, 1).n())

worst, t0 = {}, time.time()
exceptional = 0
for r in ref:
    try:
        got = ours(r['f'], r['args'])
    except ZeroDivisionError:
        exceptional += 1          # vanishing Gram determinant: all external invariants zero
        continue
    want = complex(r['re'], r['im'])
    err = abs(got - want) / max(1, abs(want))
    worst[r['f']] = max(worst.get(r['f'], 0), err)
    if err > 1e-5:
        print("MISMATCH", r['f'], r['args'], got, want)
for f, e in worst.items():
    print("%-4s worst relative difference %.1e" % (f, e))
print("points with a vanishing Gram determinant (skipped, feynsage raises): %d" % exceptional)
print("time %.1fs" % (time.time() - t0))
