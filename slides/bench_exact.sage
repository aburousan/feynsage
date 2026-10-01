# Kite, seeds up to 9 dots and 4 numerator powers, target F(3,3,2,3,3):
#   sage slides/bench_exact.sage ff | trimmed | exact
import sys, time; sys.path.insert(0, '.')
from feynsage import family, symmetries
from feynsage.laporta import Reducer
from feynsage.ff import reduce_ff, reduce_exact_trimmed
mode = sys.argv[1]
RMAX, SMAX = (int(sys.argv[2]), int(sys.argv[3])) if len(sys.argv) > 3 else (9, 4)
kite = family(["l1", "l1 + q", "l1 + l2", "l1 + l2 + q", "l2"], kin={"q^2": 1})
target = (3, 3, 2, 3, 3) if RMAX >= 9 else (2, 2, 1, 2, 2)
red = Reducer(kite, symmetries=symmetries(kite))
t0 = time.time()
if mode == "ff":
    res = reduce_ff(red, [target], rmax=RMAX, smax=SMAX, verbose=True)[target]
elif mode == "trimmed":
    res = reduce_exact_trimmed(red, [target], rmax=RMAX, smax=SMAX, verbose=True)[target]
else:
    red.run(rmax=RMAX, smax=SMAX, verbose=True)
    res = red.reduce(target)
T = time.time() - t0
print("MODE %s  time %.1f s" % (mode, T))
open('slides/out/bench_%s_%d_%d.txt' % (mode, RMAX, SMAX), 'w').write("%.1f\n%s\n" % (T, {str(k): str(v) for k, v in res.items()}))
