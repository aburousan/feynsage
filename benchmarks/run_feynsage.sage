# feynsage side of one benchmark point:
#     sage run_feynsage.sage <case> <r> <s> <nproc> <method> <outdir>
# Writes <outdir>/feynsage.txt (the reduction), preferred (the masters, for Kira) and fs_time.json.
import sys, os, time, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), '..'))
from feynsage import *
from feynsage.laporta import Reducer
from feynsage.ff import reduce_ff, reduce_exact_trimmed
# (no Python variable may be called like a kinematic symbol: Sage would read it into the strings below)
case, rmax, smax, nproc, method, outdir = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5], sys.argv[6]
exec(open(os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), 'cases.py')).read())
c = CASES[case]
one = SR.var(c["one"])
props = [(q, "1" if m2 == c["one"] else m2) for q, m2 in c["props"]]       # mass^2 = 1 means mass 1
kin = {(a + "^2" if a == b else a + "." + b): str(SR(v).subs({one: 1})) for (a, b), v in c["sp"].items()}
fam = family(props, kin=kin, loops=c["loops"], name=case)
tg = targets(case, rmax, smax)
t0 = time.time()
red = Reducer(fam, symmetries=symmetries(fam))
if method == "ff":
    table = reduce_ff(red, tg, rmax=rmax, smax=smax, nproc=nproc)
else:
    table = reduce_exact_trimmed(red, tg, rmax=rmax, smax=smax)
t1 = time.time()
masters = sorted({m for row in table.values() for m in row})
os.makedirs(outdir, exist_ok=True)
with open(os.path.join(outdir, "feynsage.txt"), "w") as f:
    for t, row in table.items():
        f.write("%s|%s\n" % (",".join(map(str, t)), ";".join("%s:%s" % (",".join(map(str, m)), cc) for m, cc in row.items())))
open(os.path.join(outdir, "preferred"), "w").write("".join("%s[%s]\n" % (case, ",".join(map(str, m))) for m in masters))
open(os.path.join(outdir, "targets"), "w").write("".join("%s[%s]\n" % (case, ",".join(map(str, t))) for t in tg))
json.dump({"reduce_seconds": t1 - t0, "targets": len(tg), "masters": len(masters)},
          open(os.path.join(outdir, "fs_time.json"), "w"))
print("feynsage %s r=%d s=%d: %d targets, %d masters, %.2f s" % (case, rmax, smax, len(tg), len(masters), t1 - t0))
