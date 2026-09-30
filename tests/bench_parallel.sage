# Finite-field reduction of the two-loop kite with several worker processes.
import sys, time; sys.path.insert(0, '.')
from feynsage import *
kite = family(["l1", "l1 + q", "l1 + l2", "l1 + l2 + q", "l2"], kin={"q^2": 1})
targets = [(3,3,2,3,3), (4,2,2,2,2), (2,2,2,2,-3)]
ref = None
for nproc in (1, 4, 8):
    t0 = time.time()
    res = ibp_reduce(kite, targets, method="ff", nproc=nproc)
    dt = time.time() - t0
    same = ref is None or all(res[t] == ref[t] for t in targets)
    ref = ref or res
    print("nproc = %d   %.1f s   same result: %s" % (nproc, dt, same), flush=True)
