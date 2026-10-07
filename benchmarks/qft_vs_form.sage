# Overhead of the Pythonic layer (feynsage.qft) against FORM run by hand on the same program.
# Tr(p1/ ... pn/) in 4 dimensions (trace4) and in d dimensions (tracen), n distinct momenta.
#   sage benchmarks/qft_vs_form.sage [nmax]
import sys, time, subprocess, os, tempfile
sys.path.insert(0, '.')
from feynsage import *
from feynsage import form as _form
from feynsage.qft.parser import Reader, expression_text
from feynsage.qft.compiler import Program

nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 12
d = var('d')
rows = []


def form_by_hand(n, dim4):
    decl = ("" if dim4 else "Symbol d;\nDimension d;\n") + "Vectors q1,...,q%d;\n" % n
    code = decl + "Off statistics;\nLocal F = g_(1,q1,...,q%d);\n%s,1;\nFormat nospaces;\nPrint;\n.end\n" % (n, "trace4" if dim4 else "tracen")
    with tempfile.TemporaryDirectory() as tmp:
        f = os.path.join(tmp, 'b.frm')
        open(f, 'w').write(code)
        t = time.time()
        out = subprocess.run([_form.FORM, '-q', f], capture_output=True, text=True, cwd=tmp).stdout
        return time.time() - t, len(out)


for n in range(4, nmax + 1, 2):
    qs = momenta(" ".join("q%d" % i for i in range(1, n + 1)))
    qs = qs if isinstance(qs, tuple) else (qs,)
    for dim in (4, d):
        tf, size = form_by_hand(n, dim == 4)
        expr = prod(slash(x) for x in qs)
        set_backend("subprocess")
        t = time.time(); r1 = dirac_trace(expr, dim=dim); t1 = time.time() - t
        set_backend("persistent")
        dirac_trace(slash(qs[0]) * slash(qs[0]))                         # warm the worker
        t = time.time(); r2 = dirac_trace(expr, dim=dim); t2 = time.time() - t
        same = (r1 - r2).is_trivial_zero() if n < 12 else bool(r1.nops() == r2.nops())
        rows.append((n, "4" if dim == 4 else "d", r1.nops() if n > 2 else 1, size, tf, t1, t2, same))
        print("n=%2d dim=%s terms=%7d  FORM by hand %7.3f s | dirac_trace: new process %7.3f s, persistent %7.3f s  %s"
              % (rows[-1][:3] + rows[-1][4:7] + ("" if same else "MISMATCH",)), flush=True)
print("\n| n | dim | terms | FORM by hand (s) | dirac_trace, process per call (s) | dirac_trace, persistent (s) |")
print("|---|---|---|---|---|---|")
for n, dm, nt, size, tf, t1, t2, same in rows:
    print("| %d | %s | %d | %.3f | %.3f | %.3f |" % (n, dm, nt, tf, t1, t2))
