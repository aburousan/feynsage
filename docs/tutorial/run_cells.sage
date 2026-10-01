# Runs the code cells of the tutorial (cells.py) one after another in one namespace, like a
# Jupyter notebook, and stores what each cell prints, the value of its last line (as LaTeX
# or as text) and any figure in out/cells.json.  Run from the repository root:
#     sage docs/tutorial/run_cells.sage
import sys, os, io, ast, json, contextlib, time, re
sys.path.insert(0, '.')
from sage.repl.preparse import preparse
import matplotlib; matplotlib.use('Agg')
import matplotlib.figure
HERE = 'docs/tutorial/'
exec(open(HERE + 'cells.py').read(), globals())            # defines CELLS = [(id, code, opts), ...]
ns = {}
exec("from sage.all import *", ns)
out = {}
for n, (cid, code, opts) in enumerate(CELLS, 1):
    src = preparse(code)
    tree = ast.parse(src)
    last = tree.body[-1] if tree.body and isinstance(tree.body[-1], ast.Expr) and not code.rstrip().endswith(';') else None
    body = tree.body[:-1] if last is not None else tree.body
    buf = io.StringIO()
    t0 = time.time()
    with contextlib.redirect_stdout(buf):
        exec(compile(ast.Module(body=body, type_ignores=[]), cid, 'exec'), ns)
        value = eval(compile(ast.Expression(body=last.value), cid, 'eval'), ns) if last is not None else None
    rec = {"code": code.strip("\n"), "stdout": buf.getvalue().rstrip("\n"), "latex": None, "text": None, "fig": None,
           "seconds": float(round(float(time.time() - t0), 2)), "n": int(n)}
    fig = opts.get("fig")
    if fig is not None:
        f = ns[fig] if isinstance(fig, str) else value
        f.savefig(HERE + 'out/%s.svg' % cid, bbox_inches='tight')
        rec["fig"] = 'out/%s.svg' % cid
        value = None
    if value is not None:
        if isinstance(value, matplotlib.figure.Figure):
            value.savefig(HERE + 'out/%s.svg' % cid, bbox_inches='tight')
            rec["fig"] = 'out/%s.svg' % cid
        elif opts.get("text") or isinstance(value, (str, bool, int)):
            rec["text"] = str(value)
        elif hasattr(value, '_latex_') or hasattr(value, 'parent') or isinstance(value, (tuple, list, dict)):
            tex = str(latex(value)) if not hasattr(value, '_latex_') else value._latex_()
            # matrices as pmatrix (Sage writes \left(\begin{array}{rr..}...\end{array}\right))
            tex = re.sub(r'\\left\(\\begin\{array\}\{[lcr]+\}', r'\\begin{pmatrix}', tex)
            tex = tex.replace('\\end{array}\\right)', '\\end{pmatrix}')
            rec["latex"] = tex
        else:
            rec["text"] = repr(value)
    out[cid] = rec
    print("[%2d] %-22s %.2f s" % (n, cid, rec["seconds"]))
json.dump(out, open(HERE + 'out/cells.json', 'w'), indent=1)
