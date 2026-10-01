# Runs every computation shown on the slides and saves the outputs (slides/out/*.txt) and figures.
# Run from the repository root:  sage slides/compute.sage
import sys, time, json; sys.path.insert(0, '.')
from feynsage import *
from feynsage import plotting
import matplotlib; matplotlib.use('Agg')
O = 'slides/out/'
res = {}
def save(name, text):
    open(O + name + '.txt', 'w').write(str(text))
    print('==', name, '\n', text)

# 1. graph polynomials of the two-loop kite
g = diagram("kite")
save('kite_U', g.U()); save('kite_F', g.F())
fig = g.plot(); fig.savefig(O + 'kite.svg', bbox_inches='tight')

# 2. IBP reduction of the kite (finite fields)
kite = family(["l1", "l1 + q", "l1 + l2", "l1 + l2 + q", "l2"], kin={"q^2": 1})
t = time.time(); red = ibp_reduce(kite, ["F(2,2,1,2,2)"]); res['ibp_time'] = time.time() - t
save('kite_red', red)
save('ibp_time', "%.1f" % res['ibp_time'])

# 3. one-loop tensor reduction, symbols
t = time.time()
r = loop("l^mu l^nu", ["l", "m"], ["l + p", "m"], kin={"p^2": "s"})
save('bubble_tensor', r); save('bubble_time', "%.2f" % (time.time() - t))

# 4. closed forms
c = C0(0, 0, 5, 1, 1, 1)
save('c0_num', c.n(digits=25))
cc = explicit(c); save('c0_closed', cc)
d = D0(1, 2, 3, 4, -5, -6, 1, 1, 1, 1)
save('d0_num', d.n(digits=20).real())
save('d0_len', len(str(explicit(d))))

# 5. IR poles
save('ir_tri', C0(0, 5, 0, 0, 0, 0))
save('ir_box', D0(0, 0, 0, 0, -2, -3, 0, 0, 0, 0))
v = C0(1, -3, 1, 0, 1, 1)
save('ir_vertex', "1/eps: %s   finite: %s" % (pole_parts(v)[1].n(digits=12), finite_part(v, 1).n(digits=12)))

# 6. vanishing Gram determinant: g-2 at q^2 = 0
m = var('m')
save('gram_c1', PVC(0, 1, 0, m^2, 0, m^2, 0, m, m))
save('gram_c00', PVC(1, 0, 0, m^2, 0, m^2, 0, m, m).log_expand().expand())

# 7. a plot through thresholds
s = var('s')
fig = plotting.quick_plot([B0(s, 1, 1), C0(0, 0, s, 1, 1, 1)], (s, -4, 12), labels=[r"$B_0(s;1,1)$", r"$C_0(0,0,s;1,1,1)$"],
                          points=120, parts="both")
fig.savefig(O + 'threshold.svg', bbox_inches='tight')

# 8. the disputed C0
save('bug_c0', C0(-4, -3, -1, 0, 0, 1).n(digits=18))
json.dump(res, open(O + 'res.json', 'w'))

# ---- LaTeX versions of the outputs, typeset on the slides with mitex
def clean(t):
    t = t.replace(r"\mathit{eps}", r"\epsilon").replace(r"\, ", " ")
    t = t.replace(r"+ \left(-s\right)", "- s ").replace(r"\left(-s\right)", "-s ")
    return t.replace("+ -", "- ")
def tex(name, e):
    open(O + name + '.tex', 'w').write(clean(latex(e)))
tex('kite_U', g.U()); tex('kite_F', g.F())
row = list(red.table.values())[0]
def frac(c):
    return r"\frac{%s}{%s}" % (latex(c.numerator().factor()), latex(c.denominator().factor()))
open(O + 'kite_red.tex', 'w').write(clean(r"F(2,2,1,2,2) = " + r" \\ + ".join(
    r"%s\, F%s" % (frac(c), str(tuple(mm)).replace(' ', '')) for mm, c in row.items())))
for st, e in r.parts:
    tex('bubble_' + ('g' if 'g' in st else 'pp'), e)
def split_sum(e, per=3):
    ops = SR(e).operands() if SR(e).operator() is not None and 'add' in str(SR(e).operator()) else [e]
    rows = [ops[i:i + per] for i in range(0, len(ops), per)]
    body = r" \\ & + ".join(" + ".join(clean(latex(t)) for t in row) for row in rows)
    return r"\begin{aligned} & " + body.replace("+ -", "- ") + r"\end{aligned}"
open(O + 'c0_closed.tex', 'w').write(split_sum(cc))
def laurent(e):
    e = SR(e).subs(mu=1).expand()
    lead = latex(e.coefficient(eps, -2).simplify_full() / eps ** 2)
    c = [latex(e.coefficient(eps, j).simplify_full()) for j in (-1, 0)]
    return clean(r"%s + \frac{1}{\epsilon} \left(%s\right) + %s" % (lead, c[0], c[1]))
open(O + 'ir_tri.tex', 'w').write(laurent(C0(0, 5, 0, 0, 0, 0)))
open(O + 'ir_box.tex', 'w').write(laurent(D0(0, 0, 0, 0, -2, -3, 0, 0, 0, 0)))
tex('gram_c1', PVC(0, 1, 0, m^2, 0, m^2, 0, m, m))
tex('gram_c00', PVC(1, 0, 0, m^2, 0, m^2, 0, m, m).log_expand().expand())
