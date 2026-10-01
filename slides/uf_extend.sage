# The speaker's notebook (uf-new-short.nb) finds U and F by completing the square, loop by loop,
# for four diagrams.  Here the same polynomials are found for bigger diagrams in three independent
# ways inside feynsage (graph rules: spanning trees and 2-forests; Kirchhoff matrix-tree theorem;
# matrix method U = det M) and written out for FeynCalc's FCFeynmanPrepare as an outside check.
import sys, time, json; sys.path.insert(0, '.')
from feynsage import *
from feynsage.family import IntegralFamily
O = 'slides/out/'
diagrams = {
    "kite": (graph("L-B, L-T, B-R, T-R, T-B", {"L": "p", "R": "-p"}, kin={"p^2": "s"}),
             {"p.p": "s"}),
    "double box": (graph("A-B, B-E, E-C, C-D, D-F, F-A, E-F",
                         {"A": "p1", "B": "p2", "C": "p3", "D": "-p1 - p2 - p3"},
                         kin={"p1^2": 0, "p2^2": 0, "p3^2": 0, "p1.p2": "s/2", "p2.p3": "t/2", "p1.p3": "-(s+t)/2"}),
                   {"p1.p1": "0", "p2.p2": "0", "p3.p3": "0", "p1.p2": "s/2", "p2.p3": "t/2", "p1.p3": "-(s+t)/2"}),
    "banana (3 loops)": (graph("A-B:m, A-B:m, A-B:m, A-B:m", {"A": "p", "B": "-p"}, kin={"p^2": "s"}),
                         {"p.p": "s"}),
    "triple box (3 loops)": (graph("T1-T2, T2-T3, T3-T4, B1-B2, B2-B3, B3-B4, T1-B1, T2-B2, T3-B3, T4-B4",
                         {"B1": "p1", "T1": "p2", "T4": "p3", "B4": "-p1 - p2 - p3"},
                         kin={"p1^2": 0, "p2^2": 0, "p3^2": 0, "p1.p2": "s/2", "p2.p3": "t/2", "p1.p3": "-(s+t)/2"}),
                   {"p1.p1": "0", "p2.p2": "0", "p3.p3": "0", "p1.p2": "s/2", "p2.p3": "t/2", "p1.p3": "-(s+t)/2"}),
}
rows, wl = [], ['$LoadAddOns = {}; Get["FeynCalc`"];', 'out = OpenWrite["%s/slides/out/feyncalc_uf.txt"];' % __import__('os').getcwd()]
def fmom(d):
    s = " + ".join("(%s)*%s" % (c, k) for k, c in d.items())
    return s if s else "0"
for name, (g, sp) in diagrams.items():
    t = time.time(); U, F = g.U(), g.F(); tg = time.time() - t
    t = time.time(); Uk = g.U_kirchhoff(); tk = time.time() - t
    props, loops = g.family()
    fam = IntegralFamily('f', loops, g.kin, props)
    t = time.time(); Um, Fm = fam.UF(); tm = time.time() - t
    same_k = (Uk == U)
    same_m = (g.R(str(Um)) == U) and (g.R(str(Fm)) == F)
    rows.append({"name": name, "L": len(loops), "N": len(props), "Uterms": len(U.monomials()), "Fterms": len(F.monomials()),
                 "trees": len(g.spanning_trees()), "kirchhoff": bool(same_k), "matrix": bool(same_m),
                 "t_graph": tg, "t_matrix": tm, "U": str(U), "F": str(F)})
    print(name, rows[-1]["Uterms"], rows[-1]["Fterms"], same_k, same_m, "%.2f %.2f" % (tg, tm))
    # FeynCalc input: same propagators, same order
    fads = ", ".join("FAD[{%s, %s}]" % (fmom(q), ("Sqrt[%s]" % m2) if str(m2) != "0" else 0) for q, m2 in props)
    sps = "; ".join("SPD[%s, %s] = %s" % (k.split('.')[0], k.split('.')[1], v) for k, v in sp.items())
    wl.append('ClearScalarProducts[]; %s;' % sps)
    wl.append('r = FCFeynmanPrepare[{%s}, {%s}, Names -> x, Indexed -> False];' % (fads, ", ".join(loops)))
    wl.append('WriteString[out, "%s|" <> ToString[InputForm[Expand[r[[1]]]]] <> "|" <> ToString[InputForm[Expand[r[[2]]]]] <> "\\n"];' % name)
wl.append('Close[out];')
open(O + 'feyncalc_uf.wl', 'w').write("\n".join(wl))
json.dump(rows, open(O + 'uf_rows.json', 'w'), indent=1)

# compare with FeynCalc (run slides/out/feyncalc_uf.wl in Mathematica first; this part is skipped if
# its output is missing)
import os
if os.path.exists(O + 'feyncalc_uf.txt'):
    ours = {r["name"]: r for r in rows}
    same = {}
    for line in open(O + 'feyncalc_uf.txt'):
        name, Ufc, Ffc = line.strip().split('|')
        r = ours[name]
        same[name] = bool((SR(Ufc) - SR(r["U"])).expand() == 0 and (SR(Ffc) - SR(r["F"])).expand() == 0)
        print(name, "FeynCalc identical:", same[name])
    json.dump(same, open(O + 'uf_feyncalc.json', 'w'))
