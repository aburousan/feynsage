# The four examples of the UF notebook shared by the speaker (uf-new-short.nb, a Mathematica routine
# that completes the square loop by loop).  feynsage must give the same U and F, both from the
# propagators (matrix method) and from the graph (spanning trees and 2-forests).
import sys; sys.path.insert(0, '.')
from feynsage import *

def check(name, fam, graph, U_ref, F_ref):
    U, F = fam.UF()
    R = U.parent()
    ok_fam = (U == R(U_ref)) and (F == R(F_ref))
    ok_graph = None
    if graph is not None:
        Rg = graph.R
        ok_graph = (graph.U() == Rg(U_ref)) and (graph.F() == Rg(F_ref))
    print("%-8s matrix method = notebook: %s   graph rules = notebook: %s" % (name, ok_fam, ok_graph))
    return ok_fam and ok_graph is not False

ok = True
# triangle, massless, Minkowski: p1^2 = P1Sq, p2^2 = P2Sq, (p1+p2)^2 = QSq
kin = Kinematics(['p1', 'p2'], ['P1Sq', 'P2Sq', 'QSq'],
                 {('p1','p1'): 'P1Sq', ('p2','p2'): 'P2Sq', ('p1','p2'): '(QSq - P1Sq - P2Sq)/2'}, euclidean=False)
fam = IntegralFamily('tri', ['k'], kin, [(mom(k=1), 0), (mom(k=1, p1=1), 0), (mom(k=1, p1=1, p2=1), 0)])
# line i is propagator i: 1 = k (Z-X), 2 = k + p1 (X-Y), 3 = k + p1 + p2 (Y-Z); p1 enters at X, p2 at Y
g = FeynmanGraph([('Z','X',0), ('X','Y',0), ('Y','Z',0)],
                 {'X': mom(p1=1), 'Y': mom(p2=1), 'Z': mom(p1=-1, p2=-1)}, kin)
ok &= check("triangle", fam, g, "x1 + x2 + x3", "-P1Sq*x1*x2 - QSq*x1*x3 - P2Sq*x2*x3")

# bubble with two masses
kin = Kinematics(['p1'], ['p1sq', 'm1', 'm2'], {('p1','p1'): 'p1sq'}, euclidean=False)
fam = IntegralFamily('bub', ['k'], kin, [(mom(k=1), 'm1^2'), (mom(k=1, p1=1), 'm2^2')])
g = FeynmanGraph([('A','B','m1^2'), ('A','B','m2^2')], {'A': mom(p1=1), 'B': mom(p1=-1)}, kin)
ok &= check("bubble", fam, g, "x1 + x2", "m1^2*x1^2 + m1^2*x1*x2 + m2^2*x1*x2 - p1sq*x1*x2 + m2^2*x2^2")

# massless box: p_i^2 = 0, p1.p2 = s/2, p1.p3 = t/2, p2.p3 = -(s+t)/2
kin = Kinematics(['p1', 'p2', 'p3'], ['s', 't'],
                 {('p1','p1'): 0, ('p2','p2'): 0, ('p3','p3'): 0, ('p1','p2'): 's/2', ('p1','p3'): 't/2', ('p2','p3'): '-(s+t)/2'},
                 euclidean=False)
fam = IntegralFamily('box', ['k'], kin, [(mom(k=1), 0), (mom(k=1, p1=1), 0), (mom(k=1, p1=1, p2=1), 0), (mom(k=1, p3=-1), 0)])
# lines 1 = k, 2 = k + p1, 3 = k + p1 + p2, 4 = k - p3 around the square
g = FeynmanGraph([('W','X',0), ('X','Y',0), ('Y','Z',0), ('Z','W',0)],
                 {'X': mom(p1=1), 'Y': mom(p2=1), 'Z': mom(p1=-1, p2=-1, p3=-1), 'W': mom(p3=1)}, kin)
ok &= check("box", fam, g, "x1 + x2 + x3 + x4", "-s*x1*x3 - t*x2*x4")

# two-loop sunset, three equal masses, p^2 = x m^2
kin = Kinematics(['p'], ['m', 'x'], {('p','p'): 'x*m^2'}, euclidean=False)
fam = IntegralFamily('sun', ['l', 'k'], kin, [(mom(k=1), 'm^2'), (mom(k=1, l=-1), 'm^2'), (mom(l=1, p=1), 'm^2')])
g = FeynmanGraph([('A','B','m^2'), ('A','B','m^2'), ('A','B','m^2')], {'A': mom(p=1), 'B': mom(p=-1)}, kin)
ok &= check("sunset", fam, g, "x2*x3 + x1*(x2 + x3)",
            "m^2*(x1^2*(x2 + x3) + x2*x3*(x2 + x3) + x1*(x2^2 - (-3 + x)*x2*x3 + x3^2))")
print("ALL AGREE WITH THE LECTURE NOTEBOOK" if ok else "MISMATCH")
