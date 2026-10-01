# Graph polynomials at two and three loops: spanning trees and 2-forests against the
# matrix-tree theorem (U) and the momentum-space matrix method (U and F).
import sys; sys.path.insert(0, '.')
from feynsage import *
from feynsage.easy import graph, diagram

def check(name, g, U_ref=None, F_ref=None):
    U, F = g.U(), g.F()
    props, loops = g.family()
    Um, Fm = IntegralFamily(name, loops, g.kin, props).UF()
    Um, Fm = g.R(str(Um)), g.R(str(Fm))
    ok = g.U_kirchhoff() == U and Um == U and Fm == F
    if U_ref is not None:
        ok = ok and U == g.R(U_ref) and F == g.R(F_ref)
    print("%-10s L=%d trees=%3d 2-forests=%3d  %s" % (name, g.L, len(g.spanning_trees()), len(g.two_forests()), "PASS" if ok else "FAIL"))
    assert ok

# Minkowski (default): F = -F_0 + U sum x m^2, so the vertex of the notes has F = -q2 [...]
check("vertex2", diagram("vertex2"), "(x2+x3)*(x1+x4+x5) + (x1+x2+x3+x4+x5)*x6",
      "-q2*(x1*x3*x4 + x1*x2*(x3+x4) + x2*x3*(x4+x5) + (x1+x2)*(x3+x4)*x6)")
check("banana3", diagram("banana3"), "x2*x3*x4 + x1*x3*x4 + x1*x2*x4 + x1*x2*x3",
      "-s*x1*x2*x3*x4 + (x2*x3*x4 + x1*x3*x4 + x1*x2*x4 + x1*x2*x3)*m^2*(x1+x2+x3+x4)")
check("ladder3", diagram("ladder3"))
check("triplebox", diagram("triplebox"))
check("sunbub3", graph("A-B:m, A-C:m, B-C:m, B-C:m, B-C:m", {"A": "p", "C": "-p"}, kin={"p^2": "s"}))
