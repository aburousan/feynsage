import sys; sys.path.insert(0, '.')
from feynsage import *

def check(name, g, U_ref, F_ref):
    U, F = g.U(), g.F()
    Uk = g.U_kirchhoff()
    props, loops = g.family()
    fam = IntegralFamily(name, loops, g.kin, props)
    Um = Fm = None
    if fam is not None:
        Um, Fm = fam.UF()
        Um = g.R(str(Um)); Fm = g.R(str(Fm))
    print(name, "| trees:", len(g.spanning_trees()), "| U tree = ref:", U == g.R(U_ref), "| F tree = ref:", F == g.R(F_ref),
          "| Kirchhoff U = tree U:", Uk == U, "| matrix method:", (Um == U and Fm == F) if fam else "n/a (needs ISPs)")

# two-loop kite, Euclidean, p^2 = pp.  Vertices L, T, B, R.  Lines 1 L-B, 2 L-T, 3 B-R, 4 T-R, 5 T-B
kin = Kinematics(['p'], ['pp'], {('p','p'): 'pp'}, euclidean=True)
g = FeynmanGraph([('L','B',0), ('L','T',0), ('B','R',0), ('T','R',0), ('T','B',0)], {'L': mom(p=1), 'R': mom(p=-1)}, kin)
check('kite', g, "(x1+x2)*(x3+x4) + x5*(x1+x2+x3+x4)",
      "pp*(x1*x2*(x3+x4) + x3*x4*(x1+x2) + x5*(x1+x3)*(x2+x4))")

# massless box, Euclidean: p1..p3 basis, p4 = -p1-p2-p3, all p_i^2 = 0, s = (p1+p2)^2, t = (p2+p3)^2.
# cyclic lines: 1 left (V1-V2), 2 top (V2-V3), 3 right (V3-V4), 4 bottom (V4-V1); p1 at V1, p2 at V2, p3 at V3, p4 at V4
kb = Kinematics(['p1','p2','p3'], ['s','t'], {('p1','p1'):0, ('p2','p2'):0, ('p3','p3'):0,
      ('p1','p2'):'s/2', ('p2','p3'):'t/2', ('p1','p3'):'-(s+t)/2'}, euclidean=True)
gb = FeynmanGraph([('V1','V2',0), ('V2','V3',0), ('V3','V4',0), ('V4','V1',0)],
                  {'V1': mom(p1=1), 'V2': mom(p2=1), 'V3': mom(p3=1), 'V4': mom(p1=-1, p2=-1, p3=-1)}, kb)
check('box', gb, "x1+x2+x3+x4", "s*x2*x4 + t*x1*x3")

# phi^4 two-loop diagram of the notes: A (p1, p2 in), B (p3 out), C (p4 out); lines 1 A-B, 2 A-C, 3 and 4 B-C.
# basis {p3, p4}: incoming at A is p3 + p4, at B it is -p3, at C it is -p4; p3^2 = p4^2 = 0, 2 p3.p4 = Q2
kf = Kinematics(['p3','p4'], ['Q2'], {('p3','p3'):0, ('p4','p4'):0, ('p3','p4'):'Q2/2'}, euclidean=True)
gf = FeynmanGraph([('A','B',0), ('A','C',0), ('B','C',0), ('B','C',0)],
                  {'A': mom(p3=1, p4=1), 'B': mom(p3=-1), 'C': mom(p4=-1)}, kf)
check('phi4-two', gf, "x3*x4 + (x1+x2)*(x3+x4)", "Q2*x1*x2*(x3+x4)")

# bubble with two masses
km = Kinematics(['p'], ['pp','ma','mb'], {('p','p'): 'pp'}, euclidean=True)
gm = FeynmanGraph([('A','B','ma'), ('A','B','mb')], {'A': mom(p=1), 'B': mom(p=-1)}, km)
check('bubble m1 m2', gm, "x1+x2", "pp*x1*x2 + (ma*x1+mb*x2)*(x1+x2)")
