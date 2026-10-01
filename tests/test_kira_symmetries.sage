# Families that need sector symmetries, against Kira 3.1 (reference results in tests/data, made on
# ADRISHTA with the masters below given to Kira as preferred masters).  Kira found the same number
# of masters by itself.  Integrals with numerators included.
#  - the planar two-loop vertex of the notes plus the numerator line (l2 + p1)^2, q^2 = 1
#  - the two-loop sunset with three equal masses (m^2 = 1) and the numerator lines k + p, l
import sys, re; sys.path.insert(0, '.')
from feynsage import *

def kira_table(path, name, R, names):
    txt = open(path).read().strip().strip('{}')
    out = {}
    for block in txt.split('\n,\n'):
        lhs, rhs = block.split('->', 1)
        t = tuple(int(x) for x in lhs.split('[')[1].split(']')[0].split(','))
        row = {}
        for ln in rhs.strip().split('\n'):
            mm = re.match(r'\s*\+\s*' + name + r'\[([-\d,]+)\]\*\((.*)\)\s*$', ln)
            row[tuple(int(x) for x in mm.group(1).split(','))] = R(sage_eval(mm.group(2), locals=names))
        out[t] = row
    return out

def check(label, fam, targets, path, name, nmasters):
    res = ibp_reduce(fam, targets)
    K = fam.kin.K
    names = {str(g): K(g) for g in fam.kin.R.gens()}
    kira = kira_table(path, name, K, names)
    same = len(res.masters) == nmasters and all(
        set(row) == set(kira[t]) and all(K(row[m]) == kira[t][m] for m in row) for t, row in res.table.items())
    print("%-34s %d masters, every coefficient identical to Kira: %s" % (label, len(res.masters), same))
    assert same

check("two-loop vertex + numerator line",
      family(["l1", "l2", "-l2 + p1 + p2", "-l1 + p1 + p2", "l1 - p1", "-l1 + l2", "l2 + p1"],
             kin={"p1^2": 0, "p2^2": 0, "p1.p2": "1/2"}),
      ["F(1,1,1,1,1,1,-1)", "F(2,1,1,1,1,1,0)", "F(1,1,1,1,1,1,-2)", "F(1,2,1,1,1,1,-1)", "F(1,1,1,1,1,2,0)"],
      'tests/data/kira_vertex2.m', 'v2', 3)
check("sunset, three equal masses",
      family([("k", "1"), ("k - l", "1"), ("l + p", "1"), "k + p", "l"], kin={"p^2": "pp"}),
      ["F(2,1,1,0,0)", "F(2,2,1,0,0)", "F(1,1,1,-1,0)", "F(1,1,1,-1,-1)", "F(2,1,1,-1,0)", "F(1,1,3,0,0)"],
      'tests/data/kira_sunset.m', 'sun', 3)
print("ALL AGREE WITH KIRA")
