# Compare feynsage's reduction with Kira's, coefficient by coefficient.
#     sage compare.sage <case> <outdir> <kira workdir>     (prints IDENTICAL or the differences)
import sys, os, re
case, outdir, work = sys.argv[1], sys.argv[2], sys.argv[3]
exec(open(os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), 'cases.py')).read())
c = CASES[case]
names = ['d'] + [n for n, _ in c["invariants"] if n != c["one"]]
R = PolynomialRing(QQ, names).fraction_field()
loc = {n: R.gen(i) for i, n in enumerate(names)}
loc[c["one"]] = R(1)
ours = {}
for line in open(os.path.join(outdir, "feynsage.txt")):
    t, rest = line.strip().split("|")
    row = {}
    if rest:
        for item in rest.split(";"):
            m, cc = item.split(":", 1)
            row[tuple(int(x) for x in m.split(","))] = R(sage_eval(cc, locals=loc))
    ours[tuple(int(x) for x in t.split(","))] = row
path = os.path.join(work, "results", case, "kira_targets.m")
if not os.path.exists(path):
    print("NO KIRA RESULT"); sys.exit(0)
txt = open(path).read().strip().strip("{}")
kira = {}
for block in txt.split("\n,\n"):
    lhs, rhs = block.split("->", 1)
    t = tuple(int(x) for x in lhs.split("[")[1].split("]")[0].split(","))
    row = {}
    for ln in rhs.strip().split("\n"):
        ln = ln.strip()
        if not ln or ln == "0":
            continue
        mm = re.match(r"\+?\s*" + case + r"\[([-\d,]+)\]\*\((.*)\)\s*$", ln)
        row[tuple(int(x) for x in mm.group(1).split(","))] = R(sage_eval(mm.group(2), locals=loc))
    kira[t] = row
for t in ours:                       # Kira does not list a target that is itself a master
    if t not in kira and ours[t] == {t: R(1)}:
        kira[t] = {t: R(1)}
bad = [t for t in ours if t not in kira or set(ours[t]) != set(kira[t]) or any(ours[t][m] != kira[t][m] for m in ours[t])]
print("IDENTICAL" if not bad else "DIFFERENT: %d of %d targets, e.g. %s" % (len(bad), len(ours), bad[:3]))
