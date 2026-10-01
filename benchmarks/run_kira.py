# Kira side of one benchmark point:
#     python3 run_kira.py <case> <r> <s> <threads> <mode: fermat|firefly> <outdir> <kira binary>
# Uses <outdir>/targets and <outdir>/preferred written by run_feynsage.sage, so that both programs
# reduce the same integrals to the same masters.
import sys, os, shutil, subprocess, json
case, r, s, threads, mode, outdir, kira = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5], sys.argv[6], sys.argv[7]
here = os.path.dirname(os.path.abspath(sys.argv[0]))
exec(open(os.path.join(here, "cases.py")).read())
c = CASES[case]
work = os.path.join(outdir, "kira_%s_%d" % (mode, threads))
shutil.rmtree(work, ignore_errors=True)
os.makedirs(os.path.join(work, "config"))
top_bits = sum(1 << i for i in range(c["top"]))
fam = ["integralfamilies:", '  - name: "%s"' % case, "    loop_momenta: [%s]" % ", ".join(c["loops"]),
       "    top_level_sectors: [%d]" % top_bits, "    propagators:"]
fam += ['      - [ "%s", %s ]' % (q, m2) for q, m2 in c["props"]]
open(os.path.join(work, "config", "integralfamilies.yaml"), "w").write("\n".join(fam) + "\n")
kin = ["kinematics:", "  incoming_momenta: [%s]" % ", ".join(c["ext"]), "  kinematic_invariants:"]
kin += ["    - [%s, %d]" % (n, dim) for n, dim in c["invariants"]]
kin += ["  scalarproduct_rules:"] + ["    - [[%s,%s], %s]" % (a, b, v) for (a, b), v in c["sp"].items()]
kin += ["  symbol_to_replace_by_one: %s" % c["one"]]
open(os.path.join(work, "config", "kinematics.yaml"), "w").write("\n".join(kin) + "\n")
src = os.path.join(outdir, "fs_ff_1")                   # targets and masters from feynsage's run
if not os.path.exists(os.path.join(src, "targets")):
    os.makedirs(src, exist_ok=True)
    open(os.path.join(src, "targets"), "w").write("".join("%s[%s]\n" % (case, ",".join(map(str, t))) for t in targets(case, r, s)))
shutil.copy(os.path.join(src, "targets"), os.path.join(work, "targets"))
have_pref = os.path.exists(os.path.join(src, "preferred"))
if have_pref:
    shutil.copy(os.path.join(src, "preferred"), os.path.join(work, "preferred"))
run = ["      run_initiate: true"]
run += ["      run_firefly: true"] if mode == "firefly" else ["      run_triangular: true", "      run_back_substitution: true"]
jobs = ["jobs:", "  - reduce_sectors:", "      reduce:",
        "        - {topologies: [%s], sectors: [%d], r: %d, s: %d}" % (case, top_bits, c["top"] + r, s),
        "      select_integrals:", "        select_mandatory_list:", "          - [%s, targets]" % case,
        ] + (["      preferred_masters: preferred"] if have_pref else []) + run + ["  - kira2math:", "      target:", "        - [%s, targets]" % case]
open(os.path.join(work, "jobs.yaml"), "w").write("\n".join(jobs) + "\n")
log = open(os.path.join(work, "kira.log"), "w")
rc = subprocess.call([kira, "--parallel=%d" % threads, "jobs.yaml"], cwd=work, stdout=log, stderr=subprocess.STDOUT)
print("kira %s %s r=%d s=%d threads=%d: rc %s" % (mode, case, r, s, threads, rc))
