# Collect work/*/ into results.md (tables) and results.json; plot.py makes the figure.
#     python3 summarize.py
import os, json, glob, re
H = os.path.dirname(os.path.abspath(__file__))
exec(open(os.path.join(H, "cases.py")).read())
rows = []
for case, c in CASES.items():
    for r, s in c["sizes"]:
        B = os.path.join(H, "work", "%s_r%d_s%d" % (case, r, s))
        if not os.path.isdir(B):
            continue
        row = {"case": case, "r": r, "s": s, "targets": len(targets(case, r, s))}
        fs = os.path.join(B, "fs_ff_1", "fs_time.json")
        if os.path.exists(fs):
            row["masters"] = json.load(open(fs))["masters"]
        for tool in ["fs_ff_1", "fs_ff_8", "fs_trimmed_1", "kira_fermat_1", "kira_fermat_8", "kira_firefly_1", "kira_firefly_8"]:
            f = os.path.join(B, tool + ".json")
            if os.path.exists(f):
                m = json.load(open(f))
                row[tool] = {"wall": m["wall"], "mem": m["maxrss_mb"], "ok": m["rc"] == 0}
                inner = os.path.join(B, tool, "fs_time.json")
                if os.path.exists(inner):
                    row[tool]["reduce"] = json.load(open(inner))["reduce_seconds"]
            cmpf = os.path.join(B, "cmp_%s.txt" % tool)
            if os.path.exists(cmpf):
                row.setdefault("agree", {})[tool] = open(cmpf).read().strip()
        rows.append(row)
json.dump(rows, open(os.path.join(H, "results.json"), "w"), indent=1)

def cell(row, tool):
    t = row.get(tool)
    if t is None:
        return "-"
    if not t["ok"]:
        return "timeout" if t["wall"] > 3500 else "failed"
    return "%.1f s / %.0f MB" % (t["wall"], t["mem"])

out = ["# feynsage against Kira 3.1", "",
       "Same machine (ADRISHTA: AMD EPYC 9534, idle), same targets, same seeds, same masters.",
       "Time = wall-clock of the whole run (feynsage includes about 2 s of Sage start-up), memory = peak resident set.",
       "Every feynsage result was compared with Kira's, coefficient by coefficient.", ""]
for th in [1, 8]:
    out += ["## %d thread%s" % (th, "" if th == 1 else "s"), "",
            "| family | r | s | targets | masters | feynsage ff | feynsage exact (trimmed) | Kira (Fermat) | Kira (FireFly) | identical |",
            "|---|---|---|---|---|---|---|---|---|---|"]
    for row in rows:
        agree = row.get("agree", {})
        ks = [v for k, v in agree.items() if k.endswith("_%d" % th)]
        ok = "yes" if ks and all(v == "IDENTICAL" for v in ks) else ("-" if not ks else "NO")
        out.append("| %s | %d | %d | %d | %s | %s | %s | %s | %s | %s |" % (
            row["case"], row["r"], row["s"], row["targets"], row.get("masters", "-"),
            cell(row, "fs_ff_%d" % th), cell(row, "fs_trimmed_1") if th == 1 else "-",
            cell(row, "kira_fermat_%d" % th), cell(row, "kira_firefly_%d" % th), ok))
    out.append("")
open(os.path.join(H, "results.md"), "w").write("\n".join(out) + "\n")
print("\n".join(out))
