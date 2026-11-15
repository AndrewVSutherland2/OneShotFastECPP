#!/usr/bin/env python3
"""Summarise the prover runs with the stopping-rule scheduler (logs in rule768_local/) into rule_runs.json.

One record per run: the top level's search line (orders bought, scan and ECM core-seconds, wall, the reward R of
the rule at the end, the wins the model expected from the work done), the build, and the whole chain.
"Core-seconds" on the 16-core driver box are thread-seconds of its 32 hardware threads."""
import glob, json, os, re
HERE = os.path.dirname(os.path.abspath(__file__))
rows = []
for fn in sorted(glob.glob(os.path.join(HERE, "rule768_local", "p*.log")) + glob.glob(os.path.join(HERE, "rule1030_local", "*.log")),
                 key=lambda f: (os.path.dirname(f), int(re.search(r"p(\d+)", os.path.basename(f)).group(1)) if re.search(r"p(\d+)", os.path.basename(f)) else 0)):
    txt = open(fn).read()
    bits = int(re.search(r"proving p0 \((\d+) bits\)", txt).group(1))
    m = re.search(r"search: (\d+) orders, scan (\d+) \+ ecm (\d+) core-s(?: \((\d+) above the rule\))?, ([\d.]+)s wall, (\d+) units, (\d+) splits, "
                  r"rule R = (\d+) core-s(?:; the model expected (-?[\d.]+) wins from this work)?", txt)
    if not m:
        continue
    b = re.search(r"level %d bits: built D=-(\d+) h~(\d+) in ([\d.]+)s" % bits, txt)
    c = re.search(r"chain complete: (\d+) levels, ([\d.]+)s", txt)
    f = re.search(r"cost factor ([\d.]+)", txt)
    fd = re.search(r"dscan\[0,[^)]*\)[^:]*: \+(\d+) D", txt)
    hd = re.search(r"H_D", txt)
    rows.append({"file": os.path.relpath(fn, HERE), "bits": bits, "orders": int(m.group(1)), "scan_s": int(m.group(2)),
                 "ecm_s": int(m.group(3)), "over_s": int(m.group(4)) if m.group(4) else None, "search_wall": float(m.group(5)),
                 "units": int(m.group(6)), "splits": int(m.group(7)), "R": int(m.group(8)),
                 "expwin": float(m.group(9)) if m.group(9) else None,
                 "D": -int(b.group(1)) if b else None, "build_s": float(b.group(3)) if b else None,
                 "levels": int(c.group(1)) if c else None, "chain_s": float(c.group(2)) if c else None,
                 "cost_factor": float(f.group(1)) if f else None, "first_chunk_D": int(fd.group(1)) if fd else None,
                 "verified": "verified twice" in txt})
json.dump(rows, open(os.path.join(HERE, "rule_runs.json"), "w"), indent=1)
for r in rows:
    print("%-34s %4d bits  %6d orders  scan %5d + ecm %6d core-s  search %6.1fs  E[wins] %s  build %5.1fs  chain %6.1fs  %s"
          % (r["file"], r["bits"], r["orders"], r["scan_s"], r["ecm_s"], r["search_wall"],
             "%5.2f" % r["expwin"] if r["expwin"] is not None else "  n/a", r["build_s"] or 0, r["chain_s"] or 0,
             "verified" if r["verified"] else "NOT VERIFIED"))
a = [r for r in rows if r["bits"] < 800]
if a:
    n = len(a)
    w = sum(r["scan_s"] + r["ecm_s"] for r in a) / n
    e = [r["expwin"] for r in a if r["expwin"] is not None]
    print("%d runs below 800 bits: mean search work %.0f core-s (%.2f core-h), mean search wall %.0f s, mean chain %.0f s, median chain %.0f s"
          % (n, w, w / 3600, sum(r["search_wall"] for r in a) / n, sum(r["chain_s"] for r in a) / n, sorted(r["chain_s"] for r in a)[n // 2]))
    print("expected wins at the win: mean %.2f over %d runs (sum %.2f)" % (sum(e) / len(e), len(e), sum(e)))
