#!/usr/bin/env python3
"""frontier_table.py -- the top level of the runs nextprime(10^c), c = 330..410, on one c8a.48xlarge
with yafu (logs in yafu_c8a_5/ and yafu_c8a_8/): search and build of the first level, the wins the
model expected at the first win, the whole chain.  Writes frontier.json for make_report.py."""
import json
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = ([(c, "yafu_c8a_5/c%d.log" % c) for c in range(330, 401, 10)] + [(410, "yafu_c8a_8/c410.log")]
        + [(c, "yafu_big/c%d.log" % c) for c in range(420, 451, 10)])
RESUMED = re.compile(r"search resumed from the checkpoint of \S+ \S+ on (\S+): \d+ orders bought, \d+ live candidates, "
                     r"B = \S+, cap = \S+, scan (\d+) \+ ecm (\d+) core-s, (\d+) s of search")
START = re.compile(r"\[(\d\d):(\d\d):(\d\d)\] proving p0")
STAMP = re.compile(r"\[(\d\d):(\d\d):(\d\d)\]")


def legs_account(txt):
    """Wall seconds and CPU seconds of a run that was written by several processes (a box taken
    away, a restart from a checkpoint).  A leg that ended in its own RUN line counts with what it
    measured.  A leg that was cut off counts with what the next leg took over from its last
    checkpoint: the search seconds and the scan and ECM core-seconds credited there, less those
    it had started from.  What a cut-off leg did after its last checkpoint is lost and not counted.
    Returns (wall, cpu, legs, cut-off legs)."""
    parts = re.split(r"(?=\[\d\d:\d\d:\d\d\] revised format)", txt)
    parts = [p for p in parts if "proving p0" in p]
    wall = cpu = 0.0
    cut = 0
    start = (0.0, 0.0)                      # (search seconds, core-seconds) the leg started from
    for i, p in enumerate(parts):
        r = RESUMED.search(p)
        if r:
            start = (float(r.group(4)), float(r.group(2)) + float(r.group(3)))
        m = RUN.search(p)
        if m:
            wall += float(m.group(1))
            cpu += float(m.group(2)) + float(m.group(3))
        else:
            cut += 1
            nxt = RESUMED.search(parts[i + 1]) if i + 1 < len(parts) else None
            if nxt:
                wall += float(nxt.group(4)) - start[0]
                cpu += float(nxt.group(2)) + float(nxt.group(3)) - start[1]
    return wall, cpu, len(parts), cut

SEARCH = re.compile(r"search: (\d+) orders, scan (\d+) \+ ecm (\d+) core-s \((\d+) above the rule\), ([\d.]+)s wall, .*?"
                    r"rule R = (\d+) core-s; the model expected ([\d.]+) wins.*?won at B1=(\d+)")
BUILT = re.compile(r"level (\d+) bits: built D=-(\d+) h~\d+ in ([\d.]+)s")
RUN = re.compile(r"RUN c\d+ wall=([\d.]+) user=([\d.]+) sys=([\d.]+)")
SCAN = re.compile(r"dscan\[[\d.e+]+,([\d.e+]+)\) (?:maxfb=|cap [\d.e+]+ -> )([\d.e+]+)")


def main():
    rows = []
    for c, f in RUNS:
        if not os.path.exists(os.path.join(HERE, f)):
            continue
        txt = open(os.path.join(HERE, f)).read()
        bits = int(re.search(r"proving p0 \((\d+) bits\)", txt).group(1))
        top = txt[:txt.index("level %d bits: built" % bits) + 200]
        s = SEARCH.search(top)
        b = [m for m in BUILT.finditer(txt) if int(m.group(1)) == bits][0]
        wall, cpu, nlegs, ncut = legs_account(txt)
        B = cap = 0.0
        for m in SCAN.finditer(top):
            B, cap = max(B, float(m.group(1))), max(cap, float(m.group(2)))
        rows.append({"c": c, "bits": bits, "wall": wall, "cpu": cpu, "legs": nlegs, "legs_cut": ncut,
                     "orders": int(s.group(1)), "scan_s": int(s.group(2)), "ecm_s": int(s.group(3)),
                     "search_wall": float(s.group(5)), "R": int(s.group(6)), "expwin": float(s.group(7)),
                     "win_b1": int(s.group(8)), "D": -int(b.group(2)), "build_s": float(b.group(3)),
                     "B": B, "cap": cap, "levels": len(BUILT.findall(txt)),
                     "verified": "verified twice" in txt})
    with open(os.path.join(HERE, "frontier.json"), "w") as f:
        json.dump(rows, f, indent=1)
    for r in rows:
        print("%d %4d bits %d legs (%d cut) chain %6.0f s %6.0f core-h | search %6.0f s, %7d orders, scan %5.0f + ecm %5.0f core-h, "
              "R %.1e, exp %.2f, won at %.0e | build %6.0f s D=%d | B %.1e cap %.1e %s"
              % (r["c"], r["bits"], r["legs"], r["legs_cut"], r["wall"], r["cpu"] / 3600, r["search_wall"], r["orders"], r["scan_s"] / 3600,
                 r["ecm_s"] / 3600, r["R"], r["expwin"], r["win_b1"], r["build_s"], r["D"], r["B"], r["cap"],
                 "ok" if r["verified"] else "NOT VERIFIED"))
    ew = [r["expwin"] for r in rows]
    print("expected wins at the first win: mean %.2f over %d primes" % (sum(ew) / len(ew), len(ew)))


if __name__ == "__main__":
    main()
