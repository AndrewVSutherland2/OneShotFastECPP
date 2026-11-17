#!/usr/bin/env python3
"""The entries of the ShortPrimalityProofs table, nextprime(10^c) for c = 10, 20, ..., 320, as proved by
shortECPP.py on one c8a.48xlarge (192 cores, yafu as the ECM engine, revised format), against the times
recorded in that repository's README (upstream main at 12aa5d7).

Three passes were run on 2026-11-16 (aws-spot-7):
  yafu_c8a_1/  the code of commit b445862: c = 160..320 and the 768- and 1030-bit primes of section 10
  yafu_c8a_2/  yafu sweeps completed, small levels sized to their work, rule from 280 bits (b0f14e4);
               the full set of modular polynomials was on the box from c = 250 on; also c = 10..150
  yafu_c8a_3/  rule table computed in worker processes, full phi set throughout: c = 10..220 and the
               768-bit primes
and on 2026-11-17 (the same box), with a costly winner kept in reserve (22eae52):
  yafu_c8a_4/  the 768-bit primes
  yafu_c8a_5/  c = 330..400, the entries of the frontier campaign, then c = 410, which is not in the table
bench_c8a.json holds, per entry, the recorded time and each pass; "best" is the latest pass that ran it."""
import json, os, re

HERE = os.path.dirname(os.path.abspath(__file__))
# README of ShortPrimalityProofs: (CPU seconds, wall seconds or None, cores or None, method)
H = 3600.0
RECORDED = {10: (1, None, None, "short.gp"), 20: (1, None, None, "short.gp"), 30: (1, None, None, "short.gp"),
            40: (2, None, None, "short.gp"), 50: (3, None, None, "short.gp"), 60: (3, None, None, "short.gp"),
            70: (45, None, None, "short.gp"), 80: (180, None, None, "short.gp"), 90: (350, None, None, "short.gp"),
            100: (800, None, None, "short.gp"), 110: (6300, None, None, "short.gp"), 120: (1900, None, None, "short.gp"),
            130: (1800, None, None, "short.gp"), 140: (13700, None, None, "short.gp"), 150: (11400, None, None, "short.gp"),
            160: (3500, 6 * 60, 10, "shortECPP.py"), 170: (25 * H, 4.2 * H, 6, "shortECPP.py"), 180: (5200, 9 * 60, 10, "shortECPP.py"),
            190: (15300, 26 * 60, 10, "shortECPP.py"), 200: (38 * H, 3.8 * H, 10, "shortECPP.py"), 210: (130 * H, 13 * H, 10, "shortECPP.py"),
            220: (21 * H, 1.3 * H, 16, "shortECPP.py"), 230: (59 * H, 3.7 * H, 16, "shortECPP.py"), 240: (310 * H, 1.6 * H, 192, "shortECPP.py"),
            250: (290 * H, 1.5 * H, 192, "shortECPP.py"), 260: (170 * H, 0.9 * H, 192, "shortECPP.py"), 270: (520 * H, 2.7 * H, 192, "shortECPP.py"),
            280: (500 * H, 2.6 * H, 192, "shortECPP.py"), 290: (380 * H, 2.0 * H, 192, "shortECPP.py"), 300: (380 * H, 2.0 * H, 192, "shortECPP.py"),
            310: (960 * H, 5.0 * H, 192, "shortECPP.py"), 320: (390 * H, 2.4 * H, 192, "shortECPP.py"),
            330: (1080 * H, 7.0 * H, 192, "shortECPP.py"), 340: (420 * H, 2.7 * H, 192, "shortECPP.py"),
            350: (900 * H, 5.6 * H, 192, "shortECPP.py"), 360: (990 * H, 6.0 * H, 192, "shortECPP.py"),
            370: (1650 * H, 10.5 * H, 192, "shortECPP.py"), 380: (70 * H, 1.1 * H, 192, "shortECPP.py"),
            390: (840 * H, 5.3 * H, 192, "shortECPP.py"), 400: (2340 * H, 12.8 * H, 192, "shortECPP.py")}
PASSES = ["yafu_c8a_1", "yafu_c8a_2", "yafu_c8a_3", "yafu_c8a_4", "yafu_c8a_5"]


def runs(d):
    out = {}
    if not os.path.exists(os.path.join(HERE, d, "progress")):
        return out
    for ln in open(os.path.join(HERE, d, "progress")):
        m = re.search(r"(\S+) done: RUN \S+ wall=([\d.]+) user=([\d.]+) sys=([\d.]+) verified=(\d)", ln)
        if m:
            tag = m.group(1)
            log = open(os.path.join(HERE, d, tag + ".log")).read()
            lv = re.search(r"chain complete: (\d+) levels", log)
            bits = re.search(r"proving p0 \((\d+) bits\)", log)
            out[tag] = {"wall": float(m.group(2)), "cpu": float(m.group(3)) + float(m.group(4)), "verified": m.group(5) == "1",
                        "levels": int(lv.group(1)) if lv else None, "bits": int(bits.group(1)) if bits else None,
                        "failed_builds": log.count("build failed")}
    return out


def main():
    P = [runs(d) for d in PASSES]
    rows = []
    for c in sorted(RECORDED):
        tag = "c%d" % c
        if not any(tag in p for p in P):
            continue                              # not run (yet)
        rec = RECORDED[c]
        row = {"c": c, "recorded_cpu": rec[0], "recorded_wall": rec[1], "recorded_cores": rec[2], "recorded_by": rec[3], "passes": {}}
        for name, p in zip(PASSES, P):
            if tag in p:
                row["passes"][name] = p[tag]
                row["best"] = dict(p[tag], set=name)
        rows.append(row)
    sets = {}
    for name, p in zip(PASSES, P):
        for pre in ("y768_", "y1030_"):
            v = [p[t] for t in sorted(p) if t.startswith(pre)]
            if v:
                w = sorted(x["wall"] for x in v)
                sets["%s/%s" % (name, pre.rstrip("_"))] = {"n": len(v), "wall_mean": sum(w) / len(w), "wall_median": w[len(w) // 2], "wall_max": w[-1],
                                                           "cpu_mean": sum(x["cpu"] for x in v) / len(v), "failed_builds": sum(x["failed_builds"] for x in v)}
    json.dump({"entries": rows, "sets": sets}, open(os.path.join(HERE, "bench_c8a.json"), "w"), indent=1)
    print("  c  bits lev |   recorded CPU      wall on cores |  now: wall      CPU  | CPU ratio  wall ratio")
    tw = tc = rc = 0.0
    for r in rows:
        b = r["best"]
        rw = ("%8.0f s on %3d" % (r["recorded_wall"], r["recorded_cores"])) if r["recorded_wall"] else "               -"
        print("%3d %5d %3d | %10.0f s  %s | %8.1f s %8.0f s | %8.0fx %s  %s%s"
              % (r["c"], b["bits"], b["levels"], r["recorded_cpu"], rw, b["wall"], b["cpu"], r["recorded_cpu"] / max(b["cpu"], 0.01),
                 ("%8.0fx" % (r["recorded_wall"] / b["wall"])) if r["recorded_wall"] else "        -", b["set"][-1], "" if b["verified"] else " NOT VERIFIED"))
        if r["c"] >= 160:
            tw += b["wall"]; tc += b["cpu"]; rc += r["recorded_cpu"]
    print("c = 160..320: %.0f s wall (%.1f min) and %.1f core-h in all; recorded %.0f core-h" % (tw, tw / 60, tc / 3600, rc / 3600))
    small = [r for r in rows if r["c"] <= 150]
    print("c = 10..150: %.1f s wall, %.0f CPU-s in all; recorded %.0f CPU-s" % (sum(r["best"]["wall"] for r in small), sum(r["best"]["cpu"] for r in small), sum(r["recorded_cpu"] for r in small)))
    for k, v in sets.items():
        print("%-22s %2d runs: wall mean %6.1f s, median %6.1f s, max %6.1f s; CPU mean %.2f core-h; failed builds %d" % (k, v["n"], v["wall_mean"], v["wall_median"], v["wall_max"], v["cpu_mean"] / 3600, v["failed_builds"]))


if __name__ == "__main__":
    main()
