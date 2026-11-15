#!/usr/bin/env python3
"""Merge runs of yafu_pi_calib.py (independent curves; cells are sums of counts) into one table.

usage: yafu_pi_merge.py <out.json> <in1.json> [<in2.json> ...]
"""
import json
import sys

cells = {}
for fn in sys.argv[2:]:
    for c in json.load(open(fn)):
        k = (c["b1"], c["bits"])
        if k not in cells:
            cells[k] = dict(c, runs=1)
            continue
        d = cells[k]
        for f in ("primes", "found", "curves", "hits", "hits_stage1", "cpu_s"):
            d[f] += c[f]
        d["cap"] = max(d["cap"], c["cap"])
        d["runs"] += 1
        if d["b2"] != c["b2"]:
            sys.exit("B2 differs between runs at B1=%d" % c["b1"])
out = [cells[k] for k in sorted(cells) if cells[k]["curves"] > 0]
json.dump(out, open(sys.argv[1], "w"), indent=0)
b1s = sorted({k[0] for k in cells})
for b1 in b1s:
    row = [c for c in out if c["b1"] == b1]
    good = [c for c in row if c["hits"] >= 8]
    print("B1=%-9d %2d cells, %2d with 8 hits or more (%d..%d bits), %8d curves, %6.0f core-s"
          % (b1, len(row), len(good), good[0]["bits"] if good else 0, good[-1]["bits"] if good else 0,
             sum(c["curves"] for c in row), sum(c["cpu_s"] for c in row)))
