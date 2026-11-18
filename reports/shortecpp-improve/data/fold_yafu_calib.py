#!/usr/bin/env python3
"""fold_yafu_calib.py -- the prover's success table for yafu's ECM, ecpp/ecm_pi_yafu.json,
from the two calibrations in this directory:

  yafu_pi_calib.json       182 cells, B1 = 2e3 .. 3e7, 6.06M curves (yafu_pi_calib.py)
  yafu_pi_calib_deep.json  44 cells at B1 = 3e6, 1e7, 3e7, 598,344 curves on factors of
                           102 .. 174 bits (yafu_pi_calib_deep.py), where the first
                           calibration had 768 curves a cell and few or no hits

Cells of the same (method, B1, factor bits) are pooled: curves, hits, primes and CPU
seconds add up.  shortsim.EcmModel uses a cell once it has 8 hits and follows the
two-stage Dickman curve, rescaled to the last such cell, beyond them.

usage: fold_yafu_calib.py [write]     (without "write": print the comparison only)
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "..", "..", "ecpp", "ecm_pi_yafu.json")
ADD = ("primes", "found", "curves", "hits", "hits_stage1", "cpu_s", "runs")


def main():
    base = json.load(open(os.path.join(HERE, "yafu_pi_calib.json")))
    deep = json.load(open(os.path.join(HERE, "yafu_pi_calib_deep.json")))
    cells = {}
    order = []
    for src in (base, deep):
        for c in src:
            k = (c["method"], c["b1"], c["bits"])
            if k not in cells:
                cells[k] = dict(c)
                cells[k].setdefault("runs", 1)
                order.append(k)
                continue
            m = cells[k]
            if m["b2"] != c["b2"]:
                sys.exit("B2 differs in cell %r" % (k,))
            for f in ADD:
                m[f] = m.get(f, 0) + c.get(f, 1 if f == "runs" else 0)
            m["cap"] = max(m.get("cap", 0), c.get("cap", 0))
    out = [cells[k] for k in sorted(order, key=lambda k: (k[1], k[2]))]
    old = {(c["method"], c["b1"], c["bits"]): c for c in base}
    print("%9s %5s  %21s  %21s" % ("B1", "bits", "first: hits/curves", "pooled: hits/curves"))
    for c in out:
        k = (c["method"], c["b1"], c["bits"])
        o = old.get(k)
        if o is None or o["curves"] != c["curves"]:
            print("%9d %5d  %21s  %21s" % (c["b1"], c["bits"],
                                           "%d/%d" % (o["hits"], o["curves"]) if o else "-",
                                           "%d/%d" % (c["hits"], c["curves"])))
    print("%d cells, %d curves" % (len(out), sum(c["curves"] for c in out)))
    if "write" in sys.argv[1:]:
        with open(OUT, "w") as f:
            json.dump(out, f, indent=0)
        print("wrote", os.path.normpath(OUT))


if __name__ == "__main__":
    main()
