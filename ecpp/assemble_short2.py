#!/usr/bin/env python3
"""assemble_short2.py -- merge the v2 migration outputs into certs/short2/certs.csv
(original chain order), validating every chain under BOTH verifiers.

usage: python3 ecpp/assemble_short2.py
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
from vshort2 import _verify as v2_verify
from repair_short2 import first_bad_level
import importlib.util
_sp = importlib.util.spec_from_file_location(
    "v1", os.path.join(os.path.dirname(__file__), "..", "verifier-batching-pr", "vsmallECPP.py"))
v1mod = importlib.util.module_from_spec(_sp); _sp.loader.exec_module(v1mod)

ROOT = os.path.join(os.path.dirname(__file__), "..")
V1CSV = os.path.join(ROOT, "certs/short2/certs_v1.csv")
REST = os.path.join(ROOT, "work/short2repair2/certs.csv")
P698 = os.path.join(ROOT, "work/short2repair/nextprime1e210v2.txt")
OUT = os.path.join(ROOT, "certs/short2/certs.csv")


def load_csv(path):
    return [[int(t) for t in line.strip().split(',')] for line in open(path) if line.strip()]


def main():
    v1 = load_csv(V1CSV)
    rest = {c[0]: c for c in load_csv(REST)} if os.path.exists(REST) else {}
    p698chain = None
    if os.path.exists(P698):
        p698chain = [int(t) for t in open(P698).read().split()]
    out, bad = [], 0
    for ints in v1:
        n = ints[0].bit_length()
        if first_bad_level(ints) is None:
            chain = ints
            src = "kept"
        elif n == 698:
            chain = p698chain
            src = "short_prove v2"
        else:
            chain = rest.get(ints[0])
            src = "repaired"
        if chain is None:
            print(f"n={n:4d}: MISSING ({src})")
            bad += 1
            continue
        okv2 = v2_verify(chain)
        okv1 = v1mod.verify(chain)
        levs = (len(chain) - 1) // 3
        print(f"n={n:4d}: {src:14s} levels={levs} v2={okv2[0]} v1={okv1}"
              + ("" if okv2[0] and okv1 else f"  <-- {okv2}"))
        if not (okv2[0] and okv1):
            bad += 1
            continue
        out.append(chain)
    if bad:
        print(f"{bad} chains missing/invalid -- certs.csv NOT written")
        sys.exit(1)
    with open(OUT, "w") as f:
        for c in out:
            f.write(",".join(map(str, c)) + "\n")
    print(f"wrote {len(out)} chains to {OUT}")


if __name__ == "__main__":
    main()
