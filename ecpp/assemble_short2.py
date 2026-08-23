#!/usr/bin/env python3
"""assemble_short2.py -- validate (default) or re-assemble the migrated short-ECPP table.

Default (no arguments): deterministically validate the TRACKED table
certs/short2/certs.csv -- every chain must verify under the revised-format verifier
(ecpp/vshort2.py) and contain no level violating the format constraints.  Exits
nonzero on any failure, so it doubles as a consistency gate.

--merge: the historical migration mode (August 2026): merge the repair outputs from
work/short2repair*/ (produced by repair_short2.py and shortECPP.py runs; those
work directories are transient and not tracked) into certs/short2/certs.csv,
validating every chain along the way.  Kept for provenance; a fresh checkout
cannot run it without re-running the migration searches.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
from vshort2 import _verify as v2_verify
from repair_short2 import first_bad_level

ROOT = os.path.join(os.path.dirname(__file__), "..")
V1CSV = os.path.join(ROOT, "certs/short2/certs_v1.csv")
RESTS = [os.path.join(ROOT, d, "certs.csv")
         for d in ("work/short2repair3", "work/short2repair4", "work/short2repair2")]
PROVED = [os.path.join(ROOT, "work/short2repair", f)
          for f in ("nextprime1e210v2.txt", "nextprime1e170v2.txt")]
OUT = os.path.join(ROOT, "certs/short2/certs.csv")


def load_csv(path):
    return [[int(t) for t in line.strip().split(',')] for line in open(path) if line.strip()]


def validate():
    bad = 0
    chains = load_csv(OUT)
    expected = [c[0] for c in load_csv(V1CSV)]    # the authoritative prime list, in order
    got = [c[0] for c in chains]
    if got != expected:
        print(f"table mismatch: {len(got)} chains present, expected the "
              f"{len(expected)} primes of certs_v1.csv in the same order")
        for p in expected:
            if p not in got:
                print(f"  missing: 10^{len(str(p))-1} chain")
        for i, p in enumerate(got):
            if p not in expected or got.index(p) != i:
                print(f"  extra/duplicate/misplaced at row {i}: 10^{len(str(p))-1}")
        sys.exit(1)
    for chain in chains:
        n = chain[0].bit_length()
        okv2 = v2_verify(chain)
        lev = first_bad_level(chain)
        print(f"n={n:4d}: levels={(len(chain)-1)//3} verify={okv2[0]} format-clean={lev is None}"
              + ("" if okv2[0] and lev is None else f"  <-- {okv2} bad-level={lev}"))
        if not (okv2[0] and lev is None):
            bad += 1
    print(f"{len(chains)} chains checked, {bad} failures")
    sys.exit(1 if bad else 0)


def merge():
    v1 = load_csv(V1CSV)
    rest = {}
    for path in RESTS:
        if os.path.exists(path):
            for c in load_csv(path):
                rest.setdefault(c[0], c)
    for path in PROVED:                       # shortECPP outputs win over drivers
        if os.path.exists(path):
            c = [int(t) for t in open(path).read().split()]
            rest[c[0]] = c
    out, bad = [], 0
    for ints in v1:
        n = ints[0].bit_length()
        if first_bad_level(ints) is None:
            chain = ints
            src = "kept"
        else:
            chain = rest.get(ints[0])
            src = "repaired"
        if chain is None:
            print(f"n={n:4d}: MISSING ({src})")
            bad += 1
            continue
        okv2 = v2_verify(chain)
        levs = (len(chain) - 1) // 3
        print(f"n={n:4d}: {src:14s} levels={levs} v2={okv2[0]}"
              + ("" if okv2[0] else f"  <-- {okv2}"))
        if not okv2[0]:
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
    if len(sys.argv) > 1 and sys.argv[1] == "--merge":
        merge()
    else:
        validate()
