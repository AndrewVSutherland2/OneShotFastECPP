#!/usr/bin/env python3
"""
Correctness test for the Cornacchia CM-discriminant search.

For a fixed prime p, run `cmsearch ... dump` to get, for every fundamental D=-d
with d <= MAXD, whether 4p = t^2 + d v^2 was found solvable (and the (t,v)).
Cross-check against PARI/GP's principal-form oracle: 4p = t^2 + d v^2 is solvable
with the CM parity iff p is represented by the principal form of discriminant -d
(qfbsolve(Qfb(1,1,(d+1)/4) or Qfb(1,0,d/4), p)).  Also verify 4p == t^2+d v^2
exactly for every reported solution.

Usage: python3 test_oracle.py [--p <decimal> | --pbits N] [--maxd 3000]
"""
import argparse, subprocess, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
CMSEARCH = HERE / "cmsearch"


def gp(script):
    r = subprocess.run(["gp", "-q"], input=script, capture_output=True, text=True, timeout=3600)
    if r.returncode != 0:
        sys.stderr.write(r.stderr); raise RuntimeError("gp failed")
    return r.stdout


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--p", type=str, default=None)
    ap.add_argument("--pbits", type=int, default=192)
    ap.add_argument("--maxd", type=int, default=3000)
    args = ap.parse_args()

    if args.p:
        primes = [int(args.p)]
    else:
        # one prime of each residue class mod 4 (the p==3 mod 4 fast-path and the
        # p==1 mod 4 Tonelli path are different code), both >> maxd
        b = args.pbits
        primes = [int(gp(f"p=2^{b-1}+314159265; while(p%4!={r}||!isprime(p),p++); print(p)\n").strip())
                  for r in (1, 3)]

    overall_ok = True
    for p in primes:
        overall_ok &= check_one(p, args.maxd)
    return 0 if overall_ok else 1


def check_one(p, maxd):
    print(f"\np = {p}  ({p.bit_length()} bits, {p % 4} mod 4)")

    # run the search in dump mode
    out = subprocess.run([str(CMSEARCH), f"p={p}", f"maxd={maxd}", "dump"],
                         capture_output=True, text=True, timeout=600).stdout
    mine = {}   # d -> (solvable, t, v)
    for line in out.splitlines():
        f = line.split()
        d = int(f[0]); sol = int(f[1])
        mine[d] = (sol, int(f[2]), int(f[3])) if sol else (0, 0, 0)

    # arithmetic check of every reported solution
    arith_bad = []
    for d, (sol, t, v) in mine.items():
        if sol and t * t + d * v * v != 4 * p:
            arith_bad.append(d)

    # PARI oracle for solvability, one statement per line
    lines = []
    for d in mine:
        if d % 4 == 3:
            q = f"Qfb(1,1,{(d+1)//4})"
        else:
            q = f"Qfb(1,0,{d//4})"
        lines.append(f'print("O {d} ", if(qfbsolve({q},{p})==[],0,1))')
    oracle = {}
    for line in gp("\n".join(lines) + "\n").splitlines():
        if line.startswith("O "):
            _, d, s = line.split()
            oracle[int(d)] = int(s)

    mism = [d for d in mine if mine[d][0] != oracle.get(d)]

    n = len(mine)
    n_solv = sum(1 for d in mine if mine[d][0])
    print(f"fundamental D tested: {n}   solvable (mine): {n_solv}")
    print(f"arithmetic check failures (4p != t^2+d v^2): {len(arith_bad)}")
    print(f"solvability mismatches vs PARI oracle: {len(mism)}")
    if arith_bad:
        print("  arith-bad d:", arith_bad[:20])
    if mism:
        for d in mism[:20]:
            print(f"  d={d}: mine={mine[d][0]} pari={oracle.get(d)}")
    ok = not arith_bad and not mism
    print("OK" if ok else "FAIL")
    return ok


if __name__ == "__main__":
    sys.exit(main())
