#!/usr/bin/env python3
"""bench_pm1_ratio.py -- measure the per-B1 cost of one P-1 run against one
ECM curve, the ratio behind the b1 // 8 effort weight in shortECPP.py's
ladder_rounds.

Inputs are deterministic (seed 7): a prime semiprime per bit size, so neither
method finds a factor and both run their full stage budgets.  The ECM binary
and environment are resolved exactly as the prover resolves them (imported
from shortECPP).  Reference results, gmp-ecm 7.0.6 on the 16-core AMD Ryzen
AI Max+ 395 development box, 2026-08-25:

     900 bits, B1=500000:   ecm 1 curve  1.3 s,  pm1 0.2 s  ->  5.6x
     900 bits, B1=10000000: ecm 1 curve 24.4 s,  pm1 3.6 s  ->  6.9x
    1330 bits, B1=500000:   ecm 1 curve  3.3 s,  pm1 0.4 s  ->  8.2x
    1330 bits, B1=10000000: ecm 1 curve 67.1 s,  pm1 6.4 s  -> 10.5x

usage: python3 ecpp/bench_pm1_ratio.py [bits ...]   (default: 900 1330)
"""
import random, subprocess, sys, time

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from shortECPP import ECM_BIN, ECM_ENV, is_prp


def randprime(bits, rng):
    while True:
        c = rng.getrandbits(bits) | (1 << (bits - 1)) | 1
        if is_prp(c):
            return c


def bench(n, args):
    """Time one run and validate it: a successful no-factor run exits 0 and,
    under -q, echoes exactly the input composite.  gmp-ecm exit codes are
    semantic (e.g. 14 = factor found), so check=True would be wrong; any
    deviation here -- an error, an unsupported option, or an (astronomically
    unlikely) factor of the prime semiprime -- aborts loudly instead of
    contaminating the ratio."""
    t = time.time()
    r = subprocess.run([ECM_BIN, "-q"] + args, input=str(n) + "\n",
                       capture_output=True, text=True, env=ECM_ENV)
    dt = time.time() - t
    if r.returncode != 0 or r.stdout.split() != [str(n)]:
        sys.exit("bench: %r did not complete unfactored (exit %d, stdout %r, "
                 "stderr %r)" % (args, r.returncode, r.stdout[:200],
                                 r.stderr[:200]))
    return dt


def main():
    sizes = [int(a) for a in sys.argv[1:]] or [900, 1330]
    rng = random.Random(7)
    print(f"ecm binary: {ECM_BIN}")
    for bits in sizes:
        n = randprime(bits // 2, rng) * randprime(bits - bits // 2, rng)
        for b1 in (500000, 10000000):
            te = bench(n, ["-c", "1", str(b1)])
            tp = bench(n, ["-pm1", "-c", "1", str(b1)])
            print(f"{bits} bits, B1={b1}: ecm 1 curve {te:.1f}s, "
                  f"pm1 {tp:.1f}s -> ratio {te / tp:.1f}x")


if __name__ == "__main__":
    main()
