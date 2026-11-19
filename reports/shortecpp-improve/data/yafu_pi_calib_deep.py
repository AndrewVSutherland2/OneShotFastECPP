#!/usr/bin/env python3
"""Measure the per-curve success probability of yafu's AVX-ECM by factor size: the deep sweeps
(B1 = 3e6, 1e7, 3e7) at 102 to 166 bits, with ten times the budget of yafu_pi_calib.py, whose cells there
have 0 to 11 hits each.  Everything else as in that script:

Same design as pi_calib.py (gmp-ecm): for each B1 and each factor size b in bits,
take 96 random b-bit primes f, form N = f*g with g a fixed 192-bit prime that is
never found, and run curves until f appears or the cap is reached.  The success
probability of a curve does not depend on the size of N, so a small N keeps the
measurement cheap.

yafu runs curves in batches of 8 (one AVX-512 vector) and stops after the batch
in which the factor appears.  With -v -v -v it reports every lane that found the
factor, so the estimate is (lanes that hit) / (curves run), the maximum
likelihood estimate under this stopping rule.

Engine: yafu 3.1.9, `-avxecm_stg2 1` (FFT stage 2, B2 a function of B1 only),
default curve family (twisted Edwards, a = -1).

usage: yafu_pi_calib.py <yafu binary> <out.json> [core-seconds per cell = 1000] [seed = 20261115]
"""
import concurrent.futures as cf
import json
import os
import random
import re
import subprocess
import sys
import time

YAFU = sys.argv[1]
OUT = sys.argv[2]
BUDGET = float(sys.argv[3]) if len(sys.argv) > 3 else 1000.0
SEED = int(sys.argv[4]) if len(sys.argv) > 4 else 20261115
NCPU = int(os.environ.get("CALIB_THREADS", os.cpu_count()))
NF = 96
MAXBATCH = 750
FOUND = re.compile(r"found factor (\d+) in stage (\d) in thread \d+, vec position (\d+)")
B2RE = re.compile(r"FFT stage 2: B2 = (\d+)")


def is_prime(n, rng):
    if n < 2:
        return False
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % q == 0:
            return n == q
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for _ in range(20):
        a = rng.randrange(2, n - 1)
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def rand_prime(bits, rng):
    while True:
        n = rng.getrandbits(bits) | (1 << (bits - 1)) | 1
        if is_prime(n, rng):
            return n


def run(b1, n, f, batches):
    """(curves run, lanes that found f, lanes that found it in stage 1, seconds, B2)."""
    t = time.time()
    r = subprocess.run([YAFU, "-B1ecm", str(b1), "-threads", "1", "-avxecm_stg2", "1", "-v", "-v", "-v",
                        "-logfile", "/dev/null", "-session", "/dev/null"],
                       input="ecm(%d,%d)\n" % (n, 8 * batches), capture_output=True, text=True, cwd="/tmp")
    dt = time.time() - t
    out = r.stdout.replace("\r", "\n")
    ran = 8 * out.count("Commencing curves ")
    lanes, lanes1 = set(), set()
    for m in FOUND.finditer(out):
        if f and int(m.group(1)) % f == 0:
            lanes.add(int(m.group(3)))
            if m.group(2) == "1":
                lanes1.add(int(m.group(3)))
    b2 = B2RE.search(out)
    return ran, len(lanes), len(lanes1), dt, int(b2.group(1)) if b2 else 0


def main():
    rng = random.Random(SEED)
    g = rand_prime(192, rng)
    g2 = rand_prime(160, rng)
    # the deep sweeps only, at the factor sizes where the first calibration had a handful of hits per cell
    tiers = [(3000000, range(102, 151, 4)), (10000000, range(102, 159, 4)), (30000000, range(106, 167, 4))]
    tasks = []
    cells = {}
    for b1, brange in tiers:
        ran, _, _, dt, b2 = run(b1, g * g2, 0, 1)        # one batch on a number with no factor in reach
        assert ran == 8, (b1, ran)
        cap = max(1, min(MAXBATCH, int(BUDGET / (NF * max(dt, 1e-3)))))
        print("B1=%d: B2=%d, %.3fs per batch of 8, %d primes x up to %d batches per cell" % (b1, b2, dt, NF, cap), flush=True)
        for b in brange:
            cells[(b1, b)] = {"method": "ecm", "engine": "yafu", "b1": b1, "b2": b2, "bits": b, "primes": 0, "found": 0,
                              "curves": 0, "hits": 0, "hits_stage1": 0, "cpu_s": 0.0, "cap": 8 * cap}
            for _ in range(NF):
                tasks.append((b1, b, rand_prime(b, rng), cap))
    rng.shuffle(tasks)
    tasks.sort(key=lambda t: -t[0])                       # longest tiers first, shuffled within a tier
    print("%d tasks on %d threads" % (len(tasks), NCPU), flush=True)
    t0 = time.time()
    done = 0

    def work(t):
        b1, b, f, cap = t
        return t, run(b1, f * g, f, cap)

    def save():
        json.dump(list(cells.values()), open(OUT + ".tmp", "w"), indent=0)
        os.replace(OUT + ".tmp", OUT)

    with cf.ThreadPoolExecutor(max_workers=NCPU) as ex:
        for t, (ran, hits, hits1, dt, _) in ex.map(work, tasks):
            c = cells[t[:2]]
            c["primes"] += 1
            c["found"] += 1 if hits else 0
            c["curves"] += ran
            c["hits"] += hits
            c["hits_stage1"] += hits1
            c["cpu_s"] += dt
            done += 1
            if done % 500 == 0:
                print("[%6.0fs] %d/%d" % (time.time() - t0, done, len(tasks)), flush=True)
                save()
    save()
    for c in cells.values():
        print("B1=%-9d %3d bits: %4d hits / %7d curves = %.5f  (%d of %d primes, %.0f s)"
              % (c["b1"], c["bits"], c["hits"], c["curves"], c["hits"] / max(c["curves"], 1), c["found"], c["primes"], c["cpu_s"]))
    print("CALIB_DONE %.0fs" % (time.time() - t0), flush=True)


if __name__ == "__main__":
    main()
