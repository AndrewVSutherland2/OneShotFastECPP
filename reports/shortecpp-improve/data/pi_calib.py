#!/usr/bin/env python3
"""Measure gmp-ecm's per-curve success probability by factor size.

For each tier (method, B1, default B2, default parametrisation -- exactly what
shortECPP.py runs) and each factor size b in bits, take many random b-bit
primes f, form N = f*g with g a fixed 192-bit prime that is never found, and
run curves until f appears or the curve cap is hit.  The estimate is
successes / curves run.  P-1 is deterministic, so it gets one run per prime.

Output: out/pi_calib.json, a list of cells
    {method, b1, bits, primes, curves, hits, cpu_s}
"""
import concurrent.futures as cf
import json
import os
import random
import re
import resource
import subprocess
import sys
import time

R = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(R, "out")
ECM = R + "/inst/deb/usr/bin/ecm"
ENV = dict(os.environ)
ENV["LD_LIBRARY_PATH"] = R + "/inst/deb/usr/lib/x86_64-linux-gnu"
NCPU = os.cpu_count()
BUDGET = float(sys.argv[1]) if len(sys.argv) > 1 else 3000.0   # core-seconds per cell
NF = 96


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


def run(method, b1, n, f, curves):
    cmd = [ECM] + (["-pm1"] if method == "pm1" else []) + ["-c", str(curves), str(b1)]
    t = time.time()
    r = subprocess.run(cmd, input=str(n) + "\n", capture_output=True, text=True, env=ENV)
    dt = time.time() - t
    ran = len(re.findall(r"^Using B1=", r.stdout, re.M))
    hit = 0
    for m in re.finditer(r"Factor found in step \d: (\d+)", r.stdout):
        if int(m.group(1)) % f == 0:
            hit = 1
    return ran, hit, dt


def main():
    rng = random.Random(20261114)
    g = rand_prime(192, rng)
    tiers = [("ecm", 2000, range(34, 95, 4)), ("ecm", 11000, range(34, 115, 4)),
             ("ecm", 50000, range(38, 131, 4)), ("ecm", 250000, range(46, 143, 4)),
             ("ecm", 1000000, range(58, 151, 4)), ("ecm", 3000000, range(66, 159, 4)),
             ("pm1", 500000, range(34, 131, 4)), ("pm1", 10000000, range(34, 151, 4))]
    g2 = rand_prime(160, rng)
    tasks = []
    cells = {}
    for method, b1, brange in tiers:
        # time one curve on a number with no findable factor
        ran, _, dt = run(method, b1, g * g2, g, 1)
        tcurve = max(dt, 1e-3)
        if method == "pm1":
            nf, cap = min(4000, max(400, int(BUDGET / tcurve))), 1
        else:
            nf, cap = NF, max(8, min(6000, int(BUDGET / (NF * tcurve))))
        print("%s B1=%d: %.3fs/curve, %d primes x up to %d curves per cell" % (method, b1, tcurve, nf, cap), flush=True)
        for b in brange:
            cells[(method, b1, b)] = {"method": method, "b1": b1, "bits": b, "primes": 0, "curves": 0,
                                      "hits": 0, "cpu_s": 0.0, "cap": cap}
            for _ in range(nf):
                f = rand_prime(b, rng)
                tasks.append((method, b1, b, f, cap))
    rng.shuffle(tasks)
    tasks.sort(key=lambda t: -t[1])                 # longest tiers first, shuffled within a tier
    print("%d tasks" % len(tasks), flush=True)
    t0 = time.time()
    done = 0

    def work(t):
        method, b1, b, f, cap = t
        return t, run(method, b1, f * g, f, cap)
    with cf.ThreadPoolExecutor(max_workers=NCPU) as ex:
        for t, (ran, hit, dt) in ex.map(work, tasks):
            c = cells[t[:3]]
            c["primes"] += 1
            c["curves"] += ran
            c["hits"] += hit
            c["cpu_s"] += dt
            done += 1
            if done % 2000 == 0:
                print("[%5.0fs] %d/%d" % (time.time() - t0, done, len(tasks)), flush=True)
                json.dump(list(cells.values()), open(os.path.join(OUT, "pi_calib.json.tmp"), "w"))
                os.replace(os.path.join(OUT, "pi_calib.json.tmp"), os.path.join(OUT, "pi_calib.json"))
    json.dump(list(cells.values()), open(os.path.join(OUT, "pi_calib.json"), "w"), indent=0)
    print("PI_DONE %.0fs" % (time.time() - t0), flush=True)


if __name__ == "__main__":
    _soft, _hard = resource.getrlimit(resource.RLIMIT_NOFILE)
    resource.setrlimit(resource.RLIMIT_NOFILE, (min(_hard, 65536), _hard))   # 192 workers x several pipes each
    main()
