#!/usr/bin/env python3
"""yafu's AVX-ECM (AVX-512 IFMA build) against gmp-ecm 7.0.5 on 1030-bit inputs.

Two measurements, both with every core of the box busy (one single-threaded
process per core, as the prover runs gmp-ecm):

  success  for each B1 and factor size b: 96 numbers N = f*g of 1030 bits with f a
           random b-bit prime; each engine runs one batch of 8 curves on each N.
           Recorded: whether f was found, and that every factor reported divides N.
  cost     for each B1: 192 numbers with no factor in reach; each engine runs
           exactly 8 curves on one of them per core.  Recorded: wall seconds.

usage: ytest.py <outdir> [stg2 mode for yafu: 1 = FFT with gmp-ecm's B2 (default), 2 = FFT with B2 from memory, 0 = pairing]
"""
import concurrent.futures as cf
import json
import os
import random
import re
import resource
import shutil
import subprocess
import sys
import tempfile
import time

OUT = sys.argv[1]
STG2 = sys.argv[2] if len(sys.argv) > 2 else "1"
YAFU = os.path.expanduser("~/yb/yafu/yafu")
BENCH = os.path.expanduser("~/ecmbench")
ECM = BENCH + "/inst/deb/usr/bin/ecm"
ENV = dict(os.environ)
ENV["LD_LIBRARY_PATH"] = BENCH + "/inst/deb/usr/lib/x86_64-linux-gnu"
NCPU = os.cpu_count()
NBITS = 1030
NF = 96
CELLS = [(50000, (54, 62, 70)), (250000, (62, 74, 86)), (1000000, (70, 82, 94))]


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
    for _ in range(16):
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


def run_yafu(n, b1, curves):
    d = tempfile.mkdtemp(prefix="y", dir=OUT + "/tmp")
    t = time.time()
    r = subprocess.run([YAFU, "-B1ecm", str(b1), "-threads", "1", "-avxecm_stg2", STG2, "-v"],
                       input="ecm(%d,%d)\n" % (n, curves), capture_output=True, text=True, cwd=d)
    dt = time.time() - t
    shutil.rmtree(d, ignore_errors=True)
    facs = [int(x) for x in re.findall(r"found (?:prp|c|p)\d+ factor = (\d+)", r.stdout)]
    m = re.search(r"(\d+)/(\d+) curves on C\d+ @ B1=(\d+), B2=([0-9e+.x*]+)", r.stdout)
    b2 = m.group(4) if m else ""
    return facs, dt, b2, r.stdout


def run_gmp(n, b1, curves):
    t = time.time()
    r = subprocess.run([ECM, "-c", str(curves), str(b1)], input="%d\n" % n, capture_output=True, text=True, env=ENV)
    dt = time.time() - t
    facs = [int(x) for x in re.findall(r"Factor found in step \d: (\d+)", r.stdout)]
    m = re.search(r"Using B1=\d+, B2=(\d+)", r.stdout)
    ran = len(re.findall(r"^Using B1=", r.stdout, re.M))
    return facs, dt, (m.group(1) if m else ""), r.stdout, ran


def main():
    os.makedirs(OUT + "/tmp", exist_ok=True)
    rng = random.Random(20261115)
    tasks = []
    for b1, bl in CELLS:
        for b in bl:
            g = rand_prime(NBITS - b, rng)
            for _ in range(NF):
                f = rand_prime(b, rng)
                for eng in ("yafu", "gmp"):
                    tasks.append(("success", eng, b1, b, f, f * g))
        hard = [rand_prime(NBITS // 2, rng) * rand_prime(NBITS - NBITS // 2, rng) for _ in range(8)]
        for i in range(NCPU):
            for eng in ("yafu", "gmp"):
                tasks.append(("cost", eng, b1, 0, 0, hard[i % 8]))
    # the cost tasks of one engine and B1 run together on all cores; success tasks fill in afterwards
    order = {"cost": 0, "success": 1}
    tasks.sort(key=lambda t: (order[t[0]], -t[2], t[1]))
    print("%d tasks, yafu stage 2 mode %s" % (len(tasks), STG2), flush=True)
    res = []
    samples = {}
    t0 = time.time()

    def work(t):
        kind, eng, b1, b, f, n = t
        if eng == "yafu":
            facs, dt, b2, out = run_yafu(n, b1, 8)
            ran = 8
        else:
            facs, dt, b2, out, ran = run_gmp(n, b1, 8)
        bad = [x for x in facs if x <= 1 or n % x]
        found = bool(f) and any(x % f == 0 for x in facs if x not in bad)
        return {"kind": kind, "eng": eng, "b1": b1, "bits": b, "found": found, "nfac": len(facs), "bad": len(bad),
                "secs": dt, "b2": b2, "ran": ran}, out

    # cost groups first, one (engine, B1) group at a time so that each runs on a full box by itself
    groups = {}
    for t in tasks:
        groups.setdefault((t[0], t[1], t[2]) if t[0] == "cost" else ("success",), []).append(t)
    for key in sorted(groups, key=lambda k: (k[0] != "cost", k)):
        grp = groups[key]
        tg = time.time()
        with cf.ThreadPoolExecutor(max_workers=NCPU) as ex:
            for r, out in ex.map(work, grp):
                res.append(r)
                k = (r["kind"], r["eng"], r["b1"])
                if k not in samples:
                    samples[k] = out[-1500:]
        print("[%5.0fs] %s: %d tasks in %.1fs" % (time.time() - t0, key, len(grp), time.time() - tg), flush=True)
        json.dump({"stg2": STG2, "res": res, "samples": {"|".join(map(str, k)): v for k, v in samples.items()}},
                  open(OUT + "/ytest_stg2_%s.json" % STG2, "w"))
    # summary
    print("cost, core-seconds per curve with all cores busy (8 curves per process):")
    for b1, _ in CELLS:
        row = []
        for eng in ("gmp", "yafu"):
            v = sorted(r["secs"] / 8 for r in res if r["kind"] == "cost" and r["eng"] == eng and r["b1"] == b1)
            b2 = [r["b2"] for r in res if r["kind"] == "cost" and r["eng"] == eng and r["b1"] == b1][0]
            row.append("%s %.3f (median of %d, B2=%s)" % (eng, v[len(v) // 2], len(v), b2))
        print("  B1=%d: %s" % (b1, "   ".join(row)))
    print("success, numbers out of %d whose factor was found in one batch of 8 curves:" % NF)
    for b1, bl in CELLS:
        for b in bl:
            row = []
            for eng in ("gmp", "yafu"):
                v = [r for r in res if r["kind"] == "success" and r["eng"] == eng and r["b1"] == b1 and r["bits"] == b]
                row.append("%s %d (bad factors %d)" % (eng, sum(r["found"] for r in v), sum(r["bad"] for r in v)))
            print("  B1=%d, %d-bit factor: %s" % (b1, b, "   ".join(row)))
    print("YTEST_DONE %.0fs" % (time.time() - t0), flush=True)


if __name__ == "__main__":
    _soft, _hard = resource.getrlimit(resource.RLIMIT_NOFILE)
    resource.setrlimit(resource.RLIMIT_NOFILE, (min(_hard, 65536), _hard))
    main()
