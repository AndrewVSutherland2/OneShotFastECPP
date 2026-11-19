#!/usr/bin/env python3
"""Cost per curve of yafu's AVX-ECM as the prover uses it, by B1 and size of the number.

The prover hands a yafu process a list of numbers (one expression ecm(n, c) per
line) so that a process lasts a second or more; the cost of a curve is then the
wall time of such a process divided by the curves it ran, start-up included.
Here every hardware thread runs one such process at the same time on numbers
with no factor in reach (products of two primes of half the size), and the
median over the processes is recorded.

usage: yafu_cost2.py <yafu binary> <out.json> [sizes=260,364,...] [deep=1]
"""
import concurrent.futures as cf
import json
import os
import random
import subprocess
import sys
import time

YAFU = sys.argv[1]
OUT = sys.argv[2]
KV = dict(a.split("=", 1) for a in sys.argv[3:] if "=" in a)
SIZES = [int(x) for x in KV.get("sizes", "260,364,520,676,780,884,1030,1196,1352").split(",")]
DEEP = KV.get("deep", "1") == "1"
NCPU = int(KV.get("threads", os.cpu_count()))
# B1 -> curves per expression (the sweep the prover's ladder uses)
SWEEP = {2000: 16, 11000: 24, 50000: 16, 250000: 16, 1000000: 16, 3000000: 16, 10000000: 8, 30000000: 8}


def is_prime(n, rng):
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


def run(nums, b1, curves):
    t = time.time()
    subprocess.run([YAFU, "-B1ecm", str(b1), "-threads", "1", "-avxecm_stg2", "1", "-aprcl_p", "1",
                    "-logfile", "/dev/null", "-session", "/dev/null"],
                   input="".join("ecm(%d,%d)\n" % (n, curves) for n in nums), capture_output=True, text=True, cwd="/tmp")
    return time.time() - t


def main():
    rng = random.Random(11)
    res = []
    t = time.time()
    subprocess.run([YAFU, "-threads", "1", "-logfile", "/dev/null", "-session", "/dev/null"], input="2+2\n",
                   capture_output=True, text=True, cwd="/tmp")
    print("yafu start-up alone: %.3f s" % (time.time() - t), flush=True)
    for bits in SIZES:
        hard = [rand_prime(bits // 2, rng) * rand_prime(bits - bits // 2, rng) for _ in range(32)]
        for b1 in sorted(SWEEP):
            if b1 > 1000000 and not (DEEP and bits == 1030):
                continue
            c = SWEEP[b1]
            # one expression alone on an idle machine sizes the list
            t1 = max(run(hard[:1], b1, c) - 0.03, 1e-3)
            k = max(1, min(32, int(1.5 / t1)))
            with cf.ThreadPoolExecutor(max_workers=NCPU) as ex:
                ts = sorted(ex.map(lambda i: run([hard[(i + j) % 32] for j in range(k)], b1, c), range(NCPU)))
            med = ts[len(ts) // 2]
            rec = {"bits": bits, "b1": b1, "curves": c, "tails_per_process": k, "process_s": med,
                   "sec_per_curve": med / (k * c), "idle_single_s": t1, "threads": NCPU}
            res.append(rec)
            print(rec, flush=True)
            json.dump(res, open(OUT, "w"), indent=0)
    print("YCOST2_DONE", flush=True)


if __name__ == "__main__":
    main()
