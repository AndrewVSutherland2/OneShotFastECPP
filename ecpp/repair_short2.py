#!/usr/bin/env python3
"""repair_short2.py -- migrate a published short-ECPP table to the revised (August
2026) format.  For each chain in certs.csv: find the first level whose DATA violates
the revised constraints (a filler prime > B = ceil(n^2/log2 n), or a filler radical
with log2 rad(m) > ceil(n/log2 n)), keep the valid prefix, and re-find the chain from
that level's modulus with short2.gp
(first-winning-worker parallelism, kills by Popen handle only).  Every output chain is
validated with BOTH vshort2.py (v2) and the original vsmallECPP.py (v1 => superset).

usage: python3 ecpp/repair_short2.py <certs_v1.csv> <outdir>
Writes <outdir>/certs.csv (same chain order), logs under <outdir>/log/.
"""
import os, sys, math, time, subprocess, random
from math import gcd, isqrt, log2

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
from vshort2 import _verify as v2_verify, sieve_primes
import importlib.util
_sp = importlib.util.spec_from_file_location(
    "v1", os.path.join(os.path.dirname(__file__), "..", "verifier-batching-pr", "vsmallECPP.py"))
v1mod = importlib.util.module_from_spec(_sp); _sp.loader.exec_module(_sp and v1mod)

GP = "gp"
SHORT2 = os.path.join(os.path.dirname(__file__), "short2.gp")


def v1_parse(ints):
    """[(A,x,o,p_next,m,facs)] per level under the ORIGINAL (n^2-smooth) parse."""
    p0 = ints[0]; n = p0.bit_length(); n2 = n * n
    primes = sieve_primes(n2)
    out = []
    for j in range(1, len(ints), 3):
        A, x, o = ints[j], ints[j+1], ints[j+2]
        m = 1; facs = []
        oo = o
        for q in primes:
            if q * q > oo and oo <= n2:
                if oo > 1:
                    facs.append(oo); m *= oo; oo = 1
                break
            if q > n2:
                break
            if oo % q == 0:
                facs.append(q)
                while oo % q == 0:
                    oo //= q; m *= q
        out.append((A, x, o, o // m, m, facs))
    return out


def first_bad_level(ints):
    """first level violating the v2 data constraints, or None."""
    from math import ceil
    n = ints[0].bit_length(); lg = log2(n)
    B = ceil(n * n / lg); radlim = ceil(n / lg)
    for lev, (A, x, o, pnext, m, facs) in enumerate(v1_parse(ints)):
        rad = 1
        for q in facs:
            rad *= q
        if (facs and max(facs) > B) or (rad > 1 and rad.bit_length() > radlim) or rad == 1:
            return lev
    return None


def run_repair(p0, ntop, pstart, workers, tlim, logdir, tag):
    os.makedirs(logdir, exist_ok=True)
    procs = []
    for w in range(workers):
        seed = random.randrange(1, 2**31)
        inp = os.path.join(logdir, f"{tag}.w{w}.in")
        out = os.path.join(logdir, f"{tag}.w{w}.out")
        with open(inp, "w") as f:
            f.write(f"\\r {SHORT2}\nsetrand({seed});SC_tlim={tlim};printshort2from({pstart},{ntop});\n")
        procs.append((subprocess.Popen([GP, "-q"], stdin=open(inp), stdout=open(out, "w"),
                                       stderr=open(out + ".err", "w")), out))
    t0 = time.time()
    tail = None
    while tail is None:
        time.sleep(2)
        alive = 0
        for pr, out in procs:
            if pr.poll() is None:
                alive += 1
            try:
                with open(out) as f:
                    for line in f:
                        toks = line.split()
                        if toks and toks[0] == str(pstart) and len(toks) >= 4:
                            tail = [int(t) for t in toks]
                            break
            except FileNotFoundError:
                pass
            if tail:
                break
        if tail is None and alive == 0:
            return None, time.time() - t0
    for pr, _ in procs:                      # kill losers by handle (PID-safe)
        if pr.poll() is None:
            pr.terminate()
    for pr, _ in procs:
        try:
            pr.wait(timeout=10)
        except subprocess.TimeoutExpired:
            pr.kill()
    return tail, time.time() - t0


def main():
    csv_in, outdir = sys.argv[1], sys.argv[2]
    logdir = os.path.join(outdir, "log")
    os.makedirs(outdir, exist_ok=True); os.makedirs(logdir, exist_ok=True)
    random.seed(20260821)
    chains = [[int(t) for t in line.strip().split(',')]
              for line in open(csv_in) if line.strip()]
    results = [None] * len(chains)
    jobs = []
    for idx, ints in enumerate(chains):
        lev = first_bad_level(ints)
        n = ints[0].bit_length()
        if lev is None:
            results[idx] = ints
            print(f"[{n:4d}] pass, kept", flush=True)
            continue
        parsed = v1_parse(ints)
        pstart = ints[0] if lev == 0 else parsed[lev - 1][3]
        prefix = []
        for A, x, o, *_ in parsed[:lev]:
            prefix += [A, x, o]
        jobs.append((pstart.bit_length(), idx, ints[0], n, pstart, prefix, lev))
    jobs.sort(reverse=True)                  # big moduli first
    print(f"{len(jobs)} chains to repair", flush=True)
    for bits, idx, p0, ntop, pstart, prefix, lev in jobs:
        workers = 14 if bits >= 600 else 8 if bits >= 450 else 4 if bits >= 300 else 2 if bits >= 150 else 1
        tlim = 30 if bits >= 450 else 20
        tag = f"n{ntop}.lev{lev}"
        print(f"[{ntop:4d}] repairing from level {lev} (modulus {bits} bits, "
              f"{workers} workers) ...", flush=True)
        tail, dt = run_repair(p0, ntop, pstart, workers, tlim, logdir, tag)
        if tail is None:
            print(f"[{ntop:4d}] FAILED (workers died), see {logdir}/{tag}.*", flush=True)
            continue
        full = [p0] + prefix + tail[1:]
        okv2 = v2_verify(full)
        okv1 = v1mod.verify(full)
        print(f"[{ntop:4d}] repaired in {dt:.0f}s: v2={okv2[0]} v1={okv1} "
              f"({(len(full)-1)//3} levels)", flush=True)
        if okv2[0] and okv1:
            results[idx] = full
        else:
            print(f"[{ntop:4d}] VALIDATION FAILED: {okv2}", flush=True)
    with open(os.path.join(outdir, "certs.csv"), "w") as f:
        for r in results:
            if r is not None:
                f.write(",".join(str(v) for v in r) + "\n")
    done = sum(1 for r in results if r is not None)
    print(f"wrote {done}/{len(chains)} chains to {outdir}/certs.csv", flush=True)


if __name__ == "__main__":
    main()
