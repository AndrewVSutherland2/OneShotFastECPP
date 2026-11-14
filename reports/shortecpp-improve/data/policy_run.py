#!/usr/bin/env python3
"""Run a stopping policy with real gmp-ecm on real candidate orders and count
the winners, tier by tier, for comparison with the model's prediction.

usage: policy_run.py <expdir> [policy.json] [out.json]
<expdir> holds p768.txt, strip.txt ("N S T s prp" per order) and strip.txt.fac
(the factorisation of each s).  policy.json: {"tiers": [[method, B1, curves],..],
"cls": [class of each order], "go": {class: [[[lo, hi), ...] per tier]}} where a
composite tail of t bits (rounded) in that class runs tier j iff t lies in one
of the intervals.
"""
import concurrent.futures as cf
import json
import math
import os
import random
import re
import resource
import subprocess
import sys
import time

E = sys.argv[1]
POL = sys.argv[2] if len(sys.argv) > 2 else os.path.join(E, "policy.json")
OUTF = sys.argv[3] if len(sys.argv) > 3 else os.path.join(E, "policy_run.json")
BENCH = os.path.expanduser("~/ecmbench")
ECM = BENCH + "/inst/deb/usr/bin/ecm"
ENV = dict(os.environ)
ENV["LD_LIBRARY_PATH"] = BENCH + "/inst/deb/usr/lib/x86_64-linux-gnu"
NCPU = os.cpu_count()
RNG = random.Random(7)


def log(*a):
    print(time.strftime("[%H:%M:%S]"), *a, flush=True)


def is_prp(n, rounds=20):
    if n < 2:
        return False
    for b in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % b == 0:
            return n == b
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for _ in range(rounds):
        a = RNG.randrange(2, n - 1)
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


def window_ok(sfac, q, L, p, bcap, radlim, qfloor, cap=500000):
    """Revised-format acceptance: some divisor m > 1 of s (primes <= bcap, log2
    rad(m) <= radlim) puts o = m*q in (L, r(m)*L); q above the recursion floor
    and below sqrt p."""
    if q <= qfloor or q * q >= p:
        return False
    primes = sorted(pr for pr in sfac if pr <= bcap)
    if not primes:
        return False
    mbound = (primes[-1] * L) // q + 1
    divs = [(1, 0, 1)]                      # (m, least prime, radical)
    for pr in primes:
        new = []
        for (m, lp, rad) in divs:
            mm = m
            for _ in range(sfac[pr]):
                mm *= pr
                if mm > mbound:
                    break
                new.append((mm, lp if lp else pr, rad * pr))
            if len(divs) + len(new) > cap:
                break
        divs += new
        if len(divs) > cap:
            break
    for (m, lp, rad) in divs:
        if m > 1 and L < m * q < lp * L and rad.bit_length() <= radlim:
            return True
    return False


_IN = re.compile(r"^Input number is (\d+)")
_FA = re.compile(r"^Found (?:prime|composite) factor of +\d+ digits: (\d+)")
_CO = re.compile(r"^(Prime|Composite) cofactor (\d+) has")


def run_ecm(nums, method, b1, curves):
    """One gmp-ecm process over several numbers -> ({n: (factors, cof, cofprime)}, seconds)."""
    cmd = [ECM] + (["-pm1"] if method == "pm1" else []) + ["-c", str(curves), str(b1)]
    t = time.time()
    r = subprocess.run(cmd, input="\n".join(map(str, nums)) + "\n", capture_output=True, text=True, env=ENV)
    dt = time.time() - t
    res = {}
    cur, facs, cof, cp = None, [], None, False

    def flush():
        if cur is not None and facs and cof is not None and math.prod(facs) * cof == cur:
            res[cur] = (facs, cof, cp)
    for line in r.stdout.splitlines():
        m = _IN.match(line)
        if m:
            flush()
            cur, facs, cof, cp = int(m.group(1)), [], None, False
            continue
        m = _FA.match(line)
        if m:
            facs.append(int(m.group(1)))
            continue
        m = _CO.match(line)
        if m:
            cof, cp = int(m.group(2)), m.group(1) == "Prime"
    flush()
    return res, dt


def main():
    p = int(open(os.path.join(E, "p768.txt")).read())
    n = p.bit_length()
    rt = math.isqrt(p)
    L = rt + 1 + math.isqrt(4 * rt)
    bcap = math.ceil(n * n / math.log2(n))
    radlim = math.ceil(n / math.log2(n))
    qfloor = bcap * bcap
    pol = json.load(open(POL))
    tiers = pol["tiers"]
    cls = pol["cls"]
    go = {int(c): v for c, v in pol["go"].items()}
    sec = pol.get("sec_per_curve", {})

    def admitted(c, j, tb):
        for lo, hi in go.get(c, [[]] * len(tiers))[j]:
            if lo <= tb < hi:
                return True
        return False

    tail, sfac = [], []
    with open(os.path.join(E, "strip.txt")) as f, open(os.path.join(E, "strip.txt.fac")) as g:
        for line, fl in zip(f, g):
            w = line.split()
            tail.append(int(w[2]))
            sfac.append({int(a): int(b) for a, b in (tok.split("^") for tok in fl.split())})
    ncand = len(tail)
    assert len(cls) == ncand, (len(cls), ncand)
    log("%d orders, %d tiers, %d cores" % (ncand, len(tiers), NCPU))
    winners = []
    stats = []

    def settle(i, q, where):
        """Tail of order i is the prime q: a win or a dead end."""
        if window_ok(sfac[i], q, L, p, bcap, radlim, qfloor):
            winners.append({"i": i, "qbits": q.bit_length(), "tier": where, "q": str(q)})
            return 1
        return 0

    # fresh classification: prime tails are settled at once
    alive = []
    w0 = 0
    with open(os.path.join(E, "strip.txt")) as f:
        for i, line in enumerate(f):
            prp = line.split()[4] == "1"
            if tail[i] == 1:
                continue
            if prp and is_prp(tail[i]):
                w0 += settle(i, tail[i], -1)
            else:
                alive.append(i)
    log("strip: %d composite tails, %d winners at the strip" % (len(alive), w0))
    stats.append({"tier": -1, "run": ncand, "wins": w0})
    t00 = time.time()
    for j, (method, b1, curves) in enumerate(tiers):
        todo = [i for i in alive if admitted(cls[i], j, int(round(math.log2(tail[i]))))]
        alive = todo                                   # the rest are abandoned for good
        if not todo:
            stats.append({"tier": j, "method": method, "b1": b1, "curves": curves, "run": 0, "split": 0,
                          "wins": 0, "cpu_s": 0.0, "wall": 0.0})
            continue
        # jobs: chunk cheap tails into one process; split the curves of an expensive
        # sweep over idle cores when few tails are left
        est = float(sec.get(str(b1), 1e-5 * b1)) * curves
        par = 1
        if len(todo) < NCPU and method == "ecm":
            par = max(1, min(curves, NCPU // len(todo)))
        chunk = 1 if par > 1 else max(1, min(32, int(1.5 / max(est, 1e-3)), len(todo) // NCPU))
        jobs = []
        if par > 1:
            base, extra = divmod(curves, par)
            for i in todo:
                for k in range(par):
                    c = base + (1 if k < extra else 0)
                    if c:
                        jobs.append(([i], c))
        else:
            for a in range(0, len(todo), chunk):
                jobs.append((todo[a:a + chunk], curves))
        t0 = time.time()
        cpu = 0.0
        found = {}                                     # i -> set of factors of the tail found this tier
        cprime = {}

        def work(job):
            ids, c = job
            res, dt = run_ecm([tail[i] for i in ids], method, b1, c)
            return ids, res, dt
        with cf.ThreadPoolExecutor(max_workers=NCPU) as ex:
            for ids, res, dt in ex.map(work, jobs):
                cpu += dt
                for i in ids:
                    if tail[i] in res:
                        facs, cof, cp = res[tail[i]]
                        found.setdefault(i, set()).update(facs)
        nsplit = wins = 0
        nxt = []
        for i in todo:
            if i in found:
                nsplit += 1
                T = tail[i]
                for f in sorted(found[i]):
                    g = math.gcd(f, T)
                    while g > 1 and T % g == 0:
                        T //= g
                tail[i] = T
                if T == 1:
                    continue
                if is_prp(T):
                    wins += settle(i, T, j)
                    continue
            nxt.append(i)
        alive = nxt
        st = {"tier": j, "method": method, "b1": b1, "curves": curves, "run": len(todo), "split": nsplit,
              "wins": wins, "cpu_s": cpu, "wall": time.time() - t0, "par": par, "chunk": chunk}
        stats.append(st)
        log("tier %2d %s B1=%d c=%d: run %d, split %d, wins %d, cpu %.0fs, wall %.0fs (par %d, chunk %d); total wins %d"
            % (j, method, b1, curves, len(todo), nsplit, wins, cpu, st["wall"], par, chunk, len(winners)))
        json.dump({"stats": stats, "winners": winners, "ncand": ncand, "elapsed": time.time() - t00},
                  open(OUTF + ".tmp", "w"))
        os.replace(OUTF + ".tmp", OUTF)
    log("RUN_DONE: %d winners, %.0f cpu-s of ECM, %.0fs wall"
        % (len(winners), sum(s.get("cpu_s", 0) for s in stats), time.time() - t00))


if __name__ == "__main__":
    _soft, _hard = resource.getrlimit(resource.RLIMIT_NOFILE)
    resource.setrlimit(resource.RLIMIT_NOFILE, (min(_hard, 65536), _hard))   # 192 workers x several pipes each
    main()
