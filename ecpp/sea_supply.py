#!/usr/bin/env python3
"""sea_supply.py -- SEA-sampled candidate supply for one-shot ECPP certificates.

Alternative to the CM/dscan discriminant scan: sample random Montgomery curves
E_A : y^2 = x^3 + A x^2 + x over F_p, compute the trace t with PARI's ellap
(SEA), and gate the two twist orders N = p+1-/+t through the batched
n^4-smoothness engine (`smoothtest gate`, fed one dummy record "1 t 1" per
curve).  A winning N (n^4-smooth part > L) yields a certificate
(p, A, x0, m, q_1..q_k) in voneshot format: x0 is found by rejection-sampling
a point on the winning twist (the Kummer-line formulas never need the twist
made explicit, so x0 just has to have f(x0) = x0^3+A*x0^2+x0 of the right
quadratic character) and pushing it through the scalar [N/m]; the assembled
tuple is accepted only if voneshot.verify() passes -- the challenge verifier
itself is the acceptance oracle.

usage:
  sea_supply.py (p=<dec> | pbits=<n> [seed=<s>]) [workers=N] [batch=1000]
                [gatesec=300] [cap=0] [maxcurves=0] [collect=0] [out=<file>]

  workers   parallel gp point-counting processes (default: os.cpu_count()-2)
  batch     run the smooth gate every `batch` new curves (default 1000)
  gatesec   ... or every `gatesec` seconds, whichever comes first (default 300)
  cap       per-curve ellap time cap in seconds (0 = none): curves whose point
            count exceeds the cap are abandoned (cost-based early abort)
  maxcurves stop after this many curves even without a winner (0 = no limit)
  collect   1 = do not stop at the first winner (density measurement mode)
  out       certificate output path (default certs/sea/oneshot_sea_<n>.txt)

Every curve sampled is appended to <out>.cands ("A t" lines, header "# p=..")
so an interrupted run loses nothing; a rerun with the same p resumes from it.
"""

import os
import re
import subprocess
import sys
import threading
import time
import random as pyrandom
from math import gcd
from queue import Queue, Empty

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SMOOTHTEST = os.path.join(HERE, "smoothtest")
PCACHE = os.environ.get("ONESHOT_PCACHE_DIR", os.path.join(ROOT, "work", "pcache"))

sys.path.insert(0, ROOT)
try:
    import voneshot  # the challenge verifier: ladder + verify (acceptance oracle)
except ModuleNotFoundError:      # fresh checkout: use the bundled copy
    sys.path.insert(0, os.path.join(ROOT, "verifier-batching-pr"))
    import voneshot


def log(msg):
    print("[%8.1fs] %s" % (time.time() - T0, msg), flush=True)


T0 = time.time()


# ---------------------------------------------------------------- p and cache
def canonical_p(pbits, seed):
    """The project-wide test prime for (pbits, seed): GMP urandomb with top bit
    set, then nextprime -- obtained from smoothtest so it matches dscan/oneshot
    exactly (GMP's generator is not reproducible from Python)."""
    r = subprocess.run([SMOOTHTEST, "gate", "pbits=%d" % pbits, "seed=%d" % seed,
                        "y=65536"], input="", capture_output=True, text=True)
    m = re.search(r"^p = (\d+)$", r.stderr, re.M)
    if not m:
        sys.exit("could not obtain p from smoothtest:\n" + r.stderr)
    return int(m.group(1))


def ensure_pcache(n4, threads):
    path = os.path.join(PCACHE, "oneshot_P_%d.bin" % n4)
    if os.path.exists(path):
        return path
    os.makedirs(PCACHE, exist_ok=True)
    log("building prime-product cache for y=%d (one-time, minutes)..." % n4)
    tmp = path + ".tmp"
    r = subprocess.run([SMOOTHTEST, "pbuild", "y=%d" % n4, "threads=%d" % threads,
                        "save=%s" % tmp], capture_output=True, text=True)
    if r.returncode != 0 or not os.path.exists(tmp):
        sys.exit("pbuild failed:\n" + r.stderr)
    os.replace(tmp, path)                 # atomic: no partial cache on interrupt
    log("cache built: %s" % path)
    return path


# ---------------------------------------------------------------- gp workers
GP_SCRIPT = r"""default(parisizemax,"1G");
p = {p}; setrand({seed});
while(1, my(A = 2 + random(p-4)); if (A == 2 || A == p-2, next); my(E = ellinit([0,A,0,1,0], p)); my(t = {count}); if (type(t) == "t_INT", print(A, " ", t)));
"""


def start_worker(i, p, seed, cap, queue, procs):
    count = "ellap(E)" if not cap else 'iferr(alarm(%d, ellap(E)), err, "X")' % cap
    script = GP_SCRIPT.format(p=p, seed=seed * 100003 + i, count=count)
    proc = subprocess.Popen(["stdbuf", "-oL", "gp", "-q"],
                            stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                            stderr=subprocess.DEVNULL, text=True)
    proc.stdin.write(script)
    proc.stdin.flush()
    procs.append(proc)

    def reader():
        for line in proc.stdout:
            parts = line.split()
            if len(parts) == 2 and parts[0].isdigit():
                try:
                    queue.put((int(parts[0]), int(parts[1])))
                except ValueError:
                    pass
    threading.Thread(target=reader, daemon=True).start()


# ---------------------------------------------------------------- smooth gate
WIN_RE = re.compile(r"^WIN D=-\d+ t=(-?\d+) order=(\d+)\s+m=(\d+)\s+smoothpart=\d+\s+q=\[([\d,]*)\]\s+OK$")


def run_gate(p, pfile, threads, traces):
    """Feed |t| values through `smoothtest gate`; return [(t_abs, N, m, qs)]."""
    inp = "".join("1 %d 1\n" % t for t in traces)
    r = subprocess.run([SMOOTHTEST, "gate", "p=%d" % p, "load=%s" % pfile,
                        "threads=%d" % threads], input=inp,
                       capture_output=True, text=True)
    if r.returncode != 0:
        log("gate FAILED (rc=%d):\n%s" % (r.returncode, r.stderr[-2000:]))
        return None                       # invocation failure: caller retries the batch
    wins = []
    for line in r.stdout.splitlines():
        mm = WIN_RE.match(line)
        if mm:
            t, N, m = int(mm.group(1)), int(mm.group(2)), int(mm.group(3))
            qs = tuple(int(q) for q in mm.group(4).split(",")) if mm.group(4) else ()
            wins.append((abs(t), N, m, qs))
    tm = re.search(r"smooth_parts: ([0-9.]+)s", r.stderr)
    log("gate: %d traces, %d winner(s)%s" %
        (len(traces), len(wins), (", %ss smooth_parts" % tm.group(1)) if tm else ""))
    return wins


# ---------------------------------------------------------------- assembly
def assemble(p, A, N, m, qs, tries=256):
    """Find x0 with voneshot.verify(p, A, x0', m, qs) True, where x0' is the
    x-coordinate of Q = [N/m]P for a sampled P on the twist of order N.
    chi = f(x0)^((p-1)/2): +1 -> P on E_A (order p+1-t), -1 -> quadratic twist
    (order p+1+t).  Rather than track which sign the gate's N corresponds to,
    try chi=+1 points first and fall back to chi=-1: points on the wrong twist
    simply fail verify() (their order does not divide N)."""
    e2 = (p - 1) // 2
    cof = N // m
    rng = pyrandom.Random(0xC0FFEE ^ A ^ N)
    fails = 0
    for want in (1, p - 1):
        for _ in range(tries // 2):
            x0 = rng.randrange(2, p)
            f = (x0 * x0 % p * x0 + A * x0 * x0 + x0) % p
            if f == 0:
                continue
            if pow(f, e2, p) != want:
                continue
            XQ, ZQ = voneshot.ladder(cof, x0, 1, A, p)
            if ZQ % p == 0:
                continue
            x0q = XQ * pow(ZQ, -1, p) % p
            if voneshot.verify(p, A, x0q, m, qs):
                return x0q, fails
            fails += 1
    return None, fails


# ---------------------------------------------------------------- main
def main():
    args = dict(a.split("=", 1) for a in sys.argv[1:])
    pbits = int(args.get("pbits", 0))
    seed = int(args.get("seed", 1))
    workers = int(args.get("workers", max(2, (os.cpu_count() or 4) - 2)))
    batch = int(args.get("batch", 1000))
    gatesec = float(args.get("gatesec", 300))
    cap = int(args.get("cap", 0))
    maxcurves = int(args.get("maxcurves", 0))
    collect = int(args.get("collect", 0))

    if "p" in args:
        p = int(args["p"])
    elif pbits:
        p = canonical_p(pbits, seed)
    else:
        sys.exit(__doc__)
    n = p.bit_length()
    n4 = (n * n) ** 2
    out = args.get("out", os.path.join(ROOT, "certs", "sea", "oneshot_sea_%d.txt" % n))
    outdir = os.path.dirname(out)
    if outdir:
        os.makedirs(outdir, exist_ok=True)
    cands_path = out + ".cands"

    log("p = %d (%d bits), n^4 = %d, workers=%d cap=%s" % (p, n, n4, workers, cap or "off"))
    pfile = ensure_pcache(n4, workers)

    # resume: reload previously sampled curves for this p
    seen = {}          # A -> t
    header = "# p=%d\n" % p
    if os.path.exists(cands_path):
        with open(cands_path) as f:
            first = f.readline()
            if first == header:
                for line in f:
                    a, t = line.split()
                    seen[int(a)] = int(t)
                log("resumed %d curves from %s" % (len(seen), cands_path))
            else:
                os.rename(cands_path, cands_path + ".old")
    cands_f = open(cands_path, "a" if seen else "w")
    if not seen:
        cands_f.write(header)

    queue, procs = Queue(), []
    for i in range(workers):
        start_worker(i, p, seed, cap, queue, procs)
    log("%d gp workers started" % workers)

    gated = set()      # |t| values already through the gate
    pending = dict(seen)  # A -> t not yet gated
    tmap = {}          # |t| -> [A, ...]
    for a, t in seen.items():
        tmap.setdefault(abs(t), []).append(a)
    ncurves = len(seen)
    winners_done = set()
    certs = 0
    last_gate = time.time()
    win_curve_count = 0

    def gate_and_assemble():
        nonlocal certs, win_curve_count
        new = [t for t in {abs(t) for t in pending.values()} if t not in gated]
        if not new:
            pending.clear()
            return False
        batch = dict(pending)
        pending.clear()
        wins = run_gate(p, pfile, workers, sorted(new))
        if wins is None:                  # gate invocation failed: keep the
            pending.update(batch)         # traces so a later pass retries them
            return False
        gated.update(new)
        got = False
        for t_abs, N, m, qs in wins:
            key = (N, m)
            if key in winners_done:
                continue
            winners_done.add(key)
            for A in tmap.get(t_abs, []):
                x0, fails = assemble(p, A, N, m, qs)
                if x0 is None:
                    log("winner t=%d A=%d: assembly failed (%d tries) -- skipped" % (t_abs, A, fails))
                    continue
                cert = [p, A, x0, m] + list(qs)
                with open(out, "w") as f:
                    f.write(" ".join(map(str, cert)) + "\n")
                certs += 1
                got = True
                win_curve_count = win_curve_count or ncurves
                log("WINNER after %d curves (%d candidate orders): |D|-free, t=%d" % (ncurves, 2 * ncurves, t_abs))
                log("  N = %d" % N)
                log("  m = %d (%d bits), %d large primes %s" % (m, m.bit_length(), len(qs), list(qs)))
                log("  cert -> %s" % out)
                break
        return got

    try:
        done = False
        while not done:
            try:
                A, t = queue.get(timeout=1.0)
                if A not in seen:
                    seen[A] = t
                    pending[A] = t
                    tmap.setdefault(abs(t), []).append(A)
                    ncurves += 1
                    cands_f.write("%d %d\n" % (A, t))
                    if ncurves % 200 == 0:
                        cands_f.flush()
                        rate = ncurves / (time.time() - T0)
                        log("%d curves sampled (%.2f curves/s, %.1f s/curve/worker)" %
                            (ncurves, rate, workers / rate if rate else 0))
            except Empty:
                pass
            if len(pending) >= batch or (time.time() - last_gate > gatesec and pending):
                last_gate = time.time()
                if gate_and_assemble() and not collect:
                    done = True
            if maxcurves and ncurves >= maxcurves:
                log("maxcurves reached")
                break
            if not any(pr.poll() is None for pr in procs):
                log("all workers died")
                break
    finally:
        for pr in procs:
            pr.kill()
        cands_f.close()

    if pending:
        gate_and_assemble()
    wall = time.time() - T0
    log("done: %d curves, %d cert(s), %.1f s wall (~%.0f core-s of SEA)" %
        (ncurves, certs, wall, wall * workers))
    if certs and win_curve_count:
        log("observed winner density: 1 per %d curves (%d candidate orders)" %
            (win_curve_count, 2 * win_curve_count))
    sys.exit(0 if certs else 1)


if __name__ == "__main__":
    main()
