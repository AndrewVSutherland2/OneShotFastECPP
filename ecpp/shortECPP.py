#!/usr/bin/env python3
"""shortECPP.py -- CM-method prover for short ECPP certificates
(github.com/AndrewVSutherland/ShortPrimalityProofs format).

At each level (modulus p, top-level n fixed for the whole chain, n2 = n^2):
  1. dscan enumerates (D,t,v) with 4p = t^2 + |D| v^2; candidate curve orders
     are N = p+1-t and p+1+t with 4 | N (Montgomery-representable side pairs).
  2. smoothtest strips every candidate's y0-smooth part S (y0 = 2^26, batched
     Bernstein remainder trees against a cached prime product); the n2-smooth
     part s | S supplies the certificate cofactor m, the rest of S is discarded
     into the untracked cofactor c (which never needs factoring), and the rough
     tail T = N/S is peeled by breadth-first rounds of gmp-ecm with escalating
     B1 and per-round admission caps (rank: smallest tail first).
  3. A candidate wins when a known prime q (a peeled factor, or a tail that
     became prime) satisfies n2 < q, q^2 < p, and some divisor m | s lands
     o = m*q in the window L < o < r(m)*L, L = (p^(1/4)+1)^2.  The winner's
     curve comes from H_D (cm_method / classpoly, PARI polclass fallback), the
     Montgomery A from the u-cubic, and the point of exact order o from the
     x-only ladder; the certificate stores (A, x, o) and descends to q.
Levels below CM_MIN_BITS are finished by short.gp (SEA is cheap there), called
with the top-level n2.  The chain is verified by vshort2.verify (v2=1, the
revised format; vsmallECPP.verify for the original) and independently
cross-checked level by level in PARI before it is written.

usage: shortECPP.py (p=<decimal> | c=<digits: p=nextprime(10^c)>)
                      [threads=N] [out=<file>] [tag=<name>] [B0=<dscan start>]
                      [Bmax=<dscan cap>] [maxfb=<factor-base cap>] [seed=1]
                      [v=1] [v2=1] [resume=1]

The scan's factor base is always capped: at maxfb if given (top level only),
otherwise at fb_cap(B), which grows with the scan bound.  Its square roots are
cached per modulus under work/rootcache/ for the life of the level.
"""

import concurrent.futures
import json
import math
import os
import random
import hashlib
import re
import shutil
import subprocess
import sys
import threading
import time
from math import gcd, isqrt

ECPP_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_DIR = os.path.dirname(ECPP_DIR)
SPP_DIR = os.environ.get("SPP_DIR", "")
if not SPP_DIR:
    for cand_dir in ("/home/claude/ShortPrimalityProofs",
                     os.path.expanduser("~/ShortPrimalityProofs")):
        if any(os.path.exists(os.path.join(cand_dir, v))
               for v in ("vshortECPP.py", "vsmallECPP.py")):
            SPP_DIR = cand_dir
            break
def _load_verifier():
    """Prefer the ShortPrimalityProofs checkout's verifier, but only when it
    exports every symbol we need (older copies predate the batched helpers);
    otherwise use the known-compatible bundled module."""
    import importlib.util
    names = ("ladder", "verify", "sieve_primes", "balanced_product", "remainder_tree")
    paths = []
    if SPP_DIR:
        paths.append(os.path.join(SPP_DIR, "vshortECPP.py"))
        paths.append(os.path.join(SPP_DIR, "vsmallECPP.py"))
    paths.append(os.path.join(ECPP_DIR, "vsmallECPP.py"))
    for path in paths:
        if not os.path.exists(path):
            continue
        spec = importlib.util.spec_from_file_location("vsmallECPP", path)
        mod = importlib.util.module_from_spec(spec)
        try:
            spec.loader.exec_module(mod)
        except Exception:
            continue
        if all(hasattr(mod, n) for n in names):
            sys.modules["vsmallECPP"] = mod
            return mod
    sys.exit("no compatible vsmallECPP.py found (need: ladder, verify, "
             "sieve_primes, balanced_product, remainder_tree)")

_verifier = _load_verifier()
ladder = _verifier.ladder
verify = _verifier.verify
sieve_primes = _verifier.sieve_primes
balanced_product = _verifier.balanced_product
remainder_tree = _verifier.remainder_tree

ECM_BIN = os.environ.get("CHAIN_ECM", "")
if not ECM_BIN:
    for cand_bin in ("/home/claude/.local/ecmpkg/usr/bin/ecm",
                     os.path.expanduser("~/.local/ecmpkg/usr/bin/ecm")):
        if os.path.exists(cand_bin):
            ECM_BIN = cand_bin
            break
    else:
        ECM_BIN = shutil.which("ecm") or ""
ECM_ENV = dict(os.environ)
_ecm_lib = os.path.join(os.path.dirname(os.path.dirname(ECM_BIN)), "lib", "x86_64-linux-gnu") if ECM_BIN else ""
ECM_ENV["LD_LIBRARY_PATH"] = _ecm_lib + ":" + ECM_ENV.get("LD_LIBRARY_PATH", "")
# Strip bound y0 = 2^26.  The strip only has to be exact up to n^2 (the
# certificate cofactor m is n^2-smooth); everything between n^2 and y0 is junk
# that ECM would remove anyway.  The batch cost is one reduction of the prime
# product modulo the product of the batch, so it scales with y0: 2^26 strips
# 27,000 orders at 1030 bits ten times faster than 2^32 did (6 s against 62 s),
# with identical n^2-smooth parts, and the primes it leaves, 27 to 32 bits, are
# found by the first ECM tier.
Y0 = 1 << 26
PCACHE = os.path.join(REPO_DIR, "work/pcache/oneshot_P_%d.bin" % Y0)
CM_MIN_BITS = 135                 # below this, short.gp (SEA) finishes the chain

# Revised (August 2026) format support: when V2_BCAP > 0, fillers are restricted
# to primes <= V2_BCAP = ceil(n^2/log2 n) with log2 rad(m) <= V2_RADLIM =
# ceil(n/log2 n), and recursion requires q > V2_QFLOOR = V2_BCAP^2 (terminal
# primes below B^2 come from the gp tail).  Set by prove_chain(v2=True); the
# final certificate is then gated on vshort2.verify.  Zero = original format.
V2_BCAP = 0
V2_RADLIM = 0
V2_QFLOOR = 0

# classpoly environment (mirrors setenv.sh); cm_method shells out to the
# classpoly binary, so it must be on PATH in the subprocess environment
ECM_ENV["CLASSPOLY_PHI_DIR"] = ECM_ENV.get("CLASSPOLY_PHI_DIR", os.path.join(REPO_DIR, "phi_files"))
ECM_ENV["CLASSPOLY_H_DIR"] = ECM_ENV.get("CLASSPOLY_H_DIR", os.path.join(REPO_DIR, "work/H_files"))
ECM_ENV["ONESHOT_PCACHE_DIR"] = ECM_ENV.get("ONESHOT_PCACHE_DIR", os.path.join(REPO_DIR, "work/pcache"))
ECM_ENV["PATH"] = (os.path.join(REPO_DIR, "classpoly_v1.0.3") + ":" + ECPP_DIR + ":"
                   + ECM_ENV.get("PATH", ""))

VERBOSE = 1
CAMP_DIR = os.path.join(REPO_DIR, "work", "shortcamp")
JOURNAL = None       # set per run: winners survive crashes / spot preemptions


def log(msg):
    if VERBOSE:
        print("[%s] %s" % (time.strftime("%H:%M:%S"), msg), file=sys.stderr, flush=True)


def journal_write(rec):
    if JOURNAL:
        with open(JOURNAL, "a") as f:
            f.write(json.dumps(rec) + "\n")


def journal_read(path):
    out = []
    if path and os.path.exists(path):
        for line in open(path):
            line = line.strip()
            if line:
                try:
                    out.append(json.loads(line))
                except ValueError:
                    pass
    return out


# ------------------------------------------------------------------ utilities
def is_prp(n, rounds=24):
    if n < 2:
        return False
    for b in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % b == 0:
            return n == b
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
    if n >= 1 << 64:
        rng = random.Random(n & ((1 << 64) - 1))
        bases += [rng.randrange(2, n - 2) for _ in range(max(0, rounds - 12))]
    for b in bases:
        x = pow(b, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def sqrt_mod(a, p):
    a %= p
    if a == 0:
        return 0
    if pow(a, (p - 1) // 2, p) != 1:
        return None
    if p % 4 == 3:
        return pow(a, (p + 1) // 4, p)
    s, q = 0, p - 1
    while q % 2 == 0:
        q //= 2
        s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1
    m, c, t, r = s, pow(z, q, p), pow(a, q, p), pow(a, (q + 1) // 2, p)
    while t != 1:
        i, tt = 0, t
        while tt != 1:
            tt = tt * tt % p
            i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c = i, b * b % p
        t, r = t * c % p, r * b % p
    return r


def gp(script, timeout=1200):
    inp = 'default(parisizemax,"12G");\n' + script
    r = subprocess.run(["gp", "-q"], input=inp, capture_output=True, text=True, timeout=timeout)
    if r.returncode != 0:
        raise RuntimeError("gp failed: %s" % r.stderr[:800])
    return r.stdout


def factor_small(s, small_primes):
    """Factor an n2-smooth integer by trial division over the sieved primes."""
    fs = {}
    for q in small_primes:
        if q * q > s:
            break
        while s % q == 0:
            fs[q] = fs.get(q, 0) + 1
            s //= q
    if s > 1:
        fs[s] = fs.get(s, 0) + 1
    return fs


# ------------------------------------------------------------ candidate state
class Cand:
    __slots__ = ("d", "t", "v", "N", "s", "sfac", "mid", "tail", "dead",
                 "splits", "effort", "prp_checked", "tcurves", "last_b1")

    def __init__(self, d, t, v, N):
        self.d, self.t, self.v, self.N = d, t, v, N
        self.s = 1          # n2-smooth part of N (certificate m divides this)
        self.sfac = None    # factorization of s (computed lazily)
        self.mid = 1        # (y0-smooth part)/s: discarded into the cofactor
        self.tail = N       # unfactored rough part
        self.dead = False
        self.splits = 0
        self.effort = 0.0
        self.prp_checked = False
        self.tcurves = 0    # failed curves at the current B1 tier (Bayes discount)
        self.last_b1 = 0


# ------------------------------------------------------------ pool components
def fb_cap(B):
    """Factor-base cap for a scan to B.  A square root mod p costs about 60
    Cornacchia tests and a prime q is in about B/(30 q) of the discriminants
    scanned, so roots of primes far above B/1000 never pay for themselves; the
    count of solvable D is flat in the cap (measured at 1030 bits, B = 2e10:
    a cap of B/4000 keeps 81% of what B/400 finds, with a ninth of the roots)."""
    return max(10 ** 6, min(10 ** 8, B // 1000))


def roots_path(p):
    """Per-modulus cache of the scan's factor-base square roots (dscan roots=)."""
    d = os.path.join(REPO_DIR, "work", "rootcache")
    os.makedirs(d, exist_ok=True)
    return os.path.join(d, "roots_%s.bin" % hashlib.sha256(str(p).encode()).hexdigest()[:20])


def run_dscan(p, B, Bmin, threads, maxfb=0, fbmin=0, roots=None):
    """Solvable (d, t, v) with Bmin <= d < B over the factor base capped at maxfb.
    With fbmin set (the cap of the scans already done below Bmin) the d < Bmin
    that need a prime in (fbmin, maxfb] are included too, so successive calls
    with a growing cap tile [0, B) at the latest cap without repeats."""
    cmd = [os.path.join(ECPP_DIR, "dscan"), "p=%d" % p, "B=%d" % B, "threads=%d" % threads, "dump"]
    if Bmin:
        cmd.append("Bmin=%d" % Bmin)
    if maxfb:
        cmd.append("maxfb=%d" % maxfb)
    if fbmin and Bmin:
        cmd.append("fbmin=%d" % fbmin)
    if roots:
        cmd.append("roots=" + roots)
    r = subprocess.run(cmd, capture_output=True, text=True, env=ECM_ENV)
    if r.returncode != 0:
        raise RuntimeError("dscan failed: %s" % r.stderr[:500])
    out = []
    for line in r.stdout.splitlines():
        w = line.split()
        if len(w) == 3:
            d, t, v = int(w[0]), int(w[1]), int(w[2])
            if d > 4:
                out.append((d, t, v))
    out.sort()
    return out


def ensure_pcache(threads):
    """Build the primorial cache on first use (about a second for y0 = 2^26);
    smoothtest pbuild writes it once and load= reuses it."""
    if os.path.exists(PCACHE):
        return
    os.makedirs(os.path.dirname(PCACHE), exist_ok=True)
    log("building primorial cache %s (one-time, y=%d)" % (PCACHE, Y0))
    tmp = "%s.tmp.%d" % (PCACHE, os.getpid())   # per-process: concurrent cold
    r = subprocess.run([os.path.join(ECPP_DIR, "smoothtest"), "pbuild", "y=%d" % Y0,
                        "save=" + tmp, "threads=%d" % threads],
                       capture_output=True, text=True, env=ECM_ENV)
    if r.returncode != 0 or not os.path.exists(tmp):
        if os.path.exists(tmp):
            os.remove(tmp)
        raise RuntimeError("smoothtest pbuild failed: %s" % r.stderr[:500])
    os.replace(tmp, PCACHE)   # atomic; racing builders publish identical content


def smooth_strip(Ns, threads, n2):
    """Batched strip: N -> (S, T, s, prp) with N = S*T, S the y0-smooth part
    (cached prime product, smoothtest parts), s the n2-smooth part of S with
    multiplicity, and prp the base-2 Fermat verdict on T.  The last two are
    computed inside smoothtest, in parallel, so nothing per-candidate is left
    for the interpreter."""
    if not Ns:
        return {}
    ensure_pcache(threads)
    inp = "\n".join(str(N) for N in Ns) + "\n"
    r = subprocess.run([os.path.join(ECPP_DIR, "smoothtest"), "parts", "load=" + PCACHE,
                        "threads=%d" % threads, "n2=%d" % n2, "prp=1"],
                       input=inp, capture_output=True, text=True, env=ECM_ENV)
    if r.returncode != 0:
        raise RuntimeError("smoothtest failed: %s" % r.stderr[:500])
    res = {}
    for line in r.stdout.splitlines():
        w = line.split()
        if len(w) == 5:
            N, S, T, s2 = int(w[0]), int(w[1]), int(w[2]), int(w[3])
            if S * T != N or S % s2 != 0:
                raise RuntimeError("smoothtest self-inconsistency")
            res[N] = (S, T, s2, w[4] == "1")
    for N in Ns:
        if N not in res:
            raise RuntimeError("smoothtest missed an input")
    return res


_ECM_INPUT = re.compile(r"^Input number is (\d+)")
_ECM_FACTOR = re.compile(r"^Found (?:prime|composite) factor of +\d+ digits: (\d+)")
_ECM_COF = re.compile(r"^(Prime|Composite) cofactor (\d+) has")


def parse_ecm(out, nums):
    """gmp-ecm's normal (not -q) report for a list of inputs -> {n: (factors,
    cofactor, cofactor_is_prp)} for the inputs it split.  The cofactor verdict is
    gmp-ecm's own probable-prime test; a split is kept only if the factors and
    the cofactor multiply back to the input."""
    want = set(nums)
    res = {}
    cur, facs, cof, cprp = None, [], None, False

    def flush():
        if cur in want and facs and cof is not None and math.prod(facs) * cof == cur:
            res[cur] = (facs, cof, cprp)
    for line in out.splitlines():
        m = _ECM_INPUT.match(line)
        if m:
            flush()
            cur, facs, cof, cprp = int(m.group(1)), [], None, False
            continue
        m = _ECM_FACTOR.match(line)
        if m:
            facs.append(int(m.group(1)))
            continue
        m = _ECM_COF.match(line)
        if m:
            cof, cprp = int(m.group(2)), m.group(1) == "Prime"
    flush()
    return res


class EcmRunner:
    """Runs gmp-ecm processes for one round and can stop them all at once.

    Each call handles a list of numbers in one process (gmp-ecm takes one input
    per line), which keeps the launch rate down at the shallow tiers where a
    tail needs a fraction of a second.  abort() kills every running process and
    makes pending calls return None, so a round ends as soon as a winner is in
    hand instead of when its slowest tail finishes."""

    def __init__(self):
        self.lock = threading.Lock()
        self.procs = set()
        self.stopped = False

    def run(self, nums, b1, curves, timeout, method="ecm"):
        if self.stopped:
            return None
        cmd = [ECM_BIN] + (["-pm1"] if method == "pm1" else []) + ["-c", str(curves), str(b1)]
        p = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                             stderr=subprocess.DEVNULL, text=True, env=ECM_ENV)
        with self.lock:
            if self.stopped:
                p.kill()
            else:
                self.procs.add(p)
        try:
            out, _ = p.communicate("\n".join(str(n) for n in nums) + "\n", timeout=timeout)
        except subprocess.TimeoutExpired:
            p.kill()
            out, _ = p.communicate()          # keep whatever it had finished
        finally:
            with self.lock:
                self.procs.discard(p)
        if self.stopped:
            return None
        return parse_ecm(out, nums)

    def abort(self):
        with self.lock:
            self.stopped = True
            for p in self.procs:
                p.kill()


def ecm_peel(n, b1, curves, timeout, method="ecm"):
    """One gmp-ecm run (ECM, or P-1 with method="pm1") on one number; returns
    (found_factors, cofactor)."""
    res = EcmRunner().run([n], b1, curves, timeout, method) or {}
    if n in res:
        return res[n][0], res[n][1]
    return [], n


def ecm_chunk(b1, curves, nbits, method, ntails, threads):
    """Tails per gmp-ecm process: enough that a process runs for about 1.5 s
    (a curve costs about 1e-5*B1 s at 1030 bits, P-1 an eighth of that), but
    never so many that fewer than `threads` processes are running."""
    per_tail = 1e-5 * b1 * curves * (nbits / 1030.0) ** 2 / (8 if method == "pm1" else 1)
    k = max(1, min(32, int(1.5 / max(per_tail, 1e-4))))
    return max(1, min(k, ntails // max(1, threads)))


# ------------------------------------------------------------ window matching
def window_hits(sfac, q, L, cap=500000):
    """All (o, m) with m | s (from factorization sfac), m > 1, L < o=m*q < r(m)*L,
    sorted by o.  r(m) = least prime divisor of m; since divisors are built over
    ascending primes, the least prime is the first one used.  Divisors above
    n2*L/q can never satisfy o < r*L (r <= least prime <= n2), so bound there."""
    if V2_BCAP:
        sfac = {pr: e for pr, e in sfac.items() if pr <= V2_BCAP}
    primes = sorted(sfac)
    if not primes:
        return []
    mbound = (primes[-1] * L) // q + 1          # r(m) <= largest prime of s
    divs = [(1, 0)]                             # (m, least prime; 0 = none yet)
    for pr in primes:
        new = []
        for (m, lp) in divs:
            mm = m
            for _ in range(sfac[pr]):
                mm *= pr
                if mm > mbound:
                    break
                new.append((mm, lp if lp else pr))
                if len(divs) + len(new) > cap:
                    break
        divs += new
        if len(divs) > cap:
            break
    hits = []
    for (m, lp) in divs:
        if m > 1:
            o = m * q
            if L < o < lp * L:
                if V2_RADLIM:
                    rad = 1
                    for pr in primes:
                        if m % pr == 0:
                            rad *= pr
                    if rad.bit_length() > V2_RADLIM:    # log2 rad(m) <= ceil(n/log2 n)
                        continue
                hits.append((o, m))
    hits.sort()
    return hits


def accept_q(cand, q, p, L, n2, small_primes):
    """If prime q makes a valid descent from p with cofactor m | s, return the
    sorted (o, m) options, else None."""
    if q <= (V2_QFLOOR if V2_QFLOOR else n2) or q * q >= p:
        return None                       # v2: recursion floor is B^2, not n^2
    if cand.sfac is None:
        cand.sfac = factor_small(cand.s, small_primes)
    if not cand.sfac:
        return None
    hits = window_hits(cand.sfac, q, L)
    # a filler prime ell > 97 that also divides v is unusable: exact-order
    # sampling needs the ell-volcano floor and cm_method walks only ell <= 97,
    # so such hits would burn an expensive j computation and then fail
    hits = [(o, m) for o, m in hits
            if not any(m % pr == 0 and cand.v % pr == 0
                       for pr in cand.sfac if pr > 97)]
    return hits or None


# ------------------------------------------------------------ winner -> level
def get_j(p, d, ells=None, jobs=0, attempts=3):
    """j-invariant of a curve with CM by -d over F_p (cm_method, polclass
    fallback).  The root-find's equal-degree splitting is randomized and an
    attempt can come back rootless (exit 0, no j line), so retry a few times
    before giving up on the discriminant."""
    cmd = [os.path.join(ECPP_DIR, "cm_method"), "D=-%d" % d, "p=%d" % p]
    if ells:
        cmd.append("ells=" + ",".join(str(e) for e in ells))
    if jobs > 1 and d >= 20000000:
        cmd.append("jobs=%d" % min(jobs, 64))
    for att in range(attempts):
        # distinct seed per attempt: cm_method's root-find is deterministically
        # seeded, so an identical rerun would repeat a rootless result exactly
        r = subprocess.run(cmd + ["seed=%d" % (att + 1)],
                           capture_output=True, text=True, env=ECM_ENV, timeout=3600)
        if r.returncode == 0:
            for line in r.stdout.splitlines():
                w = line.split()
                if len(w) == 2 and w[0] == "j":
                    return int(w[1])
        log("    cm_method D=-%d attempt %d/%d: exit %d, no j; tail: %s" % (
            d, att + 1, attempts, r.returncode,
            (r.stderr or r.stdout)[-260:].replace("\n", " | ")))
    if ells is None and d <= 3 * 10 ** 6:
        out = gp("r=polrootsmod(polclass(-%d),%d);\n"
                 "if(#r==0,print(\"NOROOT\"),print(lift(r[1])));\n" % (d, p), timeout=1800)
        w = out.split()
        if w and w[-1] != "NOROOT":
            return int(w[-1])
    return None


def mont_As(j, p):
    """Montgomery coefficients A (nonsingular) for j-invariant j: QR roots u of
    u^3 - 9u^2 + (6912-j)/256 u + (4j-6912)/256, A = sqrt(u)."""
    if j % p == 0 or j % p == 1728 % p:
        return []
    inv256 = pow(256, -1, p)
    c1 = (6912 - j) * inv256 % p
    c0 = (4 * j - 6912) * inv256 % p
    out = gp("v=polrootsmod(x^3-9*x^2+%d*x+%d,%d); print(apply(lift,Vec(v)))\n" % (c1, c0, p))
    line = out.strip().splitlines()[-1]
    if not line.startswith("["):
        return []
    roots = [int(t) for t in line.strip("[]").split(",") if t.strip()]
    As = []
    for u in roots:
        if u % p == 4:
            continue
        A = sqrt_mod(u, p)
        if A is None:
            continue
        if A in (2, p - 2):
            continue
        As.append(A)
    return As


def curve_side(p, A, N, Np, rng):
    """Is #E_A (the B=1 side) equal to N (+1) or N' (-1)?  0 if undetermined."""
    for _ in range(48):
        x = rng.randrange(2, p)
        f = (x * x % p * x + A * x % p * x + x) % p
        if f == 0 or pow(f, (p - 1) // 2, p) != 1:
            continue
        Xa, Za = ladder(N, x, 1, A, p)
        Xb, Zb = ladder(Np, x, 1, A, p)
        a, b = Za % p == 0, Zb % p == 0
        if a and not b:
            return 1
        if b and not a:
            return -1
    return 0


def find_point(p, A, N, o, odiv, chi_target, rng, tries=80):
    """x-coordinate of a point of exact order o on the chi_target twist of E_A
    (certificate coordinate; the twist's B is implicit).  None if not found."""
    c = N // o
    for _ in range(tries):
        x = rng.randrange(2, p)
        f = (x * x % p * x + A * x % p * x + x) % p
        if f == 0:
            continue
        chi = 1 if pow(f, (p - 1) // 2, p) == 1 else -1
        if chi != chi_target:
            continue
        X0, Z0 = ladder(c, x, 1, A, p)
        if Z0 % p == 0:
            continue
        Xo, Zo = ladder(o, X0, Z0, A, p)
        if Zo % p != 0 or gcd(Xo % p, p) != 1:
            continue
        ok = True
        for qd in odiv:
            Xq, Zq = ladder(o // qd, X0, Z0, A, p)
            if Zq % p == 0:
                ok = False
                break
        if not ok:
            continue
        xq = X0 * pow(Z0, -1, p) % p
        return xq
    return None


def build_level(p, cand, q, hits, rng, jobs=0):
    """Turn a winning (candidate, q) into a certificate level (A, x, o).
    For p = 3 mod 4 with N = 4 mod 8, a Montgomery model needs a cyclic
    2-Sylow; when v is even the H_D surface roots generically fail, so descend
    the 2-volcano to its floor first (cm_method ells=2)."""
    # primes shared between any hit's filler and v: their Sylow may sit split
    # above the ell-volcano floor, and exact-order cofactor sampling needs it
    # cyclic, so any descent must cover ALL of them
    need = sorted({pr for _, m_ in hits[:12]
                   for pr in cand.sfac
                   if m_ % pr == 0 and cand.v % pr == 0 and pr <= 97})
    volcano_first = (p % 4 == 3 and cand.N % 8 == 4 and cand.v % 2 == 0)
    tried_ells = volcano_first
    first_ells = sorted(set(need) | {2}) if volcano_first else None
    j = (get_j(p, cand.d, ells=first_ells, jobs=jobs) if volcano_first
         else get_j(p, cand.d, jobs=jobs))
    if volcano_first and j is None:
        tried_ells = False              # volcano failed: fall back to direct
        j = get_j(p, cand.d, jobs=jobs)
    while True:
        if j is not None:
            As = mont_As(j, p)
            Np = 2 * p + 2 - cand.N
            for A in As:
                side = curve_side(p, A, cand.N, Np, rng)
                if side == 0:
                    continue
                chi_target = side  # +1: E_A itself has order N; -1: its twist does
                for o, m in hits[:12]:
                    mdiv = [pr for pr in sorted(cand.sfac) if m % pr == 0]
                    odiv = mdiv + [q]
                    x = find_point(p, A, cand.N, o, odiv, chi_target, rng)
                    if x is not None:
                        return (A, x, o)
        # No usable point.  Cofactor sampling Q = [N/o]P yields exact order o
        # only when the ell-Sylow is CYCLIC for every ell | m: a split Sylow
        # (possible exactly when ell | v, i.e. the H_D curve sits above the
        # ell-volcano floor) loses one ell-power to the cofactor no matter the
        # exponent.  So retry once from the floor of every volcano that any
        # hit's filler shares with v -- the same walk the one-shot engine does
        # for each ell | m.
        if not tried_ells:
            tried_ells = True
            if need:
                log("    D=-%d: retrying via volcano floors ells=%s"
                    % (cand.d, ",".join(map(str, need))))
                j = get_j(p, cand.d, ells=need, jobs=jobs)
                continue
        return None


# ------------------------------------------------------------ the level prover
def ecm_timeout(b1, curves, nbits):
    return 60 + int(1.2e-10 * curves * b1 * nbits * nbits)


def ladder_rounds(nbits):
    """(method, B1, sweeps, per-sweep cap): each ECM tier's curve budget is
    split into small sweeps so the pool is round-robined -- every candidate
    gets its first visit at a tier before any candidate gets a second
    (admission sorts by effort first), which beats batching because success
    per curve has diminishing returns on a fixed candidate.  P-1 tiers run
    once per candidate (repeating P-1 at the same B1 finds nothing new) and
    catch q with q-1 smooth, independently of the ECM group draws.  Measured
    per-B1 cost vs one ECM curve (bench_pm1_ratio.py, deterministic inputs;
    gmp-ecm 7.0.6, Ryzen AI Max+ 395, 2026-08-25): 1/5.6 (900 bits, B1=5e5),
    1/6.9 (900, 1e7), 1/8.2 (1330, 5e5), 1/10.5 (1330, 1e7); the b1//8 weight
    below is the central value, and only the strict monotonicity of the
    cumulative keys matters for scheduling (visit order), not the constant.
    SHORTECPP_PM1=0 disables the P-1 tiers."""
    if nbits >= 1150:
        tiers = [("ecm", 2000, [12], 0), ("ecm", 11000, [20], 0),
                 ("pm1", 500000, [1], 20000),
                 ("ecm", 50000, [16, 16], 20000), ("ecm", 250000, [12] * 4, 6000),
                 ("pm1", 10000000, [1], 2000),
                 ("ecm", 1000000, [12] * 6, 2000), ("ecm", 3000000, [12] * 8, 700),
                 ("ecm", 10000000, [12] * 6, 200)]
        if os.environ.get("SHORTECPP_DEEP", "0") == "1":
            # opt-in deep tail: full-ladder veterans get one more escalation
            # (the 35-42 digit factor band is thinly covered by B1=1e7 x 72)
            tiers.append(("ecm", 30000000, [12] * 6, 100))
    elif nbits >= 900:
        tiers = [("ecm", 2000, [12], 0), ("ecm", 11000, [20], 0),
                 ("pm1", 500000, [1], 12000),
                 ("ecm", 50000, [16, 16], 12000), ("ecm", 250000, [12] * 4, 4000),
                 ("pm1", 10000000, [1], 1200),
                 ("ecm", 1000000, [12] * 6, 1200), ("ecm", 3000000, [12] * 8, 400)]
    elif nbits >= 700:
        tiers = [("ecm", 2000, [12], 0), ("ecm", 11000, [20], 0),
                 ("pm1", 500000, [1], 15000),
                 ("ecm", 50000, [16, 16], 15000), ("ecm", 250000, [12] * 4, 5000),
                 ("ecm", 1000000, [12] * 6, 1500), ("ecm", 3000000, [12] * 8, 400)]
    elif nbits >= 500:
        tiers = [("ecm", 2000, [12], 0), ("ecm", 11000, [20], 0),
                 ("pm1", 500000, [1], 0),
                 ("ecm", 50000, [16, 16], 0), ("ecm", 250000, [12] * 4, 12000),
                 ("ecm", 1000000, [12] * 6, 4000), ("ecm", 3000000, [12] * 8, 1000)]
    else:
        tiers = [("ecm", 2000, [10], 0), ("ecm", 11000, [16], 0),
                 ("pm1", 500000, [1], 0),
                 ("ecm", 50000, [14, 14], 0), ("ecm", 250000, [10] * 4, 0),
                 ("ecm", 1000000, [16] * 4, 0)]
    if os.environ.get("SHORTECPP_PM1", "1") == "0":
        tiers = [t for t in tiers if t[0] != "pm1"]
    out = []
    ecost = 0
    for method, b1, sweeps, cap in tiers:
        for c in sweeps:
            # effort keys are cumulative estimated cost, strictly increasing
            # through the ladder (a candidate visits each round once); the
            # P-1 weight is the measured cost ratio above (only monotonicity
            # is load-bearing)
            ecost += b1 * c if method == "ecm" else b1 // 8
            out.append((method, b1, c, cap, ecost))
    return out


def b_schedule(nbits, B0=None, Bmax=None):
    if B0 is None:
        B0 = {1: 30000000, 2: 60000000, 3: 120000000, 4: 400000000,
              5: 1000000000}[
            1 if nbits < 700 else 2 if nbits < 800 else 3 if nbits < 950
            else 4 if nbits < 1150 else 5]
    if Bmax is None:
        Bmax = 200000000000
    return B0, Bmax


def prove_level_cm(p, n2, small_primes, P2, threads, stats, seed=1, B0=None, Bmax=None, maxfb=0):
    nb = p.bit_length()
    rt = isqrt(p)
    L = rt + 1 + isqrt(4 * rt)
    rng = random.Random(seed * 1000003 + nb)
    B0, Bmax = b_schedule(nb, B0, Bmax)
    # The factor base is capped from the first scan (maxfb= fixes the cap,
    # otherwise it follows B) and its square roots are kept on disk for the
    # level, so a widening pays only for the primes the larger cap adds.
    fixed_cap = maxfb
    roots = roots_path(p)
    cap_prev = 0
    rounds = ladder_rounds(nb)
    t0 = time.time()
    pool = []
    winners = []       # (o, cand, q, hits)
    seen_winner_keys = set()
    B, Bmin = B0, 0
    fresh_med = None
    pmod4 = p % 4

    def note_prime_factor(c, f, fresh=True):
        """A newly known prime factor f (peeled or terminal tail)."""
        hits = accept_q(c, f, p, L, n2, small_primes)
        if hits:
            key = (c.d, c.N, f)
            if key not in seen_winner_keys:
                seen_winner_keys.add(key)
                winners.append((hits[0][0], c, f, hits))
                log("    WINNER D=-%d q %d bits, o=m*q with m=%d (|D| h~%d)"
                    % (c.d, f.bit_length(), hits[0][1], isqrt(c.d)))
                if fresh:
                    journal_write({"p": str(p), "d": c.d, "t": c.t, "v": c.v,
                                   "N": str(c.N), "q": str(f)})

    def process_result(c, found, cof):
        # duplicate-safe: concurrent batches on the same candidate may both
        # report factors; accept each found factor only if it still divides
        # the current tail, and recompute the tail directly.
        progressed = False
        for f in found:
            if f <= 1 or c.tail % f != 0:
                continue
            c.tail //= f
            c.splits += 1
            progressed = True
            if f < Y0:                      # tiny stragglers: fold into mid
                c.mid *= f
                continue
            if is_prp(f, 12):
                note_prime_factor(c, f)
            # composite or out-of-zone factors just join the discarded cofactor
        if progressed:
            c.prp_checked = False
        if c.tail == 1 or c.s * c.tail <= L:
            c.dead = True

    def settle_tail(c, prp):
        """Record the probable-prime verdict on c's current tail (from smoothtest
        or gmp-ecm); a prime tail is a dead end unless it is an acceptable q.
        Only a tail below sqrt p can be a q, so only those get the full
        Miller-Rabin test: 24 rounds in the interpreter on every prime tail was
        what the old classify phase spent its time on (about 1,500 prime tails
        per 27,000 orders, 60 ms each at 1030 bits).  A probable prime above
        sqrt p is dropped on the Fermat verdict alone; a base-2 pseudoprime
        there costs one dead branch."""
        c.prp_checked = True
        if c.dead or c.tail <= 1 or not prp:
            return
        if c.tail * c.tail >= p:
            c.dead = True
        elif is_prp(c.tail, 24):
            c.dead = True
            note_prime_factor(c, c.tail)

    def try_build():
        winners.sort(key=lambda w: w[1].d)
        while winners:
            o, c, q, hits = winners.pop(0)
            t1 = time.time()
            lev = build_level(p, c, q, hits, rng, jobs=threads)
            if lev:
                A, x, o = lev
                log("  level %d bits: built D=-%d h~%d in %.1fs; o has %d bits, q %d bits (total %.1fs, %d cands, %.0f ecm-s)"
                    % (nb, c.d, isqrt(c.d), time.time() - t1, o.bit_length(), q.bit_length(),
                       time.time() - t0, stats.get("cands", 0), stats.get("ecm_s", 0.0)))
                try:
                    os.remove(roots)          # the level is done with its scan
                except OSError:
                    pass
                return (A, x, o, q, {"D": -c.d, "t": c.t, "v": c.v, "N": str(c.N),
                                     "m": str(o // q), "bits": nb,
                                     "seconds": round(time.time() - t0, 1)})
            log("    build failed for D=-%d, trying next winner" % c.d)
        return None

    # replay journaled winners for this modulus (crash / preemption recovery)
    for rec in journal_read(JOURNAL):
        if rec.get("p") == str(p):
            c = Cand(rec["d"], rec["t"], rec["v"], int(rec["N"]))
            S0, rest = 1, c.N
            g = gcd(P2 % rest, rest)
            while g > 1:
                S0 *= g
                rest //= g
                g = gcd(g, rest)
            c.s, c.tail = S0, rest
            q = int(rec["q"])
            if c.N % q == 0 and is_prp(q, 24):
                note_prime_factor(c, q, fresh=False)
    if winners:
        log("  level %d bits: %d journaled winner(s) replayed" % (nb, len(winners)))
        lev = try_build()
        if lev:
            return lev

    while True:
        # ---- widen the candidate pool ----
        if Bmin < B:
            tds = time.time()
            fbc = fixed_cap if fixed_cap else fb_cap(B)
            recs = run_dscan(p, B, Bmin, threads, fbc, cap_prev if fbc > cap_prev else 0, roots)
            cap_prev = max(fbc, cap_prev)
            fresh = []
            for (d, t, v) in recs:
                for N in (p + 1 - t, p + 1 + t):
                    if N % 4 == 0:          # Montgomery model needs 4 | #E
                        fresh.append(Cand(d, t, v, N))
            tstrip = time.time()
            strip = smooth_strip(sorted(set(c.N for c in fresh)), threads, n2)
            for c in fresh:
                S, c.tail, c.s, prp = strip[c.N]
                c.mid = S // c.s
                if c.tail == 1 or c.s * c.tail <= L:
                    c.dead = True
                settle_tail(c, prp)
            live_new = [c for c in fresh if not c.dead]
            pool += live_new
            if live_new:
                rtb0 = rt.bit_length()
                idxs = sorted(math.log(max(math.log(max(c.s, 3)), 1.1))
                              - 0.03 * max(0, c.tail.bit_length() - rtb0)
                              for c in live_new)
                fresh_med = idxs[len(idxs) // 2]
            stats["cands"] = stats.get("cands", 0) + len(fresh)
            log("  level %d bits: dscan[%g,%g) maxfb=%g +%d D -> +%d orders (%d live) in %.1fs"
                " (strip %.1fs), %d winners"
                % (nb, Bmin, B, fbc, len(recs), len(fresh), len(live_new),
                   time.time() - tds, time.time() - tstrip, len(winners)))
            Bmin = B
            if winners:
                lev = try_build()
                if lev:
                    return lev
        # ---- ECM ladder over the pool ----
        ri = 0
        while ri < len(rounds):
            (method, b1, curves, cap, ekey) = rounds[ri]
            if winners:
                break
            # breadth-first: before the deepest tiers, prefer widening the
            # pool until the discriminant scan is well-explored
            if b1 >= 10 ** 6 and B < min(Bmax, max(8 * B0, 4 * 10 ** 9)):
                break
            live = [c for c in pool if not c.dead and c.effort < ekey]
            # Bayesian index policy: rank by log posterior win-rate per curve.
            # Coverage enters the win probability linearly (the window admits q
            # in (~sqrt(p)/s, sqrt(p))), hence ln(ln s); the tail's excess mass
            # costs ~0.03/bit (Dickman slope); each failed curve at this tier
            # multiplies the "findable factor remains here" posterior by
            # (1-p) with p ~ 1/100 at a tier's marginal digit class.  Equal
            # indices rotate (round-robin); a strong prior holds its slots
            # until ~0.01*curves of failures equalize it with the field.
            rtb = rt.bit_length()

            def index(c, extra=0):
                if c.last_b1 != b1:
                    c.tcurves, c.last_b1 = 0, b1
                return (math.log(max(math.log(max(c.s, 3)), 1.1))
                        - 0.03 * max(0, c.tail.bit_length() - rtb)
                        - 0.010 * (c.tcurves + extra))
            live.sort(key=lambda c: -index(c))
            if cap:
                live = live[:cap]
                # widen-vs-deepen economics: every failure discounts the
                # veterans; once the marginal admitted candidate is worth less
                # than a median fresh draw, the scan is the better buy (its
                # per-candidate cost is ~1 core-s vs 30-300 for deep tiers)
                if (fresh_med is not None and live and B < Bmax
                        and index(live[-1]) < fresh_med - 0.15):
                    log("    tier B1=%d: admission boundary %.2f < fresh median %.2f; widening"
                        % (b1, index(live[-1]), fresh_med))
                    break
            if method == "ecm" and len(live) < threads and live:
                # conditional duplicates: pad idle slots with extra concurrent
                # curve batches for the best candidates, each repeat valued as
                # if its earlier in-flight batches fail (index - 0.01*curves).
                # ECM only: a duplicate P-1 run on the same tail and B1 is the
                # same computation -- success depends on q-1's smoothness, not
                # on fresh randomness -- so P-1 tiers run once per candidate.
                ranked = sorted(live, key=lambda c: -index(c, extra=curves))
                pads = []
                r = 1
                while len(live) + len(pads) < threads and r <= 3:
                    for c in ranked:
                        if len(live) + len(pads) >= threads:
                            break
                        if index(c, extra=r * curves) > index(ranked[-1]) - 0.5:
                            pads.append(c)
                    r += 1
                live = live + pads
            if not live:
                ri += 1
                continue
            # Results are handled as each gmp-ecm process finishes, and the
            # round stops at the first winner: the tails still running are
            # killed and keep their old effort key, so if the build fails the
            # same tier is simply resumed for them.
            tr = time.time()
            tmo = ecm_timeout(b1, curves, nb)
            k = ecm_chunk(b1, curves, nb, method, len(live), threads)
            units = [[(c, c.tail) for c in live[i:i + k]] for i in range(0, len(live), k)]
            runner = EcmRunner()
            nsplit = ndone = 0
            with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as ex:
                futs = {ex.submit(runner.run, [t for (_, t) in u], b1, curves, tmo * len(u), method): u
                        for u in units}
                for fut in concurrent.futures.as_completed(futs):
                    res = None if fut.cancelled() else fut.result()
                    if res is None:
                        continue                    # stopped before it ran or finished
                    for (c, t) in futs[fut]:
                        ndone += 1
                        c.effort = ekey
                        c.tcurves += curves
                        if t not in res:
                            continue
                        found, cof, cprp = res[t]
                        nsplit += 1
                        process_result(c, found, cof)
                        if c.tail == cof:           # gmp-ecm already tested this cofactor
                            settle_tail(c, cprp)
                        elif not c.dead and c.tail > 1:   # a concurrent duplicate got there first
                            settle_tail(c, pow(2, c.tail - 1, c.tail) == 1)
                    if winners and not runner.stopped:
                        runner.abort()
                        for g in futs:
                            g.cancel()
            dt = time.time() - tr
            stats["ecm_s"] = stats.get("ecm_s", 0.0) + dt * threads
            log("    %s B1=%d c=%d on %d tails: %d split, %d winners, %.1fs%s"
                % (method, b1, curves, len(live), nsplit, len(winners), dt,
                   " (stopped at the first winner, %d tails done)" % ndone if runner.stopped else ""))
            if winners:
                lev = try_build()
                if lev:
                    return lev
                if runner.stopped:
                    continue                        # build failed: finish this tier
            ri += 1
        if winners:
            lev = try_build()
            if lev:
                return lev
        pool = [c for c in pool if not c.dead]
        if B >= Bmax:
            raise RuntimeError("level %d bits: exhausted dscan budget B=%d" % (nb, B))
        Bmin, B = B, min(B * 6, Bmax)
        log("  level %d bits: widening to B=%g (%d live)" % (nb, B, len(pool)))


# ------------------------------------------------------------ chain assembly
def gp_tail_chain2(q, ntop):
    """Finish the chain below CM_MIN_BITS with short2.gp (v2 caps)."""
    sgp = os.path.join(ECPP_DIR, "short2.gp")
    script = open(sgp).read() + (
        "\nSC_tlim=60;\nSC_branchcurves=200;\n"
        "lg=log(%d)/log(2); Bv=ceil(%d^2/lg);\n"
        "tv=scchain2(%d,Bv,Bv^2,ceil(%d/lg),0);\n"
        "if(type(tv)!=\"t_VEC\",print(\"FAIL\"),print(tv));\n"
        % (ntop, ntop, q, ntop))
    out = gp(script, timeout=3600)
    line = [l for l in out.splitlines() if l.strip()][-1].strip()
    if line == "FAIL" or not line.startswith("["):
        raise RuntimeError("short2.gp tail chain failed for q=%d" % q)
    return [int(t) for t in line.strip("[]").split(",")]


def gp_tail_chain(q, n2):
    """Finish the chain below CM_MIN_BITS with short.gp's SEA search."""
    sgp = os.path.join(SPP_DIR, "short.gp") if SPP_DIR else ""
    if not sgp or not os.path.exists(sgp) or "scchain(" not in open(sgp).read():
        sgp = os.path.join(ECPP_DIR, "short.gp")   # bundled copy defines scchain
    script = open(sgp).read() + (
        "\nSC_tlim=60;\nSC_branchcurves=200;\n"
        "tv=scchain(%d,%d,0);\n"
        "if(type(tv)!=\"t_VEC\",print(\"FAIL\"),print(tv));\n" % (q, n2))
    out = gp(script, timeout=3600)
    line = [l for l in out.splitlines() if l.strip()][-1].strip()
    if line == "FAIL" or not line.startswith("["):
        raise RuntimeError("short.gp tail chain failed for q=%d" % q)
    return [int(t) for t in line.strip("[]").split(",")]


def pari_check(seq, n2):
    """Independent PARI verification of every level (point on curve, exact order
    via ellmul with the known factorization, window, descent, smoothness)."""
    p = seq[0]
    small_primes = sieve_primes(n2)
    scripts = []
    for i in range(1, len(seq), 3):
        A, x, o = seq[i], seq[i + 1], seq[i + 2]
        rest = o
        fac = []
        for q in small_primes:
            if q * q > rest:
                break
            if rest % q == 0:
                fac.append(q)
                while rest % q == 0:
                    rest //= q
        if 1 < rest <= n2:
            fac.append(rest)
            rest = 1
        p_next = rest
        divs = fac + ([p_next] if p_next > 1 else [])
        rt = isqrt(p)
        L = rt + 1 + isqrt(4 * rt)
        r = min(divs)
        ok_window = (L < o < r * L) and (p_next == 1 or p_next * p_next < p)
        if not ok_window:
            return False
        scripts.append(
            "p=%d;A=%d;x=%d;o=%d;" % (p, A, x, o) +
            "B=Mod(x^3+A*x^2+x,p);if(B==0,print(\"0 0 0\"),"
            "E2=ellinit([0,lift(A/B),0,lift(1/B^2),0],p);"
            "P=[lift(Mod(x,p)/B),lift(Mod(1,p)/B)];"
            "onc=ellisoncurve(E2,P);"
            "z1=ellmul(E2,P,o)==[0];"
            "z2=1;foreach(%s,q,if(ellmul(E2,P,o/q)==[0],z2=0));"
            "print(onc,\" \",z1,\" \",z2));" % str(divs))
        p = p_next
    out = gp("\n".join(scripts) + "\n", timeout=3600)
    lines = [l.strip() for l in out.splitlines() if l.strip() and "Warning" not in l]
    return len(lines) == (len(seq) - 1) // 3 and all(l == "1 1 1" for l in lines)


def prove_chain(p0, threads, seed=1, B0=None, Bmax=None, tag=None, out=None, resume=False, maxfb=0,
                v2=False):
    global JOURNAL, V2_BCAP, V2_RADLIM, V2_QFLOOR
    if not is_prp(p0, 32):
        sys.exit("p0 is not prime")
    n = p0.bit_length()
    n2 = n * n
    if v2:
        from math import log2 as _log2, ceil as _ceil
        V2_BCAP = _ceil(n * n / _log2(n))
        V2_RADLIM = _ceil(n / _log2(n))
        V2_QFLOOR = V2_BCAP * V2_BCAP
        log("revised format: B = %d, log2 rad <= %d, recursion floor B^2 = %d"
            % (V2_BCAP, V2_RADLIM, V2_QFLOOR))
    log("proving p0 (%d bits), n2 = %d" % (n, n2))
    t0 = time.time()
    small_primes = sieve_primes(n2)
    P2 = balanced_product(small_primes)
    name = tag or ("p%dbits" % n)
    os.makedirs(CAMP_DIR, exist_ok=True)
    JOURNAL = os.path.join(CAMP_DIR, name + ".winners.jsonl")
    partial = os.path.join(CAMP_DIR, name + ".partial.json")
    seq = [p0]
    levels_meta = []
    p = p0
    stats = {}
    if resume and os.path.exists(partial):
        st = json.load(open(partial))
        if st.get("fmt", 1) != (2 if v2 else 1):
            log("partial state has fmt=%s, need %d -- ignoring it"
                % (st.get("fmt", 1), 2 if v2 else 1))
        elif st.get("p0") == str(p0):
            seq = [int(x) for x in st["seq"]]
            levels_meta = st.get("meta", [])
            p = int(st["p"])
            log("resumed: %d level(s) done, current p has %d bits"
                % ((len(seq) - 1) // 3, p.bit_length()))
    while p.bit_length() >= CM_MIN_BITS:
        A, x, o, q, meta = prove_level_cm(p, n2, small_primes, P2, threads, stats,
                                          seed=seed, B0=B0 if p == p0 else None,
                                          Bmax=Bmax if p == p0 else None,
                                          maxfb=maxfb if p == p0 else 0)
        seq += [A, x, o]
        levels_meta.append(meta)
        p = q
        with open(partial, "w") as f:
            json.dump({"p0": str(p0), "fmt": 2 if v2 else 1,
                       "seq": [str(v) for v in seq], "p": str(p),
                       "meta": levels_meta}, f)
    if p > 1:
        log("finishing below %d bits with %s (p has %d bits)"
            % (CM_MIN_BITS, "short2.gp" if v2 else "short.gp", p.bit_length()))
        t1 = time.time()
        tail = gp_tail_chain2(p, n) if v2 else gp_tail_chain(p, n2)
        log("  gp tail: %d levels in %.1fs" % (len(tail) // 3, time.time() - t1))
        seq += tail
    dt = time.time() - t0
    log("chain complete: %d levels, %.1fs; verifying..." % ((len(seq) - 1) // 3, dt))
    if v2:
        import importlib.util
        _s2 = importlib.util.spec_from_file_location(
            "vshort2", os.path.join(ECPP_DIR, "vshort2.py"))
        _v2mod = importlib.util.module_from_spec(_s2); _s2.loader.exec_module(_v2mod)
        ok = _v2mod.verify(seq)               # the revised format's verifier is the gate
        okp = pari_check(seq, V2_BCAP)        # parse at B so the split matches the format
        log("vshort2.verify = %s   PARI cross-check = %s" % (ok, okp))
    else:
        ok = verify(seq)
        okp = pari_check(seq, n2)
        log("vsmallECPP.verify = %s   PARI cross-check = %s" % (ok, okp))
    if not (ok and okp):
        sys.exit("VERIFICATION FAILED -- certificate not written")
    if out is None:
        out = os.path.join(ECPP_DIR, "..", "certs", "short",
                           (tag or ("p%dbits" % n)) + ".txt")
    outdir = os.path.dirname(out)
    if outdir:
        os.makedirs(outdir, exist_ok=True)
    with open(out, "w") as f:
        f.write(" ".join(map(str, seq)) + "\n")
    with open(out + ".json", "w") as f:
        json.dump({"p0": str(p0), "bits": n, "levels": (len(seq) - 1) // 3,
                   "seconds": round(dt, 1), "threads": threads,
                   "stats": {k: (round(v, 1) if isinstance(v, float) else v)
                             for k, v in stats.items()},
                   "cm_levels": levels_meta,
                   "seq_bits": sum(v.bit_length() for v in seq)}, f, indent=1)
    log("certificate written to %s (verified twice)" % out)
    return seq


def main():
    global VERBOSE
    try:                                  # one gmp-ecm process and two pipes per worker thread
        import resource
        soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
        resource.setrlimit(resource.RLIMIT_NOFILE, (min(hard, 65536) if hard > 0 else 65536, hard))
    except (ImportError, ValueError, OSError):
        pass
    args = dict(kv.split("=", 1) for kv in sys.argv[1:] if "=" in kv)
    if "p" in args:
        p0 = int(args["p"])
        tag = args.get("tag")
    elif "c" in args:
        c = int(args["c"])
        p0 = int(gp("print(nextprime(10^%d))" % c).split()[-1])
        tag = args.get("tag", "nextprime1e%d" % c)
    else:
        print(__doc__)
        sys.exit(2)
    VERBOSE = int(args.get("v", "1"))
    threads = int(args.get("threads", str(os.cpu_count())))
    seed = int(args.get("seed", "1"))
    B0 = int(float(args["B0"])) if "B0" in args else None
    Bmax = int(float(args["Bmax"])) if "Bmax" in args else None
    prove_chain(p0, threads, seed=seed, B0=B0, Bmax=Bmax, tag=tag, out=args.get("out"),
                resume=bool(int(args.get("resume", "0"))),
                maxfb=int(float(args.get("maxfb", "0"))),
                v2=bool(int(args.get("v2", "0"))))


if __name__ == "__main__":
    main()
