#!/usr/bin/env python3
"""shortECPP.py -- CM-method prover for short ECPP certificates
(github.com/AndrewVSutherland/ShortPrimalityProofs format).

At each level (modulus p, top-level n fixed for the whole chain, n2 = n^2):
  1. dscan enumerates (D,t,v) with 4p = t^2 + |D| v^2; candidate curve orders
     are N = p+1-t and p+1+t with 4 | N (Montgomery-representable side pairs).
  2. smoothtest strips every candidate's y0-smooth part S (y0 = 2^26, less at
     small levels, see level_params; batched Bernstein remainder trees against
     a cached prime product); the n2-smooth
     part s | S supplies the certificate cofactor m, the rest of S is discarded
     into the untracked cofactor c (which never needs factoring), and the rough
     tail T = N/S is peeled by gmp-ecm.  Which tail gets which curves is decided
     by the stopping rule of shortsched.py (levels of POLICY_MIN_BITS bits and
     more) or by tiered rounds with admission caps (smaller levels, sched=old).
  3. A candidate wins when a known prime q (a peeled factor, or a tail that
     became prime) satisfies n2 < q, q^2 < p, and some divisor m | s lands
     o = m*q in the window L < o < r(m)*L, L = (p^(1/4)+1)^2.  The winner's
     curve comes from H_D (cm_method / classpoly, PARI polclass fallback), the
     Montgomery A from the u-cubic, and the point of exact order o from the
     x-only ladder; the certificate stores (A, x, o) and descends to q.
Levels too small for a q above the recursion floor (cm_min_bits: 64 to 94 bits
in the revised format) are finished by short.gp, called with the top-level n2.
The first scan, the strip bound and the number of threads follow the size of
the level (level_params), so that a 170-bit prime costs half a second.  The chain is verified by vshort2.verify (v2=1, the
revised format; vsmallECPP.verify for the original) and independently
cross-checked level by level in PARI before it is written.

usage: shortECPP.py (p=<decimal> | c=<digits: p=nextprime(10^c)>)
                      [threads=N] [out=<file>] [tag=<name>] [B0=<dscan start>]
                      [Bmax=<dscan cap>] [maxfb=<factor-base cap>] [seed=1]
                      [v=1] [v2=1] [resume=1] [sched=auto|rule|old]
                      [engine=auto|gmp|yafu] [yafu=<binary>]

The scan's factor base is always capped: at maxfb if given (top level only),
otherwise at fb_cap(B) for the first scan.  From there the stopping-rule levels
get their orders the cheaper of two ways, measured as they go: a wider
discriminant bound, or a higher cap over the range already scanned (the
discriminants passed over for a prime above the cap).  The square roots are
cached per modulus under work/rootcache/ for the life of the level.

sched: levels of POLICY_MIN_BITS (280) bits and more are run by the stopping-rule
scheduler (shortsched.py: candidates are worked in the order of their
break-even reward and only while that is below what a fresh order from the scan
is worth); smaller levels, and everything with sched=old, by the tiered rounds
of ladder_rounds.  sched=rule fails if numpy is missing instead of falling back.

engine: the ECM program of the stopping-rule levels.  yafu (its AVX-512 ECM,
8 curves per call; binary from yafu= or SHORTECPP_YAFU, built with
USE_AVX512IFMA=1) costs several times less per curve than gmp-ecm and has its
own success table and cost table (shortsim.ENGINES); auto uses it when a binary
is given.  The tiered rounds and P-1 always use gmp-ecm.
"""

import concurrent.futures
import json
import math
import os
import queue
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
CM_MIN_BITS = 0                   # below this, short.gp (SEA) finishes the chain; 0: cm_min_bits()


def cm_min_bits(floor):
    """Smallest level given to the CM search.  Its q must exceed the recursion
    floor (n^2, or B^2 in the revised format: a terminal prime below it comes
    from the gp tail only) and is about half the level less the cofactor m, so
    the level needs twice the floor's bits and room for m."""
    return CM_MIN_BITS if CM_MIN_BITS else max(48, 2 * floor.bit_length() + 24)


def level_params(nb, threads, bound):
    """(strip bound in bits, first scan bound, threads) for a level of nb bits.
    The defaults of the large levels are fixed costs that swamp a small one:
    the 2^26 strip is a second of work whatever the number of orders, a scan to
    3e7 buys thousands of orders where a 200-bit level needs fifty, and
    processes with more threads than there is work spend their time waiting.
    bound: the smoothness bound of the certificate, which the strip must reach."""
    y0 = 26 if nb >= 320 else 20 if nb >= 200 else 16
    y0 = max(y0, bound.bit_length())
    B0 = (10 ** 5 if nb < 200 else 3 * 10 ** 5 if nb < 260 else 10 ** 6 if nb < 320 else
          3 * 10 ** 6 if nb < 400 else 10 ** 7 if nb < 500 else None)
    thr = max(1, min(threads, int(2.0 ** ((nb - 180) / 40.0))))
    return y0, B0, thr


def set_strip_bound(bits):
    global Y0, PCACHE
    Y0 = 1 << bits
    PCACHE = os.path.join(REPO_DIR, "work/pcache/oneshot_P_%d.bin" % Y0)

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

SCHED = os.environ.get("SHORTECPP_SCHED", "auto")    # auto | rule | old (sched= argument)
POLICY_MIN_BITS = 280             # below this a level costs less than the rule's set-up
ENGINE = os.environ.get("SHORTECPP_ECM", "auto")     # auto | gmp | yafu (engine= argument)
YAFU_BIN = os.environ.get("SHORTECPP_YAFU", "")      # yafu= argument
_YAFU_CWD = None                  # yafu reads yafu.ini from its working directory: give it an empty one
# The strip is mostly serial: 4,000 orders of 768 bits take 1.65 s on one thread,
# 0.41 s on 8 and on 16, 0.53 s on 32; 40,000 take 2.8 s on 8 and 2.7 s on 32.  The
# rule scheduler therefore gives it 8 threads and frees the scan's other threads
# for ECM while it runs (on 192 cores it used to hold 50 or more for a second per
# chunk).
STRIP_THREADS = 8

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
                 "splits", "effort", "prp_checked", "tcurves", "last_b1",
                 "cls", "tb", "st", "sent", "done", "infl", "ver")     # shortsched's fields

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
def fb_cap_max(nb):
    """Largest factor-base cap whose square roots mod a prime of nb bits fit in a
    quarter of the machine's memory (dscan keeps them all: nb/8 bytes each and
    about 64 of housekeeping, for half the primes below the cap)."""
    try:
        ram = os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")
    except (ValueError, OSError, AttributeError):
        ram = 16 << 30
    roots = 0.25 * ram / (nb // 8 + 64)
    y = 40.0 * roots
    for _ in range(6):                       # pi(y) / 2 = roots
        y = 2.0 * roots * (math.log(y) - 1.0)
    return int(min(y, 10 ** 11))


def fb_cap(B, nb=0):
    """Factor-base cap for a first scan to B, and for the scans of the tiered
    rounds.  A square root mod p costs about 60 Cornacchia tests and a prime q
    is in about B/(30 q) of the discriminants scanned, so roots of primes far
    above B/1000 never pay for themselves; the count of solvable D is flat in
    the cap (measured at 1030 bits, B = 2e10: a cap of B/4000 keeps 81% of what
    B/400 finds, with a ninth of the roots).  The stopping-rule levels start
    here and then raise the cap when that is the cheaper way to more orders
    (RuleDriver.next_chunk)."""
    return max(10 ** 6, min(fb_cap_max(nb) if nb else 10 ** 8, B // 1000))


_DSCAN_SUMMARY = re.compile(r"dscan B=\d+ Bmin=\d+: fb ([\d.]+)s \(sieve ([\d.]+), legendre ([\d.]+), tonelli ([\d.]+); "
                            r"(\d+) roots, (\d+) cached\), scan ([\d.]+)s, (\d+)/(\d+) solvable, (\d+) threads(?:, (\d+) tested)?")


def roots_path(p):
    """Per-modulus cache of the scan's factor-base square roots (dscan roots=)."""
    d = os.path.join(REPO_DIR, "work", "rootcache")
    os.makedirs(d, exist_ok=True)
    return os.path.join(d, "roots_%s.bin" % hashlib.sha256(str(p).encode()).hexdigest()[:20])


def run_cap(cmd, inp=None, reg=None):
    """Run cmd to completion and return (exit code, stdout, stderr).  While it
    runs the process sits in reg (a set), so another thread can kill it."""
    pr = subprocess.Popen(cmd, stdin=subprocess.PIPE if inp is not None else subprocess.DEVNULL,
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, env=ECM_ENV)
    if reg is not None:
        reg.add(pr)
    try:
        out, err = pr.communicate(inp)
    finally:
        if reg is not None:
            reg.discard(pr)
    return pr.returncode, out, err


def run_dscan(p, B, Bmin, threads, maxfb=0, fbmin=0, roots=None, reg=None, info=None):
    """Solvable (d, t, v) with Bmin <= d < B over the factor base capped at maxfb.
    With fbmin set (the cap of the scans already done below Bmin) the d < Bmin
    that need a prime in (fbmin, maxfb] are included too, so successive calls
    with a growing cap tile [0, B) at the latest cap without repeats; with
    Bmin = B the call only raises the cap.  info (a dict) receives dscan's own
    account: seconds for the square roots and for the scan, roots computed,
    discriminants tested."""
    cmd = [os.path.join(ECPP_DIR, "dscan"), "p=%d" % p, "B=%d" % B, "threads=%d" % threads, "dump"]
    if Bmin:
        cmd.append("Bmin=%d" % Bmin)
    if maxfb:
        cmd.append("maxfb=%d" % maxfb)
    if fbmin and Bmin:
        cmd.append("fbmin=%d" % fbmin)
    if roots:
        cmd.append("roots=" + roots)
    rc, sout, serr = run_cap(cmd, None, reg)
    if rc != 0:
        raise RuntimeError("dscan failed: %s" % serr[:500])
    m = _DSCAN_SUMMARY.search(serr)
    if info is not None and m:
        info.update(fb_s=float(m.group(1)), root_s=float(m.group(4)), roots_new=int(m.group(5)) - int(m.group(6)),
                    scan_s=float(m.group(7)), solvable=int(m.group(8)), threads=int(m.group(10)),
                    tested=int(m.group(11) or m.group(9)))
    out = []
    for line in sout.splitlines():
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


def smooth_strip(Ns, threads, n2, reg=None):
    """Batched strip: N -> (S, T, s, prp) with N = S*T, S the y0-smooth part
    (cached prime product, smoothtest parts), s the n2-smooth part of S with
    multiplicity, and prp the base-2 Fermat verdict on T.  The last two are
    computed inside smoothtest, in parallel, so nothing per-candidate is left
    for the interpreter."""
    if not Ns:
        return {}
    ensure_pcache(threads)
    inp = "\n".join(str(N) for N in Ns) + "\n"
    rc, sout, serr = run_cap([os.path.join(ECPP_DIR, "smoothtest"), "parts", "load=" + PCACHE,
                              "threads=%d" % threads, "n2=%d" % n2, "prp=1"], inp, reg)
    if rc != 0:
        raise RuntimeError("smoothtest failed: %s" % serr[:500])
    res = {}
    for line in sout.splitlines():
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


_YAFU_EXPR = re.compile(r"^(?:>> )?=== Starting work on batchfile expression #(\d+) ===")
_YAFU_PART = re.compile(r"^(P|PRP|C)\d+ = (\d+)\s*$")


_YAFU_PROG = re.compile(r"^ecm: (\d+)/(\d+) curves on ")


class EcmResult(dict):
    """{n: (factors, cofactor, cofactor_is_prp)} for the inputs that were split.
    ran: {n: curves run} for the inputs on which the engine stopped before the
    curves asked for were done."""

    def __init__(self):
        dict.__init__(self)
        self.ran = {}


def parse_yafu(out, nums):
    """yafu's report for the expressions ecm(n, curves), one per input, in the
    form parse_ecm returns: {n: (factors, cofactor, cofactor_is_prp)} for the
    inputs it split.  yafu lists the parts it ended with (P/PRP = probable
    prime, C = composite); the largest is taken as the cofactor, and a split is
    kept only if the parts multiply back to the input.

    yafu stops working on an input after the batch of 8 curves in which a factor
    appears, whatever is left of the curves asked for and even if the cofactor
    is composite (gmp-ecm goes on with the cofactor).  Its progress lines
    "ecm: K/TOTAL curves" give the curves it ran; they are returned in .ran for
    the inputs that stopped short, so the caller can credit the right number."""
    res = EcmResult()
    cur, parts, prog = None, [], None

    def flush():
        if cur is not None and len(parts) >= 2 and math.prod(v for (_, v) in parts) == cur:
            k = max(range(len(parts)), key=lambda i: parts[i][1])
            res[cur] = ([v for i, (_, v) in enumerate(parts) if i != k], parts[k][1], parts[k][0] != "C")
            if prog is not None and prog[0] < prog[1]:
                res.ran[cur] = prog[0]
    for line in out.replace("\r", "\n").splitlines():
        m = _YAFU_EXPR.match(line)
        if m:
            flush()
            k = int(m.group(1)) - 1
            cur, parts, prog = (nums[k] if 0 <= k < len(nums) else None), [], None
            continue
        m = _YAFU_PART.match(line)
        if m:
            parts.append((m.group(1), int(m.group(2))))
            continue
        m = _YAFU_PROG.match(line)
        if m:
            prog = (int(m.group(1)), int(m.group(2)))
    flush()
    return res


def yafu_cwd():
    global _YAFU_CWD
    if _YAFU_CWD is None:
        import atexit
        import tempfile
        _YAFU_CWD = tempfile.mkdtemp(prefix="shortecpp_yafu_")
        atexit.register(shutil.rmtree, _YAFU_CWD, True)
    return _YAFU_CWD


_ENGINE_OK = None


def ecm_engine():
    """The ECM program of the stopping-rule levels: "yafu" or "gmp".  yafu is
    used when a binary is given and passes a self-test (a planted 31-bit factor
    of a 450-bit number must come back, with the cofactor called prime); a
    binary that fails it -- wrong build for this CPU, changed output format --
    is refused with engine=yafu and replaced by gmp-ecm with engine=auto."""
    global YAFU_BIN, _ENGINE_OK
    if ENGINE == "gmp":
        return "gmp"
    if _ENGINE_OK is not None:
        return _ENGINE_OK
    if not YAFU_BIN and ENGINE == "yafu":
        YAFU_BIN = shutil.which("yafu") or ""
    why = "no yafu binary (yafu=<path> or SHORTECPP_YAFU)"
    if YAFU_BIN and os.access(YAFU_BIN, os.X_OK):
        f, g = 1191999629, (1 << 419) + 299          # both prime
        try:
            res = EcmRunner("yafu").run([f * g], 11000, 24, 120)
        except OSError as e:
            res, why = None, "yafu did not run: %r" % (e,)
        if res and res.get(f * g) == ([f], g, True):
            _ENGINE_OK = "yafu"
            return _ENGINE_OK
        if res is not None:
            why = "yafu failed its self-test (%s)" % YAFU_BIN
    if ENGINE == "yafu":
        sys.exit("engine=yafu: " + why)
    if YAFU_BIN:
        log("  %s; using gmp-ecm" % why)
    _ENGINE_OK = "gmp"
    return _ENGINE_OK


class EcmRunner:
    """Runs ECM processes for one round and can stop them all at once.

    Each call handles a list of numbers in one process (gmp-ecm takes one input
    per line, yafu one expression per line), which keeps the launch rate down at
    the shallow tiers where a tail needs a fraction of a second.  abort() kills
    every running process and makes pending calls return None, so a round ends
    as soon as a winner is in hand instead of when its slowest tail finishes.

    engine "yafu": yafu's AVX-512 ECM, 8 curves at a time (curves is rounded up
    to a multiple of 8 by yafu), FFT stage 2 with a B2 that depends on B1 only
    (-avxecm_stg2 1), APR-CL proofs of prime cofactors switched off (-aprcl_p 1:
    two seconds each at 1000 bits, and the prover tests them itself).  The curves
    are seeded from the clock, the process id and the system's random source, so
    processes running at once use different curves.  P-1 always goes to gmp-ecm."""

    def __init__(self, engine="gmp"):
        self.lock = threading.Lock()
        self.procs = set()
        self.stopped = False
        self.engine = engine

    def run(self, nums, b1, curves, timeout, method="ecm"):
        if self.stopped:
            return None
        if self.engine == "yafu" and method == "ecm":
            return self._run_yafu(nums, b1, curves, timeout)
        cmd = [ECM_BIN] + (["-pm1"] if method == "pm1" else []) + ["-c", str(curves), str(b1)]
        out = self._communicate(cmd, "\n".join(str(n) for n in nums) + "\n", timeout, None)
        return None if out is None else parse_ecm(out, nums)

    def _run_yafu(self, nums, b1, curves, timeout):
        """All `curves` curves on every number.  yafu drops an input after the
        batch in which a factor appears, so the curves it left are run on the
        composite cofactor in a further process, until none are left (gmp-ecm
        does the same within one call)."""
        cmd = [YAFU_BIN, "-B1ecm", str(b1), "-threads", "1", "-avxecm_stg2", "1", "-aprcl_p", "1",
               "-logfile", os.devnull, "-session", os.devnull]
        res = {}
        todo = [(n, n, curves) for n in nums]         # (input, what is left of it, curves left)
        t_end = time.time() + timeout
        while todo:
            out = self._communicate(cmd, "".join("ecm(%d,%d)\n" % (cur, left) for (_, cur, left) in todo),
                                    max(1.0, t_end - time.time()), yafu_cwd())
            if out is None:
                return None
            r = parse_yafu(out, [cur for (_, cur, _) in todo])
            nxt = []
            for (n, cur, left) in todo:
                if cur not in r:
                    continue                          # all its curves ran and found nothing
                found, cof, cprp = r[cur]
                if n in res:                          # the same input twice in one call
                    found = res[n][0] + [f for f in found if f not in res[n][0]]
                    if res[n][1] < cof:
                        cof, cprp = res[n][1], res[n][2]
                res[n] = (found, cof, cprp)
                left -= r.ran.get(cur, left)
                if left > 0 and cof > 1 and not cprp:
                    nxt.append((n, cof, left))
            todo = nxt
            if time.time() >= t_end:
                break
        return res

    def _communicate(self, cmd, inp, timeout, cwd):
        """Output of one process, or None if the round was stopped."""
        p = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                             stderr=subprocess.DEVNULL, text=True, env=ECM_ENV, cwd=cwd)
        with self.lock:
            if self.stopped:
                p.kill()
            else:
                self.procs.add(p)
        try:
            out, _ = p.communicate(inp, timeout=timeout)
        except subprocess.TimeoutExpired:
            p.kill()
            out, _ = p.communicate()          # keep whatever it had finished
        finally:
            with self.lock:
                self.procs.discard(p)
        return None if self.stopped else out

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
_CM_UNBUILDABLE = set()          # discriminants whose class polynomial cannot be computed on this installation

# What a build will cost, before it is run.  It is root-finding at the larger
# part x of the most balanced factorisation of the class number (cm_method's
# decomposition splits there): 0.0145 x^0.93 (bits/768)^2.74 seconds on a
# 192-core Zen 5 box, within a factor 1.25 over 31 builds of 500 to 1064 bits
# (reports/shortecpp-improve/data/build_fit.py).  _BUILD_SCALE carries the
# ratio seen in this run's own builds to the next estimates.
BUILD_EST = {}                   # |D| -> estimated seconds (try_build takes the cheapest winner first)
_BUILD_SCALE = [1.0]
_LAST_ROOT = {}                  # |D| -> root-finding seconds cm_method reported


def build_estimate(ds, nb):
    """Estimated build seconds for the discriminants -d, d in ds, at a level of
    nb bits ({d: seconds}; class numbers from PARI, a fraction of a millisecond
    each).  A discriminant PARI did not answer for is left out."""
    ds = [d for d in dict.fromkeys(ds)]
    if not ds:
        return {}
    script = "bal(h)=my(v=divisors(h));h/v[(#v+1)\\2];\n" + "".join("print(%d,\" \",bal(qfbclassno(-%d)));\n" % (d, d) for d in ds)
    out = {}
    try:
        for line in gp(script, timeout=300).splitlines():
            w = line.split()
            if len(w) == 2:
                out[int(w[0])] = _BUILD_SCALE[0] * (0.0145 * int(w[1]) ** 0.93 + 0.5) * (nb / 768.0) ** 2.74
    except (RuntimeError, ValueError, subprocess.TimeoutExpired):
        pass
    return out


def get_j(p, d, ells=None, jobs=0, attempts=3, plain=False):
    """j-invariant of a curve with CM by -d over F_p (cm_method, polclass
    fallback).  The root-find's equal-degree splitting is randomized and an
    attempt can come back rootless (exit 0, no j line), so retry a few times
    before giving up on the discriminant.  plain=True makes cm_method use the
    class polynomial itself instead of its class-group decomposition."""
    if d in _CM_UNBUILDABLE:
        return None
    cmd = [os.path.join(ECPP_DIR, "cm_method"), "D=-%d" % d, "p=%d" % p]
    if plain:
        cmd.append("tower=0")
    if ells:
        cmd.append("ells=" + ",".join(str(e) for e in ells))
    if jobs > 1 and d >= 20000000:
        cmd.append("jobs=%d" % min(jobs, 64))
    for att in range(attempts):
        # distinct seed per attempt: cm_method's root-find is deterministically
        # seeded, so an identical rerun would repeat a rootless result exactly
        # an hour, or several times the estimate where there is one (root-finding
        # at a prime degree of 10^5 takes an hour at 1300 bits)
        tmo = max(3600, int(6 * BUILD_EST.get(d, 0.0)) + 1800)
        try:
            r = subprocess.run(cmd + ["seed=%d" % (att + 1)],
                               capture_output=True, text=True, env=ECM_ENV, timeout=tmo)
        except subprocess.TimeoutExpired:
            log("    cm_method D=-%d: no answer after %d s, giving this discriminant up" % (d, tmo))
            _CM_UNBUILDABLE.add(d)
            return None
        if r.returncode == 0:
            for line in r.stdout.splitlines():
                w = line.split()
                if len(w) == 2 and w[0] == "j":
                    for ln in r.stderr.splitlines():      # degrees, timings, and whether a retry was needed
                        if ln.startswith("[H_D"):
                            log("    cm_method D=-%d: %s" % (d, ln.strip("[] ")))
                            m = re.search(r"roots?: ([\d.]+)s(?: \+ ([\d.]+)s)?", ln)
                            if m:
                                _LAST_ROOT[d] = float(m.group(1)) + float(m.group(2) or 0)
                    return int(w[1])
        log("    cm_method D=-%d attempt %d/%d: exit %d, no j; tail: %s" % (
            d, att + 1, attempts, r.returncode,
            (r.stderr or r.stdout)[-260:].replace("\n", " | ")))
        m = re.search(r"Unable to locate file (\S+)", r.stderr or "")
        if m:
            # nothing random about this one: the class group of D needs an isogeny
            # degree whose modular polynomial is not installed
            _CM_UNBUILDABLE.add(d)
            log("    D=-%d cannot be built here: %s is missing (the repository carries a subset of the "
                "modular polynomials, see INSTALL; CLASSPOLY_PHI_DIR should hold the full set)"
                % (d, os.path.basename(m.group(1))))
            return None
    if ells is None and d <= 3 * 10 ** 6:
        out = gp("r=polrootsmod(polclass(-%d),%d);\n"
                 "if(#r==0,print(\"NOROOT\"),print(lift(r[1])));\n" % (d, p), timeout=1800)
        w = out.split()
        if w and w[-1] != "NOROOT":
            return int(w[-1])
    return None


def j_order_ok(p, j, N):
    """Whether the curves with j-invariant j over F_p have N or 2p+2-N points:
    a random point of y^2 = x^3 + 3kx + 2k, k = j/(1728-j), is killed by one of
    the two.  (False for a wrong j, up to a negligible chance.)"""
    if j % p in (0, 1728 % p):
        return True
    out = gp("p=%d;\nk=Mod(%d,p)/(1728-Mod(%d,p));\nE=ellinit([3*k,2*k]);\nP=random(E);\n"
             "print(ellmul(E,P,%d)==[0]||ellmul(E,P,%d)==[0]);\n" % (p, j, j, N, 2 * p + 2 - N))
    w = out.split()
    return bool(w) and w[-1] == "1"


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
    ells = first_ells                   # what the current j was computed with
    j = (get_j(p, cand.d, ells=first_ells, jobs=jobs) if volcano_first
         else get_j(p, cand.d, jobs=jobs))
    if volcano_first and j is None:
        tried_ells = False              # volcano failed: fall back to direct
        ells = None
        j = get_j(p, cand.d, jobs=jobs)
    tried_plain = False
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
                    # An odd prime ell | gcd(m, v) divides pi - 1 in the maximal
                    # order (ell | N and ell | v force ell | t - 2), so the curve
                    # from H_D has full ell-torsion and no cofactor multiple has
                    # exact order o: such a hit can only work from the floor.
                    if not tried_ells and any(pr > 2 and cand.v % pr == 0 for pr in mdiv):
                        continue
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
                ells = need
                j = get_j(p, cand.d, ells=need, jobs=jobs)
                continue
        # Still nothing.  cm_method's default route is the class-group
        # decomposition, whose height bound is heuristic: a bound that was too
        # small and slipped past cm_method's own checks gives a j whose curves
        # have the wrong order.  That is cheap to see (a random point is killed
        # by neither N nor the twist's order), and a winner is too dear to drop
        # on it, so in that case compute j once more from the plain polynomial.
        if j is not None and not tried_plain and not j_order_ok(p, j, cand.N):
            tried_plain = True
            log("    D=-%d: the curves with the computed j do not have order N; "
                "recomputing j from the plain class polynomial" % cand.d)
            j = get_j(p, cand.d, ells=ells, jobs=jobs, plain=True)
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
        # a small level is cheaper to widen than to work deep: the cofactor to
        # clear has about nbits/2 bits
        top = 11000 if nbits < 230 else 50000 if nbits < 320 else 250000 if nbits < 410 else 1000000
        tiers = [t for t in tiers if t[1] <= top or (t[0] == "pm1" and nbits >= 320)]
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
    Bmax_given = Bmax
    y0b, B0s, threads = level_params(nb, threads, V2_BCAP if V2_BCAP else n2)
    set_strip_bound(y0b)
    B0, Bmax = b_schedule(nb, B0 if B0 else B0s, Bmax)
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
        if f * c.s <= L:
            return                        # no m | s lifts o = m*f above L
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
        winners.sort(key=lambda w: (BUILD_EST.get(w[1].d, 0.0), w[1].d))
        while winners:
            o, c, q, hits = winners.pop(0)
            t1 = time.time()
            lev = build_level(p, c, q, hits, rng, jobs=threads)
            est, act = BUILD_EST.get(c.d), _LAST_ROOT.get(c.d)
            if est and act and est > 2.0 and act > 2.0:
                # this machine against the fit, carried at half weight
                _BUILD_SCALE[0] = min(8.0, max(0.25, _BUILD_SCALE[0] * (act / est) ** 0.5))
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

    if use_rule(nb):
        return run_rule_level(p, nb, L, n2, small_primes, threads, stats, B0, Bmax_given, fixed_cap, roots,
                              winners, process_result, settle_tail, try_build)

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


# ------------------------------------------------------------ the stopping-rule scheduler
_HAVE_RULE = None
_MACHINE = {}


def use_rule(nb):
    """Whether this level is run by shortsched (see the module docstring)."""
    global _HAVE_RULE
    if SCHED == "old" or nb < POLICY_MIN_BITS:
        return False
    if _HAVE_RULE is None:
        try:
            import numpy                      # noqa: F401
            import shortsim                   # noqa: F401
            import shortsched                 # noqa: F401
            _HAVE_RULE = True
        except ImportError as e:
            if SCHED == "rule":
                raise
            log("stopping-rule scheduler unavailable (%s); using the tiered rounds" % e)
            _HAVE_RULE = False
    return _HAVE_RULE


def ecm_machine_factor(nb, threads, engine="gmp"):
    """Seconds per ECM curve here, with `threads` processes at once, relative to
    the engine's table in shortsim.ENGINES (a c8a.48xlarge with every core
    busy).  Measured once per run: two curves (yafu: one batch of 8) at
    B1 = 50000 on a number of nb bits with no factor in reach, one process per
    thread."""
    if engine in _MACHINE:
        return _MACHINE[engine]
    if nb < 400:
        # the benchmark would cost more than a small level; the table's costs
        # are right within a factor 1.4 on the machines seen so far
        return 1.0
    import shortsim as ss
    prof = ss.ENGINES[engine]
    rng = random.Random(12345)

    def prime(bits):
        while True:
            n = rng.getrandbits(bits) | (1 << (bits - 1)) | 1
            if is_prp(n, 8):
                return n
    N = prime(nb // 2) * prime(nb - nb // 2)
    b1, curves = 50000, 2 * prof["quantum"]

    def one(_):
        t = time.time()
        EcmRunner(engine).run([N], b1, curves, 600)
        return time.time() - t
    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as ex:
        ts = sorted(ex.map(one, range(threads)))
    per = ts[len(ts) // 2] / curves
    f = per / ss.make_cost(prof["sec"], prof["exponent"], prof["step"], scalar=True)(b1, nb)
    _MACHINE[engine] = f = min(20.0, max(0.2, f))
    log("  %s here: %.4f s per curve at B1=%d on %d bits with %d processes -> cost factor %.2f"
        % ("yafu" if engine == "yafu" else "gmp-ecm", per, b1, nb, threads, f))
    return f


class RuleDriver:
    """shortsched's driver for the prover: dscan and smoothtest supply the orders,
    gmp-ecm runs the units.  Everything that touches candidate state happens in
    the thread that runs level_loop; worker threads only run processes and post
    what came back."""

    def __init__(self, p, nb, L, threads, stats, B0, Bmax, fixed_cap, roots, nsm, small_primes,
                 winners, process_result, settle_tail, try_build, grow=2.0, engine="gmp"):
        self.p, self.nb, self.L, self.threads, self.stats = p, nb, L, threads, stats
        self.engine = engine
        self.small_primes = small_primes
        self.ecm_seen = 0.0
        self.Bmax, self.fixed_cap, self.roots, self.nsm, self.grow = Bmax, fixed_cap, roots, nsm, grow
        self.winners, self.process_result, self.settle_tail, self.try_build = \
            winners, process_result, settle_tail, try_build
        self.t0 = time.time()
        self.q = queue.Queue()
        self.ex = concurrent.futures.ThreadPoolExecutor(max_workers=threads + 4)
        self.runner = EcmRunner(engine)
        self.procs = set()
        self.epoch = 0
        # supply: the scan has two bounds, B on the discriminant and the cap on its
        # prime factors, and more orders come from moving either (next_chunk)
        self.B, self.cap = 0, 0
        self.B0 = min(B0, Bmax)
        self.capmax = fixed_cap if fixed_cap else fb_cap_max(nb)
        self.arm = ("widen", self.B0, fixed_cap if fixed_cap else fb_cap(self.B0, nb))
        self.last = None                      # the last widening: (orders, scan core-s, strip core-s, lo, hi, density then)
        self.raise_fit = (0.7, 1.0)           # found/predicted of the last raise: orders, core-seconds
        self.tot_orders, self.tot_test_s = 0, 0.0
        self.c_root = 1.0e-9 * nb * nb        # core-seconds per square root until one is timed
        self._nxt = None
        self.act, self.est = {}, {}           # core-seconds by B1: measured, modelled
        self.nsplit = 0
        self.split_b1 = {}                    # candidates that lost a factor to a unit, by B1
        self.prefetched = True                # the first chunk starts now, while the table is built
        self.build_ref = 0.0                  # estimated build of a typical candidate, seconds (see should_build)
        self.build_refs = []
        self.pass_t0, self.pass_n = None, 0
        self._launch_scan(threads)

    def now(self):
        return time.time()

    def log(self, msg):
        log(msg)

    # ---- supply
    def _density(self, B, cap):
        """Share of the discriminants below B whose primes are all below cap."""
        import shortsim as ss
        return float(ss.rho(max(1.0, math.log(B) / math.log(max(cap, 3)))))

    def next_chunk(self):
        """(orders expected, core-seconds expected) of the cheaper of the two ways
        to more orders, or None when neither is left:
          widen  the discriminant bound B -> grow * B at the present cap.  Counts
                 go like sqrt(B) (the chance that a discriminant is solvable
                 falls like one over its class number) and work like B.
          raise  the cap on the primes of a discriminant, cap -> 4 * cap, over the
                 range already scanned: the discriminants passed over so far for
                 a prime above the cap.  They cost what discriminants of their
                 size cost, plus a square root mod p for every prime added to
                 the factor base, where a wider B pays ever more per order.
        Both are predicted from what the scans so far have measured, and the
        prediction of a raise is corrected by the outcome of the last one."""
        if self.last is None:
            return (1.0, 1.0)                 # the first chunk is on its way
        if self._nxt is not None:
            return self._nxt[0]
        n, cs, cstrip, lo, hi, dens = self.last
        strip1 = cstrip / max(n, 1)
        opts = []
        if self.B < self.Bmax:
            Bn = min(self.Bmax, int(self.B * self.grow))
            f = self._density(Bn, self.cap) / dens
            nn = max(1.0, n * (math.sqrt(Bn) - math.sqrt(self.B)) / (math.sqrt(hi) - math.sqrt(lo)) * f)
            cc = cs * (Bn - self.B) / (hi - lo) * f + strip1 * nn
            opts.append((cc / nn, nn, cc, ("widen", Bn, self.cap)))
        ncap = min(self.capmax, 4 * self.cap, self.B // 30)
        if not self.fixed_cap and ncap >= 2 * self.cap and self.tot_orders:
            gain = self._density(self.B, ncap) / self._density(self.B, self.cap) - 1.0
            roots = 0.5 * (ncap / (math.log(ncap) - 1.0) - self.cap / (math.log(self.cap) - 1.0))
            nn = max(1.0, self.tot_orders * gain * self.raise_fit[0])
            cc = self.tot_test_s * gain * self.raise_fit[1] + roots * self.c_root + strip1 * nn
            opts.append((cc / nn, nn, cc, ("raise", self.B, ncap)))
        if not opts:
            return None
        _, nn, cc, arm = min(opts)
        self._nxt = ((nn, cc), arm, [(o[3][0], round(o[0], 4)) for o in opts])
        self.arm = arm
        return self._nxt[0]

    def _launch_scan(self, ns):
        ep = self.epoch
        kind, hi, cap = self.arm
        lo, cap0 = self.B, self.cap
        pred = self._nxt

        def scan():
            try:
                t0 = time.time()
                info = {}
                if kind == "raise":           # Bmin = B: only what the larger cap adds
                    recs = run_dscan(self.p, hi, hi, ns, cap, cap0, self.roots, reg=self.procs, info=info)
                else:
                    recs = run_dscan(self.p, hi, lo, ns, cap, cap0 if (cap > cap0 and lo) else 0, self.roots,
                                     reg=self.procs, info=info)
                t1 = time.time()
                Ns = sorted({N for (_, t, _) in recs for N in (self.p + 1 - t, self.p + 1 + t) if N % 4 == 0})
                k = min(ns, STRIP_THREADS)
                if k < ns:
                    self.q.put((ep, ("scanthreads", k)))
                strip = smooth_strip(Ns, k, self.nsm, reg=self.procs)
                t2 = time.time()
                self.q.put((ep, ("orders", (lo, hi, cap, recs, strip, (t1 - t0) * ns, (t2 - t1) * k, ns,
                                            kind, cap0, info, pred))))
            except Exception as e:            # killed by abort(), or a real failure
                self.q.put((ep, ("scanfail", repr(e))))
        threading.Thread(target=scan, daemon=True).start()

    def start_scan(self, ns):
        if self.prefetched:
            self.prefetched = False
            return
        self._launch_scan(ns)

    def take_orders(self, payload):
        from shortsched import sched_init
        lo, hi, fbc, recs, strip, cs, cstrip, ns, kind, cap0, info, pred = payload
        p, L = self.p, self.L
        fresh = []
        for (d, t, v) in recs:
            for N in (p + 1 - t, p + 1 + t):
                if N % 4 == 0:
                    fresh.append(Cand(d, t, v, N))
        live = []
        for c in fresh:
            S, c.tail, c.s, prp = strip[c.N]
            c.mid = S // c.s
            if c.tail == 1 or c.s * c.tail <= L:
                c.dead = True
            self.settle_tail(c, prp)
            if c.dead:
                continue
            # coverage: the divisors of s cover every delta in (0, log2 s) unless
            # the revised format's radical cap cuts in (s above 2^RADLIM only)
            W = math.log2(c.s) if c.s > 1 else 0.0
            if V2_RADLIM and c.s.bit_length() > V2_RADLIM:
                import shortsim as ss
                W = float(ss.coverage_cells(factor_small(c.s, self.small_primes), bcap=V2_BCAP or None,
                                            radlim=V2_RADLIM).sum())
            sched_init(c, int(round(W)), int(round(math.log2(c.tail))))
            live.append(c)
        if live and self.nb >= 500:
            # what a typical winner of this pool costs to build: the median
            # estimate over a sample of the candidates bought so far
            sample = sorted({c.d for c in live})
            sample = sample[::max(1, len(sample) // 48)][:64]
            self.build_refs += list(build_estimate(sample, self.nb).values())
            if self.build_refs:
                self.build_ref = sorted(self.build_refs)[len(self.build_refs) // 2]
        nf = max(len(fresh), 1)
        test_s = info.get("scan_s", 0.0) * info.get("threads", ns) if info else cs
        if info and info.get("roots_new", 0) > 1000 and info.get("root_s", 0.0) > 0.5:
            self.c_root = info["root_s"] * info["threads"] / info["roots_new"]
        if kind == "raise":
            if pred is not None:              # how the prediction of this raise did
                (pn, pc), _, _ = pred
                self.raise_fit = (min(1.5, max(0.2, self.raise_fit[0] * len(fresh) / max(pn, 1.0))),
                                  min(3.0, max(0.3, self.raise_fit[1] * cs / max(pc - cstrip, 1.0))))
        else:
            self.last = (nf, cs, cstrip, lo, hi, self._density(hi, fbc))
        self.B, self.cap = hi, fbc
        self.tot_orders += len(fresh)
        self.tot_test_s += test_s
        self._nxt = None
        self.stats["cands"] = self.stats.get("cands", 0) + len(fresh)
        self.stats["scan_s"] = self.stats.get("scan_s", 0.0) + cs + cstrip
        if kind == "raise":
            log("  level %d bits: dscan[0,%g) cap %g -> %g on %d threads: +%d D -> +%d orders (%d live), "
                "%.0f core-s (strip %.0f) = %.4f per order; %d winners"
                % (self.nb, hi, cap0, fbc, ns, len(recs), len(fresh), len(live), cs + cstrip, cstrip,
                   (cs + cstrip) / nf, len(self.winners)))
        else:
            log("  level %d bits: dscan[%g,%g) maxfb=%g on %d threads: +%d D -> +%d orders (%d live), "
                "%.0f core-s (strip %.0f) = %.4f per order; %d winners"
                % (self.nb, lo, hi, fbc, ns, len(recs), len(fresh), len(live), cs + cstrip, cstrip,
                   (cs + cstrip) / nf, len(self.winners)))
        nx = self.next_chunk()
        if nx is not None and self._nxt is not None and len(self._nxt[2]) > 1:
            log("    next orders: %s (core-s per order: %s)"
                % (self._nxt[1][0], ", ".join("%s %.3f" % o for o in self._nxt[2])))
        return live, len(fresh), cs + cstrip

    # ---- work
    def start_unit(self, u):
        tails = [c.tail for c in u.cands]
        ep, runner = self.epoch, self.runner
        tmo = ecm_timeout(u.b1, u.n, self.nb) * len(tails)

        def work():
            t = time.time()
            try:
                res = runner.run(tails, u.b1, u.n, tmo, u.method)
            except Exception as e:
                log("    ECM unit failed: %r" % (e,))
                res = {}
            if res is not None:
                self.q.put((ep, ("ecm", u, (tails, res, time.time() - t))))
        self.ex.submit(work)

    def wait(self, timeout):
        out = []
        try:
            item = self.q.get(timeout=timeout)
        except queue.Empty:
            return out
        while True:
            ep, ev = item
            if ep == self.epoch:
                if ev[0] == "scanfail":
                    raise RuntimeError("scan failed: %s" % ev[1])
                out.append(ev)
            try:
                item = self.q.get_nowait()
            except queue.Empty:
                return out

    def apply(self, u, payload):
        tails, res, dt = payload
        self.act[u.b1] = self.act.get(u.b1, 0.0) + dt
        self.est[u.b1] = self.est.get(u.b1, 0.0) + u.est
        for c, t in zip(u.cands, tails):
            if c.dead or t not in res:
                continue
            found, cof, cprp = res[t]
            self.nsplit += 1
            self.split_b1[u.b1] = self.split_b1.get(u.b1, 0) + 1
            self.process_result(c, found, cof)
            if c.tail == cof:                         # the ECM program already tested this cofactor
                self.settle_tail(c, cprp)
            elif not c.dead and c.tail > 1:           # the unit ran on a multiple of the present tail
                self.settle_tail(c, pow(2, c.tail - 1, c.tail) == 1)
            if not c.dead and c.tail > 1:
                c.tb = int(round(math.log2(c.tail)))

    def has_winner(self):
        return bool(self.winners)

    def should_build(self, R):
        """Whether to stop the search and build the cheapest winner in hand.  A
        winner whose class number has a large prime factor takes one to five
        minutes to build where the next winner is often seconds away: it is
        kept in reserve while its estimated build exceeds the expected wait for
        another winner (R core-seconds on all the threads) plus a typical
        build, and for no longer than that build would have taken."""
        ws = [w[1].d for w in self.winners]
        new = [d for d in ws if d not in BUILD_EST]
        if new:
            BUILD_EST.update(build_estimate(new, self.nb))
        tb = min(BUILD_EST.get(d, 0.0) for d in ws)
        now = time.time()
        wait = (R / self.threads if 0 < R < float("inf") else 0.0) + self.build_ref
        if tb <= max(wait, 5.0) or (self.pass_t0 is not None and now - self.pass_t0 >= tb):
            return True
        if self.pass_t0 is None:
            self.pass_t0 = now
        if len(ws) > self.pass_n:
            self.pass_n = len(ws)
            log("    %d winner(s) kept in reserve: the cheapest build is estimated at %.0f s, another winner and its build at %.0f s"
                % (len(ws), tb, wait))
        return False

    def abort(self):
        self.epoch += 1
        self.runner.abort()
        self.runner = EcmRunner(self.engine)
        for pr in list(self.procs):
            try:
                pr.kill()
            except OSError:
                pass

    def build(self, st):
        self.stats["ecm_s"] = self.stats.get("ecm_s", 0.0) + st["ecm_s"] - self.ecm_seen
        self.ecm_seen = st["ecm_s"]
        ratios = " ".join("%g:%.2f" % (b1, self.act[b1] / self.est[b1]) for b1 in sorted(self.act) if self.est[b1] > 0)
        eb = st.get("expwin_b1") or {}
        es = st.get("expsplit_b1") or {}
        log("    search: %d orders, scan %.0f + ecm %.0f core-s (%.0f above the rule), %.1fs wall, %d units, "
            "%d splits, rule R = %.0f core-s; the model expected %.2f wins from this work; "
            "measured/modelled ECM cost by B1: %s; wins expected by B1: %s; won at B1=%s; "
            "splits found/expected by B1: %s"
            % (st["orders"], st["scan_s"], st["ecm_s"], st.get("over_s", 0.0), st["wall"], st["units"],
               self.nsplit, st["R"], st["expwin"], ratios,
               " ".join("%g:%.3f" % (b1, eb[b1]) for b1 in sorted(eb)), st.get("win_b1"),
               " ".join("%g:%d/%.0f" % (b1, self.split_b1.get(b1, 0), es[b1]) for b1 in sorted(es))))
        self.pass_t0, self.pass_n = None, 0
        return self.try_build()

    def close(self):
        self.abort()
        self.ex.shutdown(wait=False)


def run_rule_level(p, nb, L, n2, small_primes, threads, stats, B0, Bmax_given, fixed_cap, roots,
                   winners, process_result, settle_tail, try_build):
    import shortsim as ss
    import shortsched as sd
    lam = math.log2(L)
    engine = ecm_engine()
    f = ecm_machine_factor(nb, threads, engine)
    # the scan is bounded by the rule (its price), not by a budget: dscan's
    # products stay within 64 bits up to 2^60, and long before that a wider scan
    # costs more per order than raising the cap on the primes (RuleDriver.next_chunk)
    Bmax = Bmax_given if Bmax_given else 10 ** 15
    nsm = V2_BCAP if V2_BCAP else n2          # fillers use primes up to B only in the revised format
    drv = RuleDriver(p, nb, L, threads, stats, B0, Bmax, fixed_cap, roots, nsm, small_primes,
                     winners, process_result, settle_tail, try_build, engine=engine)
    try:
        t0 = time.time()
        ss.set_grid(nb)
        em, sec, mcost, quantum = ss.engine_model(engine, f, ECPP_DIR)
        cost = ss.engine_model(engine, f, ECPP_DIR, scalar=True)[2]
        tiers = sd.default_tiers(engine, nb)
        lad = ss.Ladder(tiers, em, y0bits=Y0.bit_length() - 1, cost=mcost)
        # rewards of interest shrink fast with the size (a win costs 1.5 core-h
        # at 768 bits, a minute at 520 and a second at 300)
        rlo = min(30.0, 3.0 * (nb / 500.0) ** 8)
        tab = ss.IndexTable(lad, lam, rlo=rlo, rhi=3.0e7 if nb >= 500 else rlo * 1e5, threads=min(16, threads))
        log("  level %d bits: rule table in %.1fs (%d sweeps, rewards %.0f..%.0f core-s)"
            % (nb, time.time() - t0, len(tiers), tab.rgrid[0], tab.rgrid[-1]))
        return sd.level_loop(drv, tab, tiers, cost, {"quantum": quantum})
    finally:
        drv.close()



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
    cmmin = cm_min_bits(V2_QFLOOR if v2 else n2)
    while p.bit_length() >= cmmin:
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
            % (cmmin, "short2.gp" if v2 else "short.gp", p.bit_length()))
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
    if "sched" in args:
        global SCHED
        SCHED = args["sched"]
    global ENGINE, YAFU_BIN
    ENGINE = args.get("engine", ENGINE)
    YAFU_BIN = args.get("yafu", YAFU_BIN)
    if ENGINE not in ("auto", "gmp", "yafu"):
        sys.exit("engine must be auto, gmp or yafu")
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
