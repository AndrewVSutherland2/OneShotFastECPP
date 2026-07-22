#!/usr/bin/env python3
"""chain_prove.py -- descent-chain ECPP prover (prototype).

Builds a chain of Goldwasser-Kilian certificates for prime p.  At each level:
  1. dscan p=<p> dump enumerates (D,t,v) with 4p = t^2 + D v^2; candidate curve
     orders are N = p+1-t and p+1+t.
  2. Each N is split as N = c * q with q prime in ((p^(1/4)+1)^2, 2^(b*bits(p))]
     by peeling the cofactor: trial division (PARI factor(N,10^6)) then rounds
     of ECM with escalating B1, breadth-first across the candidate pool.
  3. The winning (D,t) is turned into a curve via polclass(-D) mod p; the twist
     with #E = N is identified by testing Q = [N/q]P on both twists; the step
     certificate is (p, a, b, Q, q) with [q]Q = O.
Descend to q; stop once q < 2^64 (deterministic Miller-Rabin range).

usage: chain_prove.py p=<decimal> [b=0.85] [B=262144] [threads=<n>]
                      [out=<file.json>] [seed=1] [ecm=<path-to-gmp-ecm>]

The certificate is verified with vchain.py before it is written.
"""

import concurrent.futures
import json
import math
import os
import random
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import vchain

ECPP_DIR = os.path.dirname(os.path.abspath(__file__))
TRIAL_LIMIT = 10 ** 6
# (B1, curves) escalation ladder: roughly the 15..35-digit ECM levels at
# fractional effort; breadth-first so cheap winners surface before deep ECM.
ECM_ROUNDS = [(2000, 16), (11000, 24), (50000, 40), (250000, 60), (1000000, 90)]


def log(msg):
    print(msg, file=sys.stderr, flush=True)


def is_prp(n):
    if n < 2:
        return False
    for b in vchain.MR_BASES:
        if n % b == 0:
            return n == b
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    bases = list(vchain.MR_BASES)
    if n >= 1 << 64:
        bases += [random.randrange(2, n - 2) for _ in range(8)]
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
    # Tonelli-Shanks; p an odd prime, a a QR mod p
    a %= p
    if a == 0:
        return 0
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


def gp(script, timeout=600):
    inp = 'default(parisizemax,"8G");\n' + script
    r = subprocess.run(["gp", "-q"], input=inp, capture_output=True,
                       text=True, timeout=timeout)
    if r.returncode != 0:
        raise RuntimeError("gp failed: %s" % r.stderr[:500])
    return r.stdout


def run_dscan(p, B, Bmin, threads):
    cmd = [os.path.join(ECPP_DIR, "dscan"), "p=%d" % p, "B=%d" % B,
           "threads=%d" % threads, "dump"]
    if Bmin:
        cmd.append("Bmin=%d" % Bmin)
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError("dscan failed: %s" % r.stderr[:500])
    out = []
    for line in r.stdout.splitlines():
        w = line.split()
        if len(w) == 3:
            d, t, v = int(w[0]), int(w[1]), int(w[2])
            if d > 4:  # skip j=0/1728 twist zoo
                out.append((d, t, v))
    out.sort()
    return out


def classify_batch(tails):
    """PARI trial division to 10^6: tails[i] -> (unfactored part R, R is prp)."""
    if not tails:
        return []
    lines = ["default(primelimit,1100000);",
             "cls(N)={my(F=factor(N,10^6),c=1);"
             "for(i=1,matsize(F)[1],if(F[i,1]<10^6,c*=F[i,1]^F[i,2]));"
             "my(R=N/c);print(R,\" \",ispseudoprime(R))};"]
    lines += ["cls(%d);" % n for n in tails]
    out = gp("\n".join(lines) + "\n")
    res = []
    for line in out.splitlines():
        w = line.split()
        if len(w) == 2:
            res.append((int(w[0]), w[1] == "1"))
    if len(res) != len(tails):
        raise RuntimeError("classify_batch: got %d results for %d inputs"
                           % (len(res), len(tails)))
    return res


class Candidate:
    __slots__ = ("d", "t", "v", "N", "tail", "dead")

    def __init__(self, d, t, v, N):
        self.d, self.t, self.v, self.N = d, t, v, N
        self.tail = None
        self.dead = False


def ecm_peel(ecm_bin, n, b1, curves, timeout):
    """Run ECM; gmp-ecm keeps going on the cofactor, so one invocation can
    peel several factors.  Output is 'f1 f2 ... cofactor' (or just n if no
    factor found); returns the remaining cofactor."""
    try:
        r = subprocess.run([ecm_bin, "-q", "-c", str(curves), str(b1)],
                           input=str(n) + "\n", capture_output=True, text=True,
                           timeout=timeout)
    except subprocess.TimeoutExpired:
        return n
    try:
        toks = [int(x) for x in r.stdout.split()]
    except ValueError:
        return n
    if toks and math.prod(toks) == n:
        return toks[-1]
    return n


def prove_level(p, b_max, B0, threads, ecm_bin, stats):
    nbits = p.bit_length()
    q_hi = 1 << int(b_max * nbits)
    q_lo = (vchain.introot(p, 4) + 2) ** 2
    t0 = time.time()
    B, Bmin = B0, 0
    winners = []      # (q, cand)
    queue = []        # composite tails still being peeled
    n_ecm = 0

    def classify(cands):
        res = classify_batch([c.tail for c in cands])
        for c, (R, prp) in zip(cands, res):
            c.tail = R
            if R < q_lo:
                c.dead = True
            elif prp:
                c.dead = True
                if R <= q_hi:
                    winners.append((R, c))
            # else: stays live (composite tail >= q_lo)

    pool_n = 0

    def try_build():
        # attempt winners smallest-q first; drop any that fail construction
        winners.sort(key=lambda w: w[0])
        while winners:
            q, c = winners.pop(0)
            step = build_step(p, c, q)
            if step:
                step["aux"] = {"D": -c.d, "t": c.t, "v": c.v,
                               "sign": 1 if c.N == p + 1 - c.t else -1,
                               "c": str(c.N // q), "pool": pool_n,
                               "ecm_calls": n_ecm,
                               "seconds": round(time.time() - t0, 2)}
                log("  level %d bits done in %.1fs: D=%d, q %d bits (drop %.3f)"
                    % (nbits, time.time() - t0, -c.d, q.bit_length(),
                       q.bit_length() / nbits))
                return step
            log("  winner (D=%d) failed curve build, trying next" % -c.d)
        return None

    while True:
        if Bmin < B:
            cands = [Candidate(d, t, v, N)
                     for (d, t, v) in run_dscan(p, B, Bmin, threads)
                     for N in (p + 1 - t, p + 1 + t)]
            Bmin = B
            for c in cands:
                c.tail = c.N
            classify(cands)
            queue += [c for c in cands if not c.dead]
            pool_n += len(cands)
            stats["candidates"] = stats.get("candidates", 0) + len(cands)
            log("  level %d bits: pool +%d (B=%d), %d live, %d winners after trial division"
                % (nbits, len(cands), B, len(queue), len(winners)))
        for b1, curves in ECM_ROUNDS:
            if winners:
                break
            timeout = 60 + curves * b1 // 100000
            live = [c for c in queue if not c.dead]
            with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as ex:
                found = list(ex.map(
                    lambda c: ecm_peel(ecm_bin, c.tail, b1, curves, timeout), live))
            n_ecm += len(live)
            stats["ecm_calls"] = stats.get("ecm_calls", 0) + len(live)
            n_split = 0
            for c, tail in zip(live, found):
                if tail == c.tail:
                    continue
                n_split += 1
                c.tail = tail
                if tail < q_lo:
                    c.dead = True
                elif is_prp(tail):
                    c.dead = True
                    if tail <= q_hi:
                        winners.append((tail, c))
            queue = [c for c in queue if not c.dead]
            log("    ecm round B1=%d c=%d: %d tested, %d split, %d live, %d winners"
                % (b1, curves, len(live), n_split, len(queue), len(winners)))
        if winners:
            step = try_build()
            if step:
                return step
            continue  # all current winners failed construction: keep searching
        Bmin, B = B, B * 4
        log("  no winner: widening scan to B=%d" % B)


def build_step(p, cand, q):
    """polclass(-D) root -> curve; find the twist with #E = N and a point Q of
    order q.  Everything after the class-polynomial root is plain Python."""
    out = gp("if(!ispseudoprime(%d),error(\"p not prime\"));\n"
             "r=polrootsmod(polclass(-%d),%d);\n"
             "if(#r==0,print(\"NOROOT\"),print(lift(r[1])));\n" % (p, cand.d, p))
    w = out.split()
    if not w or w[-1] == "NOROOT":
        return None
    j = int(w[-1])
    if j == 0 or j % p == 1728 % p:
        return None
    k = j * pow((1728 - j) % p, -1, p) % p
    a, b = 3 * k % p, 2 * k % p
    u = 2
    while pow(u, (p - 1) // 2, p) != p - 1:
        u += 1
    c = cand.N // q
    rng = random.Random(cand.d * 1000003 + q % 1000003)
    for twist in range(2):
        if twist:
            a, b = a * u * u % p, b * u * u * u % p
        usable = 0
        for _ in range(256):  # draws; ~half yield a point on this twist
            x = rng.randrange(p)
            z = (x * x * x + a * x + b) % p
            if z == 0 or pow(z, (p - 1) // 2, p) != 1:
                continue
            usable += 1
            P = (x, sqrt_mod(z, p))
            Q = ec_mul(c, P, a, p)
            if Q is not None:
                if ec_mul(q, Q, a, p) is None:
                    return {"p": str(p), "a": str(a), "b": str(b),
                            "Qx": str(Q[0]), "Qy": str(Q[1]), "q": str(q)}
                break  # wrong twist: [c]P has order not divisible by q
            if usable >= 3:
                break  # [c]P = O three times: q does not divide this order
    return None


def ec_mul(k, P, a, p):
    acc = None  # O; dbl(O)=O and O+P=P keep the semantics exact mid-ladder
    for bit in bin(k)[2:]:
        if acc is not None:
            acc = vchain.ec_dbl(acc, a, p)
        if bit == "1":
            acc = P if acc is None else vchain.ec_add(acc, P, a, p)
    return acc


def main():
    args = dict(kv.split("=", 1) for kv in sys.argv[1:] if "=" in kv)
    if "p" not in args:
        print(__doc__)
        sys.exit(2)
    p0 = int(args["p"])
    b_max = float(args.get("b", "0.85"))
    B0 = int(args.get("B", "262144"))
    threads = int(args.get("threads", str(os.cpu_count())))
    ecm_bin = args.get("ecm", os.environ.get("CHAIN_ECM", "ecm"))
    random.seed(int(args.get("seed", "1")))
    if not (0.5 < b_max < 1):
        sys.exit("need 0.5 < b < 1")
    if not is_prp(p0):
        sys.exit("p is not prime")

    t0 = time.time()
    steps, stats, p = [], {}, p0
    while p >= 1 << 64:
        steps.append(prove_level(p, b_max, B0, threads, ecm_bin, stats))
        p = int(steps[-1]["q"])
    cert = {"format": "ecpp-descent-chain-v1", "p": str(p0), "steps": steps}
    vchain.verify(cert, verbose=False)  # never emit an unverified certificate
    out = args.get("out", "chain_%d.json" % p0.bit_length())
    with open(out, "w") as f:
        json.dump(cert, f, indent=1)
    log("chain of %d steps for %d-bit p in %.1fs (%d candidate orders, %d ecm calls)"
        % (len(steps), p0.bit_length(), time.time() - t0,
           stats.get("candidates", 0), stats.get("ecm_calls", 0)))
    log("certificate written to %s (%d bytes), verified by vchain" %
        (out, os.path.getsize(out)))


if __name__ == "__main__":
    main()
