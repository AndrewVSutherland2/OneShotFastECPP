#!/usr/bin/env python3
"""vchain.py -- verify a descent-chain ECPP certificate (JSON).

Certificate format (chain_prove.py output):
  { "format": "ecpp-descent-chain-v1", "p": "<decimal>",
    "steps": [ { "p": "...", "a": "...", "b": "...", "Qx": "...", "Qy": "...",
                 "q": "...", ... }, ... ] }

Each step is a Goldwasser-Kilian style certificate: on E: y^2 = x^3 + ax + b
over Z/p, the point Q = (Qx,Qy) satisfies [q]Q = O with q > (p^(1/4)+1)^2.
If q is prime, every prime divisor r of p has an elliptic curve group over F_r
containing a point of order q, so #E(F_r) >= q and Hasse gives
r >= (sqrt(q)-1)^2 > sqrt(p); hence p is prime.  q's primality comes from the
next step of the chain; the final q < 2^64 is proven prime by Miller-Rabin with
the 12 prime bases 2..37 (deterministic below 3.317e24, Sorenson-Webster 2015).

Soundness with p composite: all curve arithmetic is affine with every inversion
guarded by gcd(den, p).  A gcd strictly between 1 and p exhibits a factor of p
(certificate rejected); gcd 1 means the denominator is invertible mod p, hence
nonzero mod every r | p, so the computation commutes with reduction mod r.  The
additions that produce O do so via congruences mod p (x1=x2, y1=-y2, or y=0),
which also hold mod r.  A valid Q has prime order q, so no intermediate multiple
[k]Q (k a proper prefix of q, or a power of 2 times one) can be O; the ladder
therefore rejects any premature O.

usage: vchain.py <certificate.json> [quiet]
Exit status 0 iff the certificate proves p prime.
"""

import json
import math
import sys
import time

MR_BASES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)  # deterministic < 3.317e24
BASE_BOUND = 1 << 64


class Bad(Exception):
    pass


def is_prime_small(n):
    # deterministic Miller-Rabin for n < 3.317e24 (we use it only below 2^64)
    if n < 2:
        return False
    for b in MR_BASES:
        if n % b == 0:
            return n == b
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for b in MR_BASES:
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


def introot(n, k):
    # floor(n^(1/k))
    if n < 0:
        raise ValueError
    if n < 2:
        return n
    x = 1 << (n.bit_length() // k + 1)
    while True:
        y = ((k - 1) * x + n // x ** (k - 1)) // k
        if y >= x:
            return x
        x = y


def inv_mod(d, p):
    # guarded inverse: gcd in (1,p) proves p composite; gcd p means d = 0 mod p
    d %= p
    g = math.gcd(d, p)
    if g == 1:
        return pow(d, -1, p)
    if g < p:
        raise Bad("p composite: nontrivial gcd %d" % g)
    return None  # d = 0 mod p


def ec_dbl(P, a, p):
    x, y = P
    lam_den = inv_mod(2 * y, p)
    if lam_den is None:
        return None  # y = 0 mod p: 2P = O
    lam = (3 * x * x + a) * lam_den % p
    x3 = (lam * lam - 2 * x) % p
    return (x3, (lam * (x - x3) - y) % p)


def ec_add(P1, P2, a, p):
    x1, y1 = P1
    x2, y2 = P2
    if (x1 - x2) % p == 0:
        if (y1 + y2) % p == 0:
            return None  # inverse points: sum is O
        return ec_dbl(P1, a, p)
    den = inv_mod(x2 - x1, p)
    if den is None:  # x1 = x2 mod p handled above; unreachable
        raise Bad("inconsistent addition")
    lam = (y2 - y1) * den % p
    x3 = (lam * lam - x1 - x2) % p
    return (x3, (lam * (x1 - x3) - y1) % p)


def mul_is_exact_order(q, Q, a, p):
    """True iff the affine ladder for [q]Q reaches O exactly at the last
    operation (never earlier).  q must be odd."""
    bits = bin(q)[2:]
    n_ops = (len(bits) - 1) + (bits.count("1") - 1)
    acc = Q
    opn = 0
    for bit in bits[1:]:
        acc = ec_dbl(acc, a, p)
        opn += 1
        if acc is None:
            return opn == n_ops
        if bit == "1":
            acc = ec_add(acc, Q, a, p)
            opn += 1
            if acc is None:
                return opn == n_ops
    return False  # never reached O


def verify_step(p, a, b, Qx, Qy, q):
    if p < 11 or p % 2 == 0:
        raise Bad("step modulus %d too small or even" % p)
    if q % 2 == 0 or q >= p:
        raise Bad("bad q")
    r4 = introot(p, 4)
    if q < (r4 + 2) ** 2:
        raise Bad("q too small: q <= (p^(1/4)+1)^2")
    if math.gcd((4 * a * a * a + 27 * b * b) % p, p) != 1:
        raise Bad("singular curve (or gcd with p nontrivial)")
    if not (0 <= Qx < p and 0 <= Qy < p):
        raise Bad("point coordinates out of range")
    if (Qy * Qy - (Qx * Qx * Qx + a * Qx + b)) % p != 0:
        raise Bad("Q not on curve")
    if not mul_is_exact_order(q, (Qx, Qy), a, p):
        raise Bad("[q]Q != O (or premature O)")


def verify(cert, verbose=True):
    if cert.get("format") != "ecpp-descent-chain-v1":
        raise Bad("unknown certificate format")
    p0 = int(cert["p"])
    steps = cert["steps"]
    if not steps:
        raise Bad("empty chain")
    if int(steps[0]["p"]) != p0:
        raise Bad("chain does not start at p")
    t0 = time.time()
    for i, s in enumerate(steps):
        p, a, b = int(s["p"]), int(s["a"]), int(s["b"])
        Qx, Qy, q = int(s["Qx"]), int(s["Qy"]), int(s["q"])
        t1 = time.time()
        verify_step(p, a, b, Qx, Qy, q)
        if i + 1 < len(steps):
            if int(steps[i + 1]["p"]) != q:
                raise Bad("chain broken at step %d" % i)
        else:
            if q >= BASE_BOUND:
                raise Bad("final q >= 2^64 but chain ends")
            if not is_prime_small(q):
                raise Bad("final q composite")
        if verbose:
            print("step %2d: p %4d bits -> q %4d bits   %.3f ms"
                  % (i, p.bit_length(), q.bit_length(), 1e3 * (time.time() - t1)))
    if verbose:
        print("chain of %d steps verified in %.3f ms: p is prime" %
              (len(steps), 1e3 * (time.time() - t0)))
    return True


def main():
    args = [a for a in sys.argv[1:] if a != "quiet"]
    quiet = "quiet" in sys.argv[1:]
    if len(args) != 1:
        print(__doc__)
        sys.exit(2)
    with open(args[0]) as f:
        cert = json.load(f)
    try:
        verify(cert, verbose=not quiet)
    except Bad as e:
        print("INVALID certificate: %s" % e)
        sys.exit(1)
    if quiet:
        print("valid")


if __name__ == "__main__":
    main()
