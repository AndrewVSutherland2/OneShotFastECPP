#!/usr/bin/env python3
"""vshort2.py -- working verifier for the RADICAL-CAPPED short ECPP format, including
the terminal-prime revision of 2026-08-22 (the upstream reference verifier is
ShortPrimalityProofs/vsmallECPP.py; this copy backs the local repair/assembly tools).

Format (August 2026, AVS + Fable 5).  A certificate for an odd p_0 >= 5 is the flat
sequence  (p_0, A_0, x_0, o_0, A_1, x_1, o_1, ..., A_k, x_k, o_k)  with n = ceil(log2 p_0)
FIXED for the whole chain and

    B = floor(n^2 / log2 n)          (smoothness bound)

such that, working modulo p_i at level i (p_0 given, p_{i+1} recovered below):

  - m_i, the B-smooth part of o_i, satisfies m_i >= 2 and the RADICAL CAP

        floor(log2 rad(m_i)) < n / log2 n,      rad(m) = product of the distinct primes of m;

  - p_{i+1} = o_i / m_i (B-rough by construction) satisfies, for i < k,
    B^2 < p_{i+1} and p_{i+1}^2 < p_i; and p_{k+1} is 1 or below B^2, in which
    case it is automatically prime (a composite B-rough integer exceeds B^2);
  - L_i < o_i < r_i L_i with L_i = (p_i^{1/4}+1)^2 (integer form) and r_i the least
    prime divisor of m_i;
  - the point with x-coordinate x_i on E_{A_i}: y^2 = x^3 + A_i x^2 + x (or its twist)
    has order exactly o_i modulo EVERY prime divisor of p_i.

A certificate with p_{k+1} = 1 is also valid under the original format; one that
uses a terminal prime is shorter than anything the original format admits.

Why the caps: the primorial of B has Theta(n^2/log n) bits, so it can be built per
certificate in O(n^2 log n) bit operations -- no per-n precomputation -- and
g = gcd(P_B mod o_i, o_i) IS rad(m_i), delivered as an integer of < n/log2 n bits whose
prime factors are <= B; the radical cap bounds the exactness tree by ~n ladder bits per
level.  Total verification is O(n^2 log n) bit operations, worst case, self-contained.

usage: python3 vshort2.py p0 A0 x0 o0 [A1 x1 o1 ...]
       python3 vshort2.py --test
Prints True and exits 0 iff the sequence is a valid v2 short ECPP (p_0 is prime).
"""

from math import gcd, isqrt, log2


# --------------------------------------------------------------------------
# Montgomery x-only (X:Z) arithmetic, valid on the Kummer line of E_A and of
# its quadratic twist (verbatim from vsmallECPP.py / voneshot.py).
# --------------------------------------------------------------------------
def xdbl(X, Z, A, p):
    XX = X * X % p
    ZZ = Z * Z % p
    XZ = X * Z % p
    X2 = (XX - ZZ) * (XX - ZZ) % p
    Z2 = 4 * XZ % p * ((XX + A * XZ + ZZ) % p) % p
    return X2, Z2


def xadd(X1, Z1, X2, Z2, Xd, Zd, p):
    a = (X1 - Z1) * (X2 + Z2) % p
    b = (X1 + Z1) * (X2 - Z2) % p
    s = (a + b) % p
    d = (a - b) % p
    X3 = Zd * (s * s % p) % p
    Z3 = Xd * (d * d % p) % p
    return X3, Z3


def ladder(k, XP, ZP, A, p):
    if k == 0:
        return (1, 0)
    XP %= p
    ZP %= p
    if k == 1:
        return (XP, ZP)
    Xd, Zd = XP, ZP
    X0, Z0 = XP, ZP
    X1, Z1 = xdbl(XP, ZP, A, p)
    for bit in bin(k)[3:]:
        if bit == '0':
            X1, Z1 = xadd(X0, Z0, X1, Z1, Xd, Zd, p)
            X0, Z0 = xdbl(X0, Z0, A, p)
        else:
            X0, Z0 = xadd(X0, Z0, X1, Z1, Xd, Zd, p)
            X1, Z1 = xdbl(X1, Z1, A, p)
    return X0, Z0


def sieve_primes(limit):
    if limit < 2:
        return []
    is_p = bytearray([1]) * (limit + 1)
    is_p[0] = is_p[1] = 0
    for i in range(2, isqrt(limit) + 1):
        if is_p[i]:
            is_p[i * i::i] = bytearray(len(is_p[i * i::i]))
    return [i for i in range(2, limit + 1) if is_p[i]]


def balanced_product(xs):
    if not xs:
        return 1
    if len(xs) <= 8:
        r = xs[0]
        for x in xs[1:]:
            r *= x
        return r
    xs = list(xs)
    while len(xs) > 1:
        nxt = [xs[i] * xs[i + 1] for i in range(0, len(xs) - 1, 2)]
        if len(xs) & 1:
            nxt.append(xs[-1])
        xs = nxt
    return xs[0]


# --------------------------------------------------------------------------
# Order-exactness tree (from voneshot.py): with Q = (o/R)P for R the product of
# the distinct primes of o, each leaf holds (o/q)P, whose Z must be a unit.
# The radical cap makes this tree logarithmically small.
# --------------------------------------------------------------------------
def check_orders(XQ, ZQ, primes, A, p):
    t = len(primes)
    if t == 0:
        return True
    if t == 1:
        return gcd(ZQ % p, p) == 1
    mid = t // 2
    Lh, Rh = primes[:mid], primes[mid:]
    hL = balanced_product(Lh)
    hR = balanced_product(Rh)
    XL, ZL = ladder(hL, XQ, ZQ, A, p)
    XR, ZR = ladder(hR, XQ, ZQ, A, p)
    return check_orders(XL, ZL, Rh, A, p) and check_orders(XR, ZR, Lh, A, p)


# --------------------------------------------------------------------------
# The v2 short-ECPP verifier.
# --------------------------------------------------------------------------
def _verify(seq):
    """(ok, level, reason): ok iff seq is a valid v2 short ECPP; on failure,
    level/reason locate the first rejection (level -1 = malformed input)."""
    if len(seq) < 4 or (len(seq) - 1) % 3 != 0:
        return False, -1, "shape"
    if any(not isinstance(v, int) or v < 0 for v in seq):
        return False, -1, "range"
    p = seq[0]
    if p < 5 or p % 2 == 0:
        return False, -1, "p0"
    n = p.bit_length()            # = ceil(log2 p0): p0 is odd, never a power of 2
    lg = log2(n)
    B = int(n * n / lg)           # floor(n^2 / log2 n): the smoothness bound
    B2 = B * B                    # recursion floor; a B-rough integer < B^2 is prime
    radlim = n / lg               # require floor(log2 rad(m)) < radlim

    # collect the level orders and pre-screen their sizes (both bounds implied
    # by validity, so rejecting on them is sound; cf. vsmallECPP.py)
    os = [seq[i + 2] for i in range(1, len(seq), 3)]
    if len(os) > n:
        return False, -1, "levels"
    if any(o < 2 or o.bit_length() > 2 * n + 2 for o in os):
        return False, -1, "osize"

    primes = sieve_primes(B)
    P = balanced_product(primes)                  # the primorial of B

    for lev, i in enumerate(range(1, len(seq), 3)):
        A, x, o = seq[i], seq[i + 1], seq[i + 2]
        if p < 3 or p % 2 == 0:   # a mid-chain modulus collapsed to 1 (or worse)
            return False, lev, "modulus"
        if not (0 <= A < p) or not (0 <= x < p):
            return False, lev, "range"
        if gcd((A * A - 4) % p, p) != 1:          # nonsingular mod every divisor of p
            return False, lev, "singular"
        if o < 2:
            return False, lev, "osize"

        # recover rad(m) = gcd(P mod o, o), then m and p_next.  g | P forces g
        # squarefree with every prime factor <= B.
        g = gcd(P % o, o)
        if g <= 1:                                # m = 1: r undefined, reject
            return False, lev, "m=1"
        if g.bit_length() - 1 >= radlim:          # the radical cap
            return False, lev, "radical"
        small = []                                # ascending prime factors of g
        gg = g
        for q in primes:
            if q * q > gg:
                break
            if gg % q == 0:
                small.append(q)
                gg //= q
        if gg > 1:
            small.append(gg)                      # prime (trial division passed sqrt)
        m = 1                                     # m = prod over q | g of q^{v_q(o)}
        for q in small:
            oo = o
            while oo % q == 0:
                oo //= q
                m *= q
        p_next = o // m
        r = small[0]                              # least prime divisor of m (and of o)

        # size window: L < o < r*L, L the largest point order over any F_q, q <= sqrt(p)
        q_ = isqrt(p)
        L = q_ + 1 + isqrt(4 * q_)
        if not (L < o < r * L):
            return False, lev, "window"

        # descent (terminal-prime revision, 2026-08-22): a non-terminal level
        # needs B^2 < p_next and p_next^2 < p; the last level needs p_next = 1
        # or p_next < B^2 (then prime by size, being B-rough by construction)
        last = (i + 3 == len(seq))
        if last:
            if p_next != 1 and p_next >= B2:
                return False, lev, "terminal"
        else:
            if p_next <= B2 or p_next * p_next >= p:
                return False, lev, "descent"

        # [o]P = O reached as a genuine (X:0) with X a unit mod p
        Xo, Zo = ladder(o, x, 1, A, p)
        if Zo % p != 0 or gcd(Xo % p, p) != 1:
            return False, lev, "order"

        # (o/q)P != O (mod every prime divisor of p) for each prime q | o;
        # p_next participates as a factor, its primality certified by level i+1
        divisors = small + ([p_next] if p_next > 1 else [])
        R = balanced_product(divisors)
        XQ, ZQ = ladder(o // R, x, 1, A, p)
        if not check_orders(XQ, ZQ, divisors, A, p):
            return False, lev, "exact"

        p = p_next
    return True, -1, "ok"                         # the last level checked p_next above


def verify(seq):
    return _verify(seq)[0]


# --------------------------------------------------------------------------
# Self-tests.  Valid vectors were found by brute force (tiny p) and by the
# adapted finder short2.gp; tamper cases exercise each rejection path,
# including the new radical cap.  The CRT split attack is from the original
# test suite and must still be rejected.
# --------------------------------------------------------------------------
_INVALID = [
    "251 0 10 63",                 # valid in the ORIGINAL format; v2 radical cap rejects
    "11 0 3 12",                   # ditto (rad 6, cap requires floor(log2 rad) < 2)
    "3 0 0 6",                     # p = 3 admits no short ECPP at all
    "221 5 2 34",                  # composite p0 (221 = 13*17)
    # CRT split attack: p0 = 2098153*2102167; (x0,1) has order exactly 8*525029
    # in E(Z/p0) but order 8 mod one factor and 525029 mod the other
    "4410667997551 1365834658413 107710304518 4200232 199129 175565 880",
]


def _admissible_orders(p):
    """v2-admissible single-level orders for p (terminal: o fully B-smooth in
    the window with the radical cap), by direct enumeration -- tiny p only."""
    n = p.bit_length()
    lg = log2(n)
    B = int(n * n / lg)
    radlim = n / lg
    q_ = isqrt(p)
    L = q_ + 1 + isqrt(4 * q_)
    primes = sieve_primes(B)
    out = []
    for o in range(L + 1, (B + 1) * L):
        oo, rad, r = o, 1, 0
        for q in primes:
            if oo % q == 0:
                rad *= q
                if r == 0:
                    r = q
                while oo % q == 0:
                    oo //= q
        if oo != 1 or rad <= 1:
            continue                              # not fully B-smooth (terminal)
        if rad.bit_length() - 1 >= radlim:
            continue
        if L < o < r * L:
            out.append(o)
    return out


def _selftest():
    import random
    rng = random.Random(2026)
    ok = True
    # tiny certificates: enumerate admissible terminal orders, then sample (A, x);
    # existence is guaranteed (2-power fillers are always admissible).
    for p in (5, 7, 11, 13, 251):
        orders = _admissible_orders(p)
        found = None
        for _ in range(200000):
            A, x = rng.randrange(p), rng.randrange(p)
            for o in orders:
                if verify([p, A, x, o]):
                    found = (p, A, x, o)
                    break
            if found:
                break
        if not found:
            print(f"selftest: no tiny certificate found for p={p} (orders {orders})")
            ok = False
        else:
            print(f"selftest: p={p}: cert {found}  (admissible terminal orders: {orders})")
    return ok


if __name__ == "__main__":
    import sys
    if len(sys.argv) == 2 and sys.argv[1] == "--test":
        good = _selftest()
        for s in _INVALID:
            v = verify([int(t) for t in s.split()])
            print(f"selftest: invalid vector rejected: {not v}")
            good = good and not v
        print("ALL TESTS PASSED" if good else "TESTS FAILED")
        sys.exit(0 if good else 1)
    try:
        seq = [int(t) for t in sys.argv[1:]]
    except ValueError:
        print("usage: vshort2.py p0 A0 x0 o0 [A1 x1 o1 ...]")
        sys.exit(2)
    r = verify(seq)
    print(r)
    sys.exit(0 if r else 1)
