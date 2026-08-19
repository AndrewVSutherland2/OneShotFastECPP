#!/usr/bin/env python3
"""
Verification of a short ECPP (github.com/AndrewVSutherland/ShortPrimalityProofs)
in quasi-quadratic time.  Montgomery-ladder core follows voneshot.py (Opus 4.8 /
A.V. Sutherland); chain logic per the ShortPrimalityProofs definition.

A short ECPP for p_0 is a flat integer sequence

    (p_0, A_0, x_0, o_0, A_1, x_1, o_1, ..., A_k, x_k, o_k)

with o_i = m_i * p_{i+1} and p_{k+1} = 1, in which, writing n = ceil(log2 p_0)
(fixed at the top level for the whole chain):

  - each p_i is odd (p_0 given; later moduli are recovered from the previous
    level as p_{i+1} = the n^2-rough part of o_i), with p_{i+1}^2 < p_i;
  - m_i, the n^2-smooth part of o_i, has least prime divisor r_i and satisfies
        L_i < o_i < r_i * L_i,
    where q_i = isqrt(p_i) and L_i = q_i + 1 + isqrt(4*q_i) is the largest
    possible order of a point on an elliptic curve over F_q for any q <= sqrt(p_i);
  - 0 <= A_i < p_i with gcd(A_i^2 - 4, p_i) = 1, and 0 <= x_i < p_i;
  - x_i is the x-coordinate of a point of order exactly o_i on the Montgomery
    curve B y^2 = x^3 + A_i x^2 + x over Z/p_i for some B -- equivalently on
    E_{A_i} or its quadratic twist; the x-only ladder never needs B or y.

Why this proves p_0 prime: by induction from the top.  At level i the verifier
establishes ord(P_i) = o_i modulo every prime divisor l of p_i, with every prime
factor of o_i certified (primes <= n^2 by trial division, p_{i+1} by the next
level, 1 trivially).  If some prime l <= sqrt(p_i) divided p_i, then o_i =
ord(P_i mod l) <= #E(F_l) <= l + 1 + floor(2*sqrt(l)) <= L_i < o_i, a
contradiction; hence p_i is prime.  The minimality window o_i < r_i * L_i keeps
the certificate size O(n) (sum of the level sizes is geometric).

Soundness for composite input is inherited from the voneshot.py machinery: the
x-only ladder is valid on the Kummer line of E and its twist so no on-curve test
is needed; [o]P = O must be reached as a genuine (X:0) with gcd(X, p) = 1
(rejecting the degenerate (0:0) collapse that arises from multiplying past the
order modulo a factor of p); each (o/q)P != O leaf requires gcd(Z, p) = 1, i.e.
nonvanishing modulo every prime divisor of p.

Cost (FFT integer multiplication assumed): one sieve to n^2 (O(n^2 loglog n) bit
operations); ONE batched trial-division pass for the whole chain -- a prime can
divide some o_i only if it divides D = prod_i o_i, an O(n)-bit integer, so a
single sweep of remainder trees against D, with all products of primes built by
balanced pairing, finds every level's small prime divisors in O(n^2 (log n)^2)
bit operations including all product building; the elliptic work is O(log o_i)
group operations per tree level with O(log log o_i) levels in arithmetic mod
p_i, and level sizes decay geometrically (p_{i+1} < sqrt(p_i)).  Total:
O(n^2 (log n)^2) bit operations and O(n^2) bits of memory, with no precomputed
tables.  (Two changes versus the previous revision, neither affecting the
accept/reject behavior: products of primes are built by balanced pairing rather
than sequential accumulation, whose cost is quadratic in the product size; and
the trial division is batched across the whole chain rather than repeated per
level, which cost an extra factor of log n.)

usage: python3 vsmallECPP.py p0 A0 x0 o0 [A1 x1 o1 ...]
       python3 vsmallECPP.py --test
Prints True and exits 0 iff the sequence is a valid short ECPP (p_0 is prime).
"""

from math import gcd, isqrt


# --------------------------------------------------------------------------
# Montgomery x-only (X:Z) arithmetic, valid on the Kummer line of E_A and of
# its quadratic twist (verbatim from voneshot.py).
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


# --------------------------------------------------------------------------
# Trial division to n^2 via sieve + batched remainder trees, in one pass for
# the whole chain.
# --------------------------------------------------------------------------
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
    """prod(xs) by power-of-two pairing: O(M(S) log k) bit operations for k
    factors of total size S, versus Theta(S^2/k) for sequential accumulation
    (each step re-scans the growing partial product)."""
    if not xs:
        return 1
    if len(xs) <= 8:                    # bounded count: sequential is O(M(S))
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


def remainder_tree(x, mods):
    if not mods:
        return []
    k = len(mods)
    size = 1
    while size < k:
        size <<= 1
    levels = [list(mods) + [1] * (size - k)]
    while len(levels[-1]) > 1:
        cur = levels[-1]
        levels.append([cur[i] * cur[i + 1] for i in range(0, len(cur), 2)])
    rems = [x % levels[-1][0]]
    for lvl in range(len(levels) - 2, -1, -1):
        cur = levels[lvl]
        rems = [rems[i >> 1] % cur[i] for i in range(len(cur))]
    return rems[:k]


def batch_prime_divisors(os, primes):
    """For each o in os, the ascending list of its distinct prime divisors that
    lie in `primes` -- computed in ONE pass over the prime list.  A prime can
    divide some o only if it divides D = prod(os), so the pass runs remainder
    trees against D (an O(n)-bit integer for a valid chain), with the primes
    grouped so each batch product is comparable to D in size; the few primes
    dividing D (at most log2 D of them) are then distributed to the o's by
    direct division.  One batch is resident at a time, so memory stays
    O(sieve) = O(n^2) bits."""
    D = balanced_product(os)
    batch_bits = max(64, D.bit_length())
    hits = []                                  # ascending primes dividing D
    batch = []
    bits = 0
    for q in primes:
        batch.append(q)
        bits += q.bit_length()
        if bits >= batch_bits:
            Q = balanced_product(batch)
            hits += [t for t, r in zip(batch, remainder_tree(D % Q, batch)) if r == 0]
            batch, bits = [], 0
    if batch:
        Q = balanced_product(batch)
        hits += [t for t, r in zip(batch, remainder_tree(D % Q, batch)) if r == 0]
    return [[q for q in hits if o % q == 0] for o in os]


# --------------------------------------------------------------------------
# Order-exactness tree (from voneshot.py): with Q = (o/R)P for R the product of
# the distinct primes of o, each leaf holds (o/q)P, whose Z must be a unit.
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
# The short-ECPP verifier.
# --------------------------------------------------------------------------
def verify(seq):
    """True iff seq = (p0, A0, x0, o0, ..., Ak, xk, ok) is a valid short ECPP."""
    if len(seq) < 4 or (len(seq) - 1) % 3 != 0:
        return False
    if any(not isinstance(v, int) or v < 0 for v in seq):
        return False
    p = seq[0]
    if p < 5 or p % 2 == 0:
        return False
    n = p.bit_length()            # = ceil(log2 p0): p0 is odd, never a power of 2

    # collect the level orders and pre-screen their sizes so that the batched
    # trial division below runs on a certificate-independent budget.  Both
    # bounds are implied by validity, so rejecting on them is sound: a valid
    # chain has at most 1 + log2(n) levels (the moduli at least halve in bit
    # length and stay above n^2 until the last), and every level order is below
    # the Hasse bound of its modulus, hence below p0^2.
    os = [seq[i + 2] for i in range(1, len(seq), 3)]
    if len(os) > n:
        return False
    if any(o < 2 or o.bit_length() > 2 * n + 2 for o in os):
        return False
    primes = sieve_primes(n * n)
    small_lists = batch_prime_divisors(os, primes)   # one pass for the whole chain

    for lev, i in enumerate(range(1, len(seq), 3)):
        A, x, o = seq[i], seq[i + 1], seq[i + 2]
        if p < 3 or p % 2 == 0:   # a mid-chain modulus collapsed to 1 (or worse)
            return False
        if not (0 <= A < p) or not (0 <= x < p):
            return False
        if gcd((A * A - 4) % p, p) != 1:      # nonsingular mod every divisor of p
            return False
        if o < 2:
            return False

        # split o into its n^2-smooth part m (with distinct primes `small`) and
        # its n^2-rough part p_next, the next modulus (1 at the end of the chain)
        small = small_lists[lev]
        rest = o
        for q in small:
            while rest % q == 0:
                rest //= q
        p_next = rest
        if o == p_next:                       # m = 1: r_i undefined, reject
            return False
        r = small[0]                          # least prime divisor of m (and of o)

        # size window: L < o < r*L, L the largest point order over any F_q, q <= sqrt(p)
        q_ = isqrt(p)
        L = q_ + 1 + isqrt(4 * q_)
        if not (L < o < r * L):
            return False

        # descent: the next modulus must satisfy p_next^2 < p
        if p_next != 1 and p_next * p_next >= p:
            return False

        # [o]P = O reached as a genuine (X:0) with X a unit mod p
        Xo, Zo = ladder(o, x, 1, A, p)
        if Zo % p != 0 or gcd(Xo % p, p) != 1:
            return False

        # (o/q)P != O (mod every prime divisor of p) for each prime q | o;
        # p_next participates as a factor, its primality certified by level i+1
        divisors = small + ([p_next] if p_next > 1 else [])
        R = balanced_product(divisors)
        XQ, ZQ = ladder(o // R, x, 1, A, p)
        if not check_orders(XQ, ZQ, divisors, A, p):
            return False

        p = p_next
    return p == 1                             # the chain must terminate exactly


# --------------------------------------------------------------------------
# Self-tests: valid vectors produced by gen_short_p256.py / the chain builder,
# plus tamper cases exercising each rejection path.
# --------------------------------------------------------------------------
_VALID = [
    # single-level certificates (k = 0), from the exhaustive p <= 256 list
    "5 1 2 8",
    "7 0 2 8",
    "11 0 3 12",
    "13 0 2 10",
    "251 0 10 63",
    # a two-level chain built by gen_short.py chain pbits=48 seed=1
    "251444687128489 109722824127413 239294303725419 18088453 37 8113 204",
]

_INVALID = [
    "251 0 10 126",                # m doubled: breaks the minimality window
    "251 0 11 63",                 # wrong point (x-order differs)
    "251 0 10 63 0 10 63",         # trailing level after the chain terminated
    "221 5 2 34",                  # composite p0 (221 = 13*17)
    "251 2 10 63",                 # singular curve (A = 2)
    "3 0 0 6",                     # p = 3 admits no short ECPP at all
    # CRT split attack: p0 = 2098153*2102167; (x0,1) has order exactly 8*525029
    # in E(Z/p0) but order 8 mod one factor and 525029 mod the other
    "4410667997551 1365834658413 107710304518 4200232 199129 175565 880",
]


def _selftest():
    for line in _VALID:
        assert verify([int(t) for t in line.split()]), line
    for line in _INVALID:
        assert not verify([int(t) for t in line.split()]), line
    print("ok")


if __name__ == "__main__":
    import sys
    args = sys.argv[1:]
    if args == ["--test"]:
        _selftest()
        sys.exit(0)
    if len(args) < 4:
        sys.stderr.write(__doc__)
        sys.exit(2)
    try:
        seq = [int(a, 0) for a in args]
    except ValueError:
        sys.stderr.write("error: all arguments must be integers\n")
        sys.exit(2)
    result = verify(seq)
    print(result)
    sys.exit(0 if result else 1)
