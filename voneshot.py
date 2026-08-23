"""
Verification of an n^4-smooth one-shot ECPP, written by Opus 4.8 Max under the
supervision of Andrew V. Sutherland.

A one-shot ECPP is a tuple (p, A, x0, m, q1, ..., qk) with:
  - p            : the integer being proved prime
  - A            : Montgomery coefficient of  E : y^2 = x^3 + A x^2 + x  over Z/pZ
  - x0           : x-coordinate of a point P on E (y-coordinate and B not needed)
  - m            : an n^4-smooth integer, n = bitlength(p), with m equal to ord(P),
                   satisfying  L < m < L*r  where L = (p^{1/4}+1)^2 (integer form
                   isqrt(p)+1+isqrt(4 isqrt p)) and r is the least prime dividing m
  - q1 < ... < qk: the prime divisors of m lying in (n^2, n^4)

It proves p prime by establishing that ord(P) = m exactly modulo every prime
divisor of p, together with the bound m > (p^{1/4}+1)^2.  If p were composite
with least prime factor l <= sqrt(p), then m | #E_l(F_l) <= (sqrt(l)+1)^2
<= (p^{1/4}+1)^2 < m, a contradiction.

The verifier trial-divides m only up to n^2; the supplied q1,...,qk account for
the prime factors of the n^2-rough part, each of which is automatically prime (a
composite divisor of the rough part below n^4 would have two factors > n^2, hence
exceed n^4).  The q_i are checked against the rough part, not against m, so that
a repeated small prime cannot be folded into a fake large prime.  The minimality
bound m < L*r keeps log m = O(n).

Target cost (FFT integer multiplication assumed):
  O(n^2 (log n)^2) bit operations and O(n^2) bits of memory.
"""

from math import gcd, isqrt


# --------------------------------------------------------------------------
# Montgomery x-only (X:Z) arithmetic on  E : y^2 = x^3 + A x^2 + x  over Z/pZ.
# Formulas depend only on A; valid on the Kummer line of E and of its twist,
# so no on-curve test of x0 is required (matching the original Pomerance set-up).
# --------------------------------------------------------------------------
def xdbl(X, Z, A, p):
    """[2](X:Z).  X' = (X^2-Z^2)^2 ;  Z' = 4 X Z (X^2 + A X Z + Z^2)."""
    XX = X * X % p
    ZZ = Z * Z % p
    XZ = X * Z % p
    X2 = (XX - ZZ) * (XX - ZZ) % p
    Z2 = 4 * XZ % p * ((XX + A * XZ + ZZ) % p) % p
    return X2, Z2


def xadd(X1, Z1, X2, Z2, Xd, Zd, p):
    """Differential addition (X1:Z1)+(X2:Z2) with known difference (Xd:Zd).
    Independent of A; the difference point may be unnormalized (Zd != 1)."""
    a = (X1 - Z1) * (X2 + Z2) % p
    b = (X1 + Z1) * (X2 - Z2) % p
    s = (a + b) % p
    d = (a - b) % p
    X3 = Zd * (s * s % p) % p
    Z3 = Xd * (d * d % p) % p
    return X3, Z3


def ladder(k, XP, ZP, A, p):
    """Montgomery ladder: returns (X:Z) of k*P for P = (XP:ZP) projective.
    Maintains (R0,R1) = (jP,(j+1)P) so the difference is always P; never
    feeds the point at infinity to xadd."""
    if k == 0:
        return (1, 0)                      # O = (1:0)
    XP %= p
    ZP %= p
    if k == 1:
        return (XP, ZP)
    Xd, Zd = XP, ZP                        # fixed difference point P
    X0, Z0 = XP, ZP                        # 1*P
    X1, Z1 = xdbl(XP, ZP, A, p)            # 2*P
    for bit in bin(k)[3:]:                 # bits below the leading 1, MSB->LSB
        if bit == '0':
            X1, Z1 = xadd(X0, Z0, X1, Z1, Xd, Zd, p)
            X0, Z0 = xdbl(X0, Z0, A, p)
        else:
            X0, Z0 = xadd(X0, Z0, X1, Z1, Xd, Zd, p)
            X1, Z1 = xdbl(X1, Z1, A, p)
    return X0, Z0


# --------------------------------------------------------------------------
# Step 1: smoothness test and prime-divisor extraction.
# --------------------------------------------------------------------------
def sieve_primes(limit):
    """All primes <= limit via a single (unsegmented) bitmap.  O(limit) bits."""
    if limit < 2:
        return []
    is_p = bytearray([1]) * (limit + 1)
    is_p[0] = is_p[1] = 0
    for i in range(2, isqrt(limit) + 1):
        if is_p[i]:
            is_p[i * i::i] = bytearray(len(is_p[i * i::i]))
    return [i for i in range(2, limit + 1) if is_p[i]]


def remainder_tree(x, mods):
    """[x % mod for mod in mods] via one product tree and a descent.
    Pads to a power of two with 1's (x % 1 == 0, sliced off at the end).
    Cost O(M(B) log k), memory O(B log k) bits, B = total size of the mods."""
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


def prime_divisors(m, primes, batch_bits):
    """Distinct primes (from `primes`) dividing m, found with batched
    remainder trees.  Primes are grouped so each batch product has ~batch_bits
    bits (comparable to m); m is reduced modulo the batch product once, then a
    remainder tree splits that down to each prime.  One batch is resident at a
    time, so memory stays O(sieve) = O(n^2) bits."""
    out = []
    batch = []
    bits = 0
    for q in primes:
        batch.append(q)
        bits += q.bit_length()
        if bits >= batch_bits:
            Q = 1
            for t in batch:
                Q *= t
            out += [q for q, r in zip(batch, remainder_tree(m % Q, batch)) if r == 0]
            batch, bits = [], 0
    if batch:
        Q = 1
        for t in batch:
            Q *= t
        out += [q for q, r in zip(batch, remainder_tree(m % Q, batch)) if r == 0]
    return out


def is_smooth(m, divisors):
    """True iff every prime factor of m lies among `divisors`.

    R = prod(divisors) is squarefree, so for any prime q | m we have
    v_q(R^N) = N if q | R and 0 otherwise.  Hence  m | R^N  iff every prime of m
    divides R and N >= max_q v_q(m).  The largest prime-power exponent of m is at
    most log2(m), so N = 2^ceil(log2 m) suffices: square R mod m  ceil(log2 m)
    times and test for 0.  No GCD.  (R = 1 gives 1 mod m, which is 0 iff m = 1,
    so the empty-divisor / m = 1 cases are handled with no special-casing.)"""
    R = 1
    for q in divisors:
        R *= q
    R %= m
    for _ in range((m - 1).bit_length()):       # ceil(log2 m) squarings
        R = R * R % m
    return R == 0


# --------------------------------------------------------------------------
# Step 2: order check.  After Q := (m/r)P, every leaf of the recursion holds
# (r/q)Q = (m/q)P, whose z-coordinate must be a unit mod p.  Each level multiplies
# by half the remaining primes, so the total scalar size per level is log r;
# with O(log t) levels the elliptic work is O(n M(n) log n) = O(n^2 (log n)^2).
# --------------------------------------------------------------------------
def check_orders(XQ, ZQ, primes, A, p):
    t = len(primes)
    if t == 0:
        return True
    if t == 1:
        return gcd(ZQ % p, p) == 1            # (m/q)P must not be O mod any l|p
    mid = t // 2
    L, Rr = primes[:mid], primes[mid:]
    hL = 1
    for q in L:
        hL *= q
    hR = 1
    for q in Rr:
        hR *= q
    XL, ZL = ladder(hL, XQ, ZQ, A, p)         # multiply in the first half, recurse on the second
    XR, ZR = ladder(hR, XQ, ZQ, A, p)         # multiply in the second half, recurse on the first
    return check_orders(XL, ZL, Rr, A, p) and check_orders(XR, ZR, L, A, p)


# --------------------------------------------------------------------------
# Top-level verifier.
# --------------------------------------------------------------------------
def verify(p, A, x0, m, qs=()):
    """Return True iff (p, A, x0, m, *qs) is a valid generalized Pomerance
    certificate (and hence p is prime).

    The certificate is (p, A, x0, m, q1, ..., qk).  qs lists the prime divisors
    of m lying in (n^2, n^4), n = bitlength(p), in strictly increasing order
    (empty when m is n^2-smooth).  Together with the primes <= n^2 found by trial
    division they account for all of m, so the verifier never factors above n^2.

    A and x0 must already be reduced: 0 <= A < p and 0 <= x0 < p (out-of-range
    inputs return False rather than being reduced mod p).  m must lie below the
    Hasse bound, m < p + 1 + 2 sqrt(p), since m = ord(P) divides #E(F_p)."""
    if p <= 3 or p % 2 == 0:       # a valid certificate requires an odd prime p > 3
        return False               # (Montgomery form needs char != 2,3; the proof needs p > 3)
    if not (0 <= A < p) or not (0 <= x0 < p):   # require canonical inputs; do not reduce silently
        return False

    # (i) size window.  L = q + 1 + floor(2 sqrt q), q = isqrt(p), is the largest
    #     order of a point on a curve over F_q for any q <= sqrt(p); m must exceed
    #     it.  m cannot exceed the Hasse bound  #E(F_p) <= p + 1 + isqrt(4p)  since
    #     m | #E.  (q + 1 + 2 sqrt q = q + 1 + isqrt(4q), increasing in q.)
    sp = isqrt(p)
    L = sp + 1 + isqrt(4 * sp)
    if m <= L:
        return False
    if m > p + 1 + isqrt(4 * p):
        return False

    # (ii) E is nonsingular modulo every prime divisor of p.
    if gcd((A * A - 4) % p, p) != 1:
        return False

    # (iii) factor m: trial-divide only up to n^2 for the small primes r_i, then
    #       splice in the supplied large primes q_i.  Trial division never goes to
    #       n^4 -- that is the whole point of carrying the q_i in the certificate.
    n = p.bit_length()
    n2 = n * n
    n4 = n2 * n2
    primes = sieve_primes(n2)
    small = prime_divisors(m, primes, batch_bits=max(64, m.bit_length()))  # ascending, each <= n^2

    #       n^2-rough part: strip every prime <= n^2 found above.  m' = rough then
    #       has every prime factor > n^2.
    rough = m
    for q in small:
        while rough % q == 0:
            rough //= q

    #       Validate the q_i against the ROUGH part m', not against m.  Each q_i
    #       must lie in (n^2, n^4), be strictly increasing, divide m', and together
    #       exhaust it (m' reduces to 1).  Because m' has every prime factor > n^2
    #       and each q_i < n^4, any q_i | m' is automatically prime -- a composite
    #       divisor would have two factors > n^2 and so exceed n^4 -- so no
    #       primality test on the q_i is required.  Validating against m' rather
    #       than m is what blocks folding a repeated small prime c <= n^2 into a
    #       fake "large prime" c*q (which divides m but not m').
    qs = list(qs)
    prev = n2
    for q in qs:
        if not (prev < q < n4):                # strictly increasing and in (n^2, n^4)
            return False
        prev = q
    rr = rough
    for q in qs:
        if rr % q != 0:                        # each listed prime must divide the rough part
            return False
        while rr % q == 0:
            rr //= q
    if rr != 1:                                # the q_i must account for all of the rough part
        return False

    divisors = small + qs                      # all distinct prime divisors of m (ascending)
    if not divisors:                           # m == 1 has no prime divisors (already excluded by m > L)
        return False

    # (iv) minimality:  m < L * r, r the least prime divisor of m.  small precedes
    #      qs and small primes are < n^2 <= q_i, so divisors[0] is the least prime.
    r = divisors[0]
    if m >= L * r:
        return False

    # (v) m*P = O, reached as the genuine point at infinity (X:0) with X a unit.
    #      If the ladder ever multiplies past ord(P) modulo some prime l | p it
    #      runs xADD(O, .) and collapses to the degenerate (0:0), which then makes
    #      Z = 0 spuriously.  Requiring gcd(X, p) = 1 rejects that collapse (and for
    #      composite p exposes a factor); this is what makes step (v) sound.  A
    #      corrupt leaf in step (vi) is already caught: (0:0) has Z = 0, failing the
    #      unit test there, so only this check needed strengthening.
    Xm, Zm = ladder(m, x0, 1, A, p)
    if Zm % p != 0 or gcd(Xm % p, p) != 1:
        return False

    # (vi) (m/q)*P != O for every prime q | m, via the divide-and-conquer tree.
    R = 1
    for q in divisors:
        R *= q
    XQ, ZQ = ladder(m // R, x0, 1, A, p)       # Q = (m/R) P
    return check_orders(XQ, ZQ, divisors, A, p)


# Valid n^4 certificates "p A x0 m q1 ... qk" (each with at least one prime in
# (n^2, n^4)), produced by the generators below and used to self-test verify().
_TEST_VECTORS = """
1000000000039 834376266027 472587544217 2240187 3539
1000000000000037 101687616200541 35220210269499 241816451 7417 32603
1000000000000000003 121042937305418467 659205001371742391 1443351938 130337
1000000000000000000000007 312714579988379121661028 899137092254762012678118 1332460718384 65141 1278439
"""


if __name__ == "__main__":
    import sys

    args = sys.argv[1:]
    if args == ["--test"]:
        # n^2-smooth certificates remain valid with no large primes (qs = ()):
        assert verify(101, 3, 24, 24) is True     # 101 prime; P=(24,*) has order 24 = 2^3*3 > 17
        assert verify(1003, 3, 24, 24) is False   # 1003 = 17 * 59
        assert verify(2, 0, 0, 1) is False        # p must be an odd prime > 3
        assert verify(3, 0, 0, 4) is False        # ditto
        assert verify(101, 104, 24, 24) is False  # A not in [0, p)  (not reduced)
        assert verify(101, 3, 125, 24) is False   # x0 not in [0, p)
        assert verify(101, 3, 24, 123) is False   # m exceeds p + 1 + floor(2 sqrt p) = 122
        # a spurious large prime out of (n^2, n^4) is rejected (n=7, n^2=49):
        assert verify(101, 3, 24, 24, (5,)) is False   # 5 <= n^2, not a valid q_i
        # fold attack: p=10000019, the point (x0=1) on A=6873344 has order 4, but the
        # cert claims m=3308=4*827.  Folding the small prime 2 into q=2*827=1654 passes
        # the literal "m mod r*q" test; validating q against the n^2-rough part (= 827)
        # rejects it, since 1654 does not divide 827.  The honest split q=827 is rejected
        # too -- (m/827)P = 4P = O exposes the true order.  Both must be False:
        assert verify(10000019, 6873344, 1, 3308, (1654,)) is False  # folded composite q
        assert verify(10000019, 6873344, 1, 3308, (827,)) is False   # true factor, wrong order
        # n^4 vectors (generated by all_smooth_pp / smooth_pp, checked end-to-end):
        for line in _TEST_VECTORS.strip().splitlines():
            vals = [int(t, 0) for t in line.split()]
            assert verify(vals[0], vals[1], vals[2], vals[3], tuple(vals[4:])) is True, line
        print("ok")
        sys.exit(0)

    if len(args) < 4:
        sys.stderr.write("usage: python voneshot.py <p> <A> <x0> <m> [q1 q2 ... qk]\n"
                         "       python voneshot.py --test\n")
        sys.exit(2)
    try:
        nums = [int(a, 0) for a in args]          # decimal or 0x-hex
    except ValueError:
        sys.stderr.write("error: p, A, x0, m, q1, ... must be integers\n")
        sys.exit(2)
    p, A, x0, m = nums[:4]
    qs = tuple(nums[4:])

    result = verify(p, A, x0, m, qs)
    print(result)
    sys.exit(0 if result else 1)
