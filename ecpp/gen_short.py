#!/usr/bin/env python3
"""gen_short.py -- generate short ECPPs (ShortPrimalityProofs format).

Modes:
  gen_short.py p256 [out=certs/short/short_ecpp_p256.txt]
      Exhaustively enumerate ALL short ECPPs (p, A, x, m) for primes p <= 256.
      For p <= 256 every short ECPP has a single level: a second level would
      need an n^2-rough modulus p_1 with n^2 < p_1 < sqrt(p_0) <= 16 while
      n^2 >= 4^2 = 16 -- impossible.  Every emitted tuple is cross-checked
      three ways: pruned enumeration (twist orders + ladder exactness), a
      naive affine group-order computation, and vsmallECPP.verify().

  gen_short.py chain pbits=48 seed=1
      Build a multi-level short ECPP for a random prime of the given size
      (PARI gp's ellap supplies the point counts), verify it with
      vsmallECPP.verify(), and print the flat sequence.

  gen_short.py tampers
      Emit tamper test vectors (each must be rejected by the verifier).
"""

import os
import subprocess
import sys
import random
from math import gcd, isqrt

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from vsmallECPP import ladder, verify, sieve_primes


# ---------------------------------------------------------------- utilities
def factorize(m):
    """Full trial-division factorization; fine for the small m used here."""
    fs = {}
    d = 2
    while d * d <= m:
        while m % d == 0:
            fs[d] = fs.get(d, 0) + 1
            m //= d
        d += 1
    if m > 1:
        fs[m] = fs.get(m, 0) + 1
    return fs


def x_order_is(p, A, x, o, prime_divs):
    """Exact x-order test via the Montgomery ladder (p prime)."""
    Xo, Zo = ladder(o, x, 1, A, p)
    if Zo % p != 0 or gcd(Xo % p, p) != 1:
        return False
    for q in prime_divs:
        Xq, Zq = ladder(o // q, x, 1, A, p)
        if Zq % p == 0:
            return False
    return True


def naive_x_order(p, A, x):
    """Order of (x, y) on B y^2 = x^3 + A x^2 + x, computed by repeated affine
    addition on the model B = f(x), y = 1 (any realization gives the same
    order).  Independent of the ladder code; p must be prime."""
    f = (x * x % p * x + A * x * x + x) % p
    if f == 0:
        return 2 if x != 0 or True else 2      # (x, 0) is 2-torsion
    B = f
    P = (x, 1)

    def add(P1, P2):
        x1, y1 = P1
        x2, y2 = P2
        if (x1 - x2) % p == 0:
            if (y1 + y2) % p == 0:
                return None
            lam = (3 * x1 * x1 + 2 * A * x1 + 1) * pow(2 * B * y1 % p, -1, p) % p
        else:
            lam = (y2 - y1) * pow((x2 - x1) % p, -1, p) % p
        x3 = (B * lam * lam - A - x1 - x2) % p
        return (x3, (lam * (x1 - x3) - y1) % p)

    Q = P
    k = 1
    limit = p + 1 + isqrt(4 * p) + 2
    while k <= limit:
        Q = add(Q, P)
        k += 1
        if Q is None:
            return k
    return -1


# ---------------------------------------------------------------- p256 mode
def window_ms(p, n2):
    """All valid m for a terminal level mod p: n^2-smooth, L < m < r*L (and
    m <= Hasse, implied by m being a point order)."""
    q = isqrt(p)
    L = q + 1 + isqrt(4 * q)
    H = p + 1 + isqrt(4 * p)
    out = []
    for m in range(L + 1, H + 1):
        fs = factorize(m)
        if max(fs) > n2:
            continue
        r = min(fs)
        if m < r * L:
            out.append((m, sorted(fs)))
    return out


def all_short_for_p(p, n2):
    """All (A, x, m) completing a terminal short-ECPP level mod prime p."""
    out = []
    Ms = window_ms(p, n2)
    if not Ms:
        return out
    sq = set(y * y % p for y in range(1, p))
    for A in range(p):
        if (A * A - 4) % p == 0:
            continue
        s = 0
        chi = [0] * p
        for x in range(p):
            f = (x * x % p * x + A * x * x + x) % p
            if f != 0:
                chi[x] = 1 if f in sq else -1
                s += chi[x]
        NE, NT = p + 1 + s, p + 1 - s
        for m, div in Ms:
            for side, N in ((1, NE), (-1, NT)):
                if N % m:
                    continue
                for x in range(p):
                    if chi[x] == side and x_order_is(p, A, x, m, div):
                        out.append((p, A, x, m))
    return sorted(set(out))


def mode_p256(out_path):
    outdir = os.path.dirname(out_path)
    if outdir:
        os.makedirs(outdir, exist_ok=True)
    primes = [p for p in range(3, 257) if all(p % d for d in range(2, isqrt(p) + 1))]
    lines = []
    summary = []
    for p in primes:
        n = p.bit_length()                      # ceil(log2 p) for odd p
        certs = all_short_for_p(p, n * n)
        for (pp, A, x, m) in certs:
            assert naive_x_order(pp, A, x) == m, (pp, A, x, m)   # oracle 2
            assert verify([pp, A, x, m]), (pp, A, x, m)          # oracle 3
        lines += ["%d %d %d %d" % c for c in certs]
        ms = sorted(set(m for (_, _, _, m) in certs))
        summary.append((p, len(certs), ms))
        print("p=%3d: %5d certs, m in %s" % (p, len(certs), ms))
    with open(out_path, "w") as f:
        f.write("# all short ECPPs (p A x0 m) for primes p <= 256\n")
        f.write("\n".join(lines) + "\n")
    none = [p for (p, c, _) in summary if c == 0]
    print("\ntotal: %d certificates for %d primes -> %s" %
          (len(lines), len(primes), out_path))
    print("primes with NO short ECPP: %s (and 2 is excluded as even)" % none)


# ---------------------------------------------------------------- chain mode
def gp_ellap(p, A):
    r = subprocess.run(["gp", "-q"], input="print(ellap(ellinit([0,%d,0,1,0],%d)))\n" % (A, p),
                       capture_output=True, text=True)
    return int(r.stdout.strip().splitlines()[-1])


def is_probable_prime(nv):
    if nv < 2:
        return False
    for b in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if nv % b == 0:
            return nv == b
    d, s = nv - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for b in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        xv = pow(b, d, nv)
        if xv in (1, nv - 1):
            continue
        for _ in range(s - 1):
            xv = xv * xv % nv
            if xv == nv - 1:
                break
        else:
            return False
    return True


def divisors_of(fs):
    ds = [1]
    for q, e in fs.items():
        ds = [d * q ** i for d in ds for i in range(e + 1)]
    return sorted(ds)


def level_from_order(p, A, side, N, o, div, rng):
    """Find x with x-order exactly o on the `side` twist of E_A mod p."""
    for _ in range(400):
        x = rng.randrange(2, p)
        f = (x * x % p * x + A * x * x + x) % p
        if f == 0:
            continue
        chi = 1 if pow(f, (p - 1) // 2, p) == 1 else -1
        if chi != side:
            continue
        Xq, Zq = ladder(N // o, x, 1, A, p)
        if Zq % p == 0:
            continue
        xq = Xq * pow(Zq, -1, p) % p
        if x_order_is(p, A, xq, o, div):
            return xq
    return None


def build_level(p, n2, rng, smallp):
    """Search A for a level mod p: order N with an n^2-rough prime cofactor
    p_next (p_next^2 < p) or fully smooth (terminal), o in the window."""
    q = isqrt(p)
    L = q + 1 + isqrt(4 * q)
    for _ in range(4000):
        A = rng.randrange(3, p - 2)
        if (A * A - 4) % p == 0:
            continue
        t = gp_ellap(p, A)
        for side, N in ((1, p + 1 - t), (-1, p + 1 + t)):
            sm, rest = 1, N
            fs = {}
            for qq in smallp:
                while rest % qq == 0:
                    rest //= qq
                    sm *= qq
                    fs[qq] = fs.get(qq, 0) + 1
            p_next = rest
            if p_next > 1:
                if p_next * p_next >= p or not is_probable_prime(p_next):
                    continue
            # choose o = m * p_next with m | sm, L < o < r(m) * L
            for m in divisors_of(fs):
                if m < 2:
                    continue
                o = m * p_next
                r = min(factorize(m))
                if L < o < r * L:
                    div = sorted(factorize(m)) + ([p_next] if p_next > 1 else [])
                    x = level_from_order(p, A, side, N, o, div, rng)
                    if x is not None:
                        return (A, x, o, p_next)
    return None


def mode_chain(pbits, seed):
    rng = random.Random(seed)
    while True:
        p0 = rng.getrandbits(pbits) | (1 << (pbits - 1)) | 1
        if is_probable_prime(p0):
            break
    n = p0.bit_length()
    n2 = n * n
    smallp = sieve_primes(n2)
    print("p0 = %d (%d bits), n^2 = %d" % (p0, n, n2))
    seq = [p0]
    p = p0
    while p > 1:
        lev = build_level(p, n2, rng, smallp)
        if lev is None:
            sys.exit("level search failed at p=%d (retry with another seed)" % p)
        A, x, o, p_next = lev
        print("  level: p=%d  A=%d  x=%d  o=%d  -> p_next=%d" % (p, A, x, o, p_next))
        seq += [A, x, o]
        p = p_next
    ok = verify(seq)
    print("chain of %d levels, vsmallECPP.verify -> %s" % ((len(seq) - 1) // 3, ok))
    print(" ".join(map(str, seq)))
    return seq


# ---------------------------------------------------------------- tampers
def mode_tampers():
    base = None
    for p in [251]:
        cs = all_short_for_p(p, p.bit_length() ** 2)
        if cs:
            base = cs[0]
            break
    p, A, x, m = base
    bad = [
        [p, A, x, m * 2],          # breaks window / exactness
        [p, A, (x + 1) % p, m],    # wrong point
        [p, A, x, m, A, x, m],     # trailing level with p_next = 1 modulus
        [221, 5, 2, 34],           # composite p0 with shape-plausible values
        [p, 2, x, m],              # singular curve A = 2
    ]
    for b in bad:
        print(" ".join(map(str, b)), "->", verify(b))


if __name__ == "__main__":
    args = dict(a.split("=", 1) for a in sys.argv[2:] if "=" in a)
    mode = sys.argv[1] if len(sys.argv) > 1 else ""
    if mode == "p256":
        mode_p256(args.get("out", "certs/short/short_ecpp_p256.txt"))
    elif mode == "chain":
        mode_chain(int(args.get("pbits", 48)), int(args.get("seed", 1)))
    elif mode == "tampers":
        mode_tampers()
    else:
        print(__doc__)
