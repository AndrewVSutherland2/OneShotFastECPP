#!/usr/bin/env python3
"""Build a COMPOSITE p0 with a sequence meeting every literal condition of the
short-ECPP definition, where "point of order m*p1 on the curve over Z/p0" is
read as the order of the point in the group E(Z/p0) = E(F_l1) x E(F_l2).

Attack: order mod l1 = m (the smooth part), order mod l2 = p1 (the certified
prime).  The lcm is exactly m*p1, in the required window -- but no single prime
divisor of p0 sees the full order, so Hasse gives no contradiction.
"""
import subprocess
import sys
import random
from math import gcd, isqrt

sys.path.insert(0, "/home/claude/OneShotFastECPP/ecpp")
from vsmallECPP import ladder, verify, sieve_primes
from gen_short import build_level, factorize, x_order_is, is_probable_prime


def gp_orders(l, A_lo, A_hi):
    """[(A, N_E, N_twist)] for A in range, from one gp call."""
    s = ("{for(A=%d, %d, if (Mod(A,%d)^2 == 4, next); "
         "my(t = ellap(ellinit([0,A,0,1,0], %d))); print(A, \" \", %d+1-t, \" \", %d+1+t))}\n"
         % (A_lo, A_hi, l, l, l, l))
    r = subprocess.run(["gp", "-q"], input=s, capture_output=True, text=True)
    out = []
    for line in r.stdout.splitlines():
        f = line.split()
        if len(f) == 3:
            out.append(tuple(int(v) for v in f))
    return out


def find_x_exact(p, A, N, o, rng, tries=3000):
    div = sorted(factorize(o))
    for _ in range(tries):
        x = rng.randrange(2, p)
        f = (x * x % p * x + A * x * x + x) % p
        if f == 0:
            continue
        Xq, Zq = ladder(N // o, x, 1, A, p)
        if Zq % p == 0:
            continue
        xq = Xq * pow(Zq, -1, p) % p
        if x_order_is(p, A, xq, o, div):
            return xq
    return None


rng = random.Random(11)
M = 8                      # smooth part: a point of order 8 modulo l1

# --- pick l1, l2 ~ 2^21 primes; p0 = l1*l2 ~ 2^42 ---
def next_prime(v):
    while not is_probable_prime(v):
        v += 2
    return v


l1 = next_prime(2 ** 21 + 1001)
l2 = next_prime(2 ** 21 + 5001)
p0 = l1 * l2
n = p0.bit_length()
n2 = n * n
q0 = isqrt(p0)
L0 = q0 + 1 + isqrt(4 * q0)
print("l1 = %d, l2 = %d, p0 = %d (%d bits), n^2 = %d, L0 = %d" % (l1, l2, p0, n, n2, L0))

# --- mod l2: want order N2 = c * p1, p1 prime, M*p1 in (L0, 2*L0), p1 > n^2 ---
lo, hi = L0 // M + 1, (2 * L0 - 1) // M
print("need p1 prime in (%d, %d), and n^2 = %d" % (lo, hi, n2))
pick2 = None
for (A, NE, NT) in gp_orders(l2, 3, 4000):
    for side, N in ((1, NE), (-1, NT)):
        for c in (4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48):
            if N % c:
                continue
            cand = N // c
            if lo < cand < hi and cand > n2 and is_probable_prime(cand):
                pick2 = (A, side, N, cand)
                break
        if pick2:
            break
    if pick2:
        break
assert pick2, "no p1 found mod l2"
A2, side2, N2, p1 = pick2
o = M * p1
print("mod l2: A=%d N=%d -> p1=%d (prime), o = %d*%d = %d" % (A2, N2, p1, M, p1, o))
assert L0 < o < 2 * L0 and n2 < p1 < isqrt(p0)

# --- mod l1: want a point of exact order M = 8 ---
pick1 = None
for (A, NE, NT) in gp_orders(l1, 3, 600):
    for side, N in ((1, NE), (-1, NT)):
        if N % M == 0:
            x = find_x_exact(l1, A, N, M, rng, tries=500)
            if x is not None:
                pick1 = (A, N, x)
                break
    if pick1:
        break
assert pick1, "no order-8 point mod l1"
A1, N1, x1 = pick1
x2 = find_x_exact(l2, A2, N2, p1, rng)
assert x2 is not None
print("mod l1: A=%d N=%d x=%d (x-order %d)" % (A1, N1, x1, M))
print("mod l2: x=%d (x-order %d)" % (x2, p1))

# --- CRT curve + point ---
inv = pow(l1, -1, l2)
A = (A1 + l1 * ((A2 - A1) * inv % l2)) % p0
x = (x1 + l1 * ((x2 - x1) * inv % l2)) % p0
assert A % l1 == A1 and A % l2 == A2 and x % l1 == x1 and x % l2 == x2
print("CRT: A=%d, x=%d;  A != +-2 mod p0: %s;  gcd(A^2-4,p0)=%d"
      % (A, x, A % p0 not in (2, p0 - 2), gcd((A * A - 4) % p0, p0)))

# --- the point's order in E(Z/p0) is exactly lcm(M, p1) = o ---
print("\norder of P in E(Z/p0):")
Xo, Zo = ladder(o, x, 1, A, p0)
print("  [o]P = O ?  Z=0: %s" % (Zo % p0 == 0))
for qd in sorted(factorize(o)):
    Xq, Zq = ladder(o // qd, x, 1, A, p0)
    print("  [o/%d]P != O ?  Z!=0: %-5s  gcd(Z,p0)=%d" % (qd, Zq % p0 != 0, gcd(Zq % p0, p0)))

# --- honest terminal level for the genuinely prime p1 ---
lev = build_level(p1, n2, rng, sieve_primes(n2))
assert lev and lev[3] == 1
A_1, x_1, o_1, _ = lev
seq = [p0, A, x, o, A_1, x_1, o_1]
print("\nsequence:", " ".join(map(str, seq)))
print("p0 = %d is COMPOSITE (= %d * %d)" % (p0, l1, l2))
print("all literal conditions:")
print("  p0 >= 5, odd:                %s" % (p0 >= 5 and p0 % 2 == 1))
print("  n^2 < p1 < sqrt(p0):         %s" % (n2 < p1 < isqrt(p0)))
print("  p1 odd, p2 = 1:              %s" % (p1 % 2 == 1))
print("  m0 = %d is n^2-smooth:        %s" % (M, max(factorize(M)) <= n2))
print("  L0 < m0*p1 < r0*L0:          %s  (%d < %d < %d)" % (L0 < o < 2 * L0, L0, o, 2 * L0))
print("  0 <= A0 < p0, A0 != +-2:     %s" % (0 <= A < p0 and A % p0 not in (2, p0 - 2)))
print("  0 <= x0 < p0:                %s" % (0 <= x < p0))
print("  level-1 (mod prime p1) ok:   %s" % verify([p1, A_1, x_1, o_1]))
print("\nvsmallECPP.verify (requires the order mod EVERY prime divisor):", verify(seq))
