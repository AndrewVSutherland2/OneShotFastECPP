#!/usr/bin/env python3
"""
Correctness test for the factor-base discriminant scan (dscan).

Two independent checks, for a prime of each residue class mod 4:
  (A) COMPLETENESS: dscan's set of scanned |D| (via `dumpscan`) must equal the
      set computed by an independent linear brute force -- every fundamental
      |D| < B all of whose prime-discriminant factors are QRs mod p.
  (B) SOLVABILITY: for each scanned D, dscan's solvable/not (via `dump`) must
      agree with PARI's principal-form oracle, and every reported (t,v) must
      satisfy 4p = t^2 + |D| v^2.
  (C) CAPPED FACTOR BASE: with maxfb=c the scanned set must be exactly the
      brute-force set restricted to |D| whose odd prime factors are <= c; a scan
      to B1 at cap c1 followed by Bmin=B1 maxfb=c2 fbmin=c1 must tile the scan to
      B2 at cap c2 with no repeats; and roots=<file> (cold, warm, extended,
      corrupted) must never change the output.

Usage: python3 test_dscan.py [--pbits 192] [--B 20000]
"""
import argparse, os, subprocess, sys, tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
DSCAN = HERE / "dscan"


def gp(s):
    r = subprocess.run(["gp", "-q"], input=s, capture_output=True, text=True, timeout=3600)
    if r.returncode: sys.stderr.write(r.stderr); raise RuntimeError("gp failed")
    return r.stdout


def jacobi(a, n):
    a %= n; r = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5): r = -r
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3: r = -r
        a %= n
    return r if n == 1 else 0


def factor(n):
    f = {}; d = 2
    while d * d <= n:
        while n % d == 0: f[d] = f.get(d, 0) + 1; n //= d
        d += 1 if d == 2 else 2
    if n > 1: f[n] = f.get(n, 0) + 1
    return f


def is_fundamental(d):                       # d = |D|, D = -d < 0
    if d % 4 == 3:
        return all(e == 1 for e in factor(d).values())
    if d % 4 == 0:
        m = d // 4
        return m % 4 in (1, 2) and all(e == 1 for e in factor(m).values())
    return False


def all_qr_factors(d, p, pmod8):
    """True iff every prime-discriminant factor of D=-d is a QR mod p."""
    D = -d
    if d % 2 == 1:
        two_part, odd_mag = 1, d
    else:
        two_part = None
        for delta in (-4, 8, -8):
            if d % abs(delta) == 0:
                q = D // delta
                if q % 2 == 1 and q % 4 == 1:
                    two_part = delta; break
        if two_part is None:
            return False
        odd_mag = d // abs(two_part)
    if two_part == -4 and pmod8 not in (1, 5): return False
    if two_part == 8 and pmod8 not in (1, 7): return False
    if two_part == -8 and pmod8 not in (1, 3): return False
    for q in factor(odd_mag):
        qstar_is_qr = (jacobi(p % q, q) == 1)   # (q*/p) = (p/q) = (p mod q / q)
        if not qstar_is_qr:
            return False
    return True


def check_one(p, B):
    pmod8 = p % 8
    print(f"\np = {p}  ({p.bit_length()} bits, {p % 4} mod 4),  B = {B}")

    # (A) completeness
    scanned = set(int(x) for x in subprocess.run(
        [str(DSCAN), f"p={p}", f"B={B}", "dumpscan"],
        capture_output=True, text=True, timeout=600).stdout.split())
    expected = set(d for d in range(3, B) if is_fundamental(d) and all_qr_factors(d, p, pmod8))
    missingA = expected - scanned
    extraA = scanned - expected
    print(f"  completeness: dscan {len(scanned)} vs brute {len(expected)}; "
          f"missing {len(missingA)}, extra {len(extraA)}")
    if missingA: print("    missing:", sorted(missingA)[:20])
    if extraA: print("    extra:", sorted(extraA)[:20])

    # (B) solvability + arithmetic
    out = subprocess.run([str(DSCAN), f"p={p}", f"B={B}", "dump"],
                         capture_output=True, text=True, timeout=600).stdout
    solv = {}
    arith_bad = []
    for line in out.splitlines():
        d, t, v = (int(x) for x in line.split())
        solv[d] = (t, v)
        if t * t + d * v * v != 4 * p: arith_bad.append(d)
    # PARI oracle on the scanned set
    lines = []
    for d in scanned:
        q = f"Qfb(1,1,{(d+1)//4})" if d % 4 == 3 else f"Qfb(1,0,{d//4})"
        lines.append(f'print("O {d} ", if(qfbsolve({q},{p})==[],0,1))')
    oracle = {}
    for line in gp("\n".join(lines) + "\n").splitlines():
        if line.startswith("O "):
            _, d, s = line.split(); oracle[int(d)] = int(s)
    mism = [d for d in scanned if (d in solv) != bool(oracle.get(d))]
    print(f"  solvable: dscan {len(solv)};  arith failures {len(arith_bad)};  oracle mismatches {len(mism)}")
    if arith_bad: print("    arith-bad:", arith_bad[:20])
    if mism: print("    mism:", [(d, d in solv, oracle.get(d)) for d in mism[:20]])

    ok = not missingA and not extraA and not arith_bad and not mism

    # (C) capped factor base, tiling across a growing cap, root cache
    def run(*args):
        return subprocess.run([str(DSCAN), f"p={p}"] + list(args), capture_output=True, text=True, timeout=600).stdout
    c1, c2, B1 = max(7, B // 400), max(31, B // 40), B // 5
    def oddmax(d):
        while d % 2 == 0: d //= 2
        return max(factor(d)) if d > 1 else 0
    capped = set(int(x) for x in run(f"B={B}", f"maxfb={c2}", "dumpscan").split())
    want = set(d for d in expected if oddmax(d) <= c2)
    full = sorted(run(f"B={B}", f"maxfb={c2}", "dump").splitlines())
    first = run(f"B={B1}", f"maxfb={c1}", "dump").splitlines()
    second = run(f"B={B}", f"Bmin={B1}", f"maxfb={c2}", f"fbmin={c1}", "dump").splitlines()
    tiled = sorted(first + second)
    with tempfile.TemporaryDirectory() as td:
        rf = os.path.join(td, "roots.bin")
        cold = sorted(run(f"B={B1}", f"maxfb={c1}", f"roots={rf}", "dump").splitlines())
        ext = sorted(run(f"B={B}", f"maxfb={c2}", f"roots={rf}", "dump").splitlines())
        warm = sorted(run(f"B={B}", f"maxfb={c2}", f"roots={rf}", "dump").splitlines())
        with open(rf, "r+b") as fh:                     # damage a root in the middle of the file
            fh.seek(os.path.getsize(rf) // 2); fh.write(b"\x5a\xa5\x5a")
        hurt = sorted(run(f"B={B}", f"maxfb={c2}", f"roots={rf}", "dump").splitlines())
    okC = (capped == want and tiled == full and len(set(first) & set(second)) == 0
           and cold == sorted(first) and ext == full and warm == full and hurt == full)
    print(f"  capped (maxfb={c2}): dscan {len(capped)} vs brute {len(want)}; tiling {len(first)}+{len(second)} "
          f"vs {len(full)}; root cache cold/extended/warm/damaged {'identical' if okC else 'MISMATCH'}")
    ok = ok and okC
    print("  " + ("OK" if ok else "FAIL"))
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pbits", type=int, default=192)
    ap.add_argument("--B", type=int, default=20000)
    a = ap.parse_args()
    primes = [int(gp(f"p=2^{a.pbits-1}+271828; while(p%4!={r}||!isprime(p),p++); print(p)\n").strip())
              for r in (1, 3)]
    ok = all(check_one(p, a.B) for p in primes)
    print("\n" + ("ALL OK" if ok else "FAILED"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
