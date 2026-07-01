#!/usr/bin/env python3
"""
Test 3 -- classpoly's large-modulus output vs its over-Z output.

classpoly can emit a class polynomial either over Z (P=0) or reduced modulo a
given P (the "explicit CRT" / ECRT path).  ECPP uses the mod-P path with P the
large prime being certified, so it needs to be exactly right.  Test 1 anchors the
over-Z output to PARI; this test anchors the mod-P output to the over-Z output:

    for p = 2^255 - 19  (a 255-bit prime -- well beyond word size, so this
    exercises the multi-word ECRT reduction):
        for every fundamental D with |D| <= MAXD and every invariant classpoly
        supports, compute H^inv_D over Z and mod p, and check that
            (over-Z coefficients) reduced mod p  ==  (mod-p coefficients)
        coefficient by coefficient.

No PARI is needed (pure classpoly self-consistency), but because Test 1 ties the
over-Z polynomials to PARI, agreement here ties the mod-p path to PARI too.

Usage:
    python3 compare_modp.py [--maxd 1000] [--jobs 24] [--verbose]
"""

import argparse
import os
import subprocess
import sys
import concurrent.futures
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CLASSPOLY = ROOT / "classpoly_v1.0.3" / "classpoly"
PHI_DIR = ROOT / "phi_files"
WORK = ROOT / "work"
OUTDIR = WORK / "test3"

P = 2**255 - 19                      # Curve25519 field prime; 255 bits

# Every invariant classpoly supports that we exercise elsewhere (Test 1 + Test 2),
# spanning j, Weber, gamma_2, the single/double-eta quotients, Ramanujan t, and
# the Atkin invariants A_3 .. A_71.
INVARIANTS = [
    (0, "j"), (1, "f"), (2, "f^2"), (5, "gamma_2"), (6, "w_{2,3}"), (9, "w_{3,3}"),
    (10, "w_{2,5}"), (11, "t"), (12, "t^2"), (14, "w_{2,7}"), (15, "w_{3,5}"),
    (21, "w_{3,7}"), (23, "w_{2,3}^2"), (24, "w_{2,5}^2"), (26, "w_{2,13}"),
    (27, "w_{2,7}^2"), (28, "w_{3,3}^2"), (403, "w_3^12"), (405, "w_5^6"),
    (407, "w_7^4"), (413, "w_13^2"), (535, "w_{5,7}"), (539, "w_{3,13}"),
    (103, "A_3"), (105, "A_5"), (107, "A_7"), (111, "A_11"), (113, "A_13"),
    (117, "A_17"), (119, "A_19"), (123, "A_23"), (129, "A_29"), (131, "A_31"),
    (141, "A_41"), (147, "A_47"), (159, "A_59"), (171, "A_71"),
]


def is_squarefree(n):
    n = abs(n)
    i = 2
    while i * i <= n:
        if n % (i * i) == 0:
            return False
        i += 1
    return True


def fundamental_discriminants(maxd):
    """Fundamental imaginary quadratic discriminants D in [-maxd, -3]."""
    out = []
    for D in range(-3, -maxd - 1, -1):
        if D % 4 == 1 and is_squarefree(D):
            out.append(D)
        elif D % 4 == 0:
            m = D // 4
            if m % 4 in (2, 3) and is_squarefree(m):
                out.append(D)
    return out


def parse_poly(path):
    """classpoly output file -> coeff vector (highest degree first), or None.

    Header lines (I=, D=, P=) have no 'X^' and are skipped; for mod-p files the
    coefficients are already in [0, p)."""
    coeffs = {}
    for line in path.read_text().splitlines():
        if "X^" not in line:
            continue
        t = line.replace("+", " ").strip()
        if not t:
            continue
        c, e = t.split("*X^")
        coeffs[int(e.strip())] = int(c.strip())
    if not coeffs:
        return None
    deg = max(coeffs)
    return [coeffs.get(e, 0) for e in range(deg, -1, -1)]


def run_classpoly(D, inv, modulus, env):
    """modulus=0 -> over Z, else mod that integer.  Returns coeff vector or None."""
    tag = "Z" if modulus == 0 else "p"
    outfile = (OUTDIR / f"H_{-D}_{inv}_{tag}.txt").resolve()
    if outfile.exists():
        outfile.unlink()
    subprocess.run([str(CLASSPOLY), str(D), str(inv), str(modulus), str(outfile), "-1"],
                   capture_output=True, text=True, env=env, cwd=str(OUTDIR), timeout=600)
    if not outfile.exists():
        return None
    return parse_poly(outfile)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--maxd", type=int, default=1000)
    ap.add_argument("--jobs", type=int, default=24)
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    if not CLASSPOLY.exists():
        sys.exit(f"classpoly binary not found at {CLASSPOLY}; run `make` first.")
    OUTDIR.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ)
    env["CLASSPOLY_PHI_DIR"] = str(PHI_DIR)
    env["CLASSPOLY_H_DIR"] = str(OUTDIR)

    all_discs = fundamental_discriminants(args.maxd)
    # classpoly only supports the mod-P path for D < -4 (for D = -3, -4 it emits
    # the linear Hilbert polynomial over Z only); skip those two here.
    discs = [D for D in all_discs if D < -4]
    skipped = [D for D in all_discs if D >= -4]
    print(f"{len(all_discs)} fundamental discriminants with |D| <= {args.maxd} "
          f"({len(discs)} testable; {len(skipped)} skipped: D>=-4 is over-Z-only in classpoly)")
    print(f"p = 2^255 - 19 = {P}")

    tasks = [(D, code) for D in discs for (code, _l) in INVARIANTS]
    print(f"{len(tasks)} (D, invariant) pairs; {2*len(tasks)} classpoly runs "
          f"({args.jobs} parallel)")

    results = {}   # (D, inv, modulus) -> vector|None
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as ex:
        futs = {}
        for (D, code) in tasks:
            futs[ex.submit(run_classpoly, D, code, 0, env)] = (D, code, 0)
            futs[ex.submit(run_classpoly, D, code, P, env)] = (D, code, P)
        done = 0
        for fut in concurrent.futures.as_completed(futs):
            results[futs[fut]] = fut.result()
            done += 1
            if done % 2000 == 0:
                print(f"    {done}/{len(futs)} runs done")

    n_compared = n_pass = 0
    mismatches = []
    z_only = []      # over-Z produced a poly, mod-p did not
    p_only = []      # mod-p produced a poly, over-Z did not
    per_inv = {}

    for (D, code) in tasks:
        label = next(l for (c, l) in INVARIANTS if c == code)
        vZ = results.get((D, code, 0))
        vP = results.get((D, code, P))
        slot = per_inv.setdefault(code, [0, 0])
        if vZ is not None and vP is not None:
            n_compared += 1
            slot[1] += 1
            reduced = [c % P for c in vZ]
            if reduced == vP:
                n_pass += 1
                slot[0] += 1
            else:
                mismatches.append((D, label, code))
        elif vZ is not None and vP is None:
            z_only.append((D, label, code))
        elif vP is not None and vZ is None:
            p_only.append((D, label, code))

    print("\n" + "=" * 72)
    print("TEST 3 RESULTS  (classpoly mod 2^255-19  vs  over-Z reduced mod 2^255-19)")
    print("=" * 72)
    print(f"compared : {n_compared}")
    print(f"  matches : {n_pass}")
    print(f"  MISMATCHES : {len(mismatches)}")
    print(f"over-Z-only (mod-p failed to produce a poly): {len(z_only)}")
    print(f"mod-p-only (over-Z failed to produce a poly): {len(p_only)}")

    print("\nper-invariant pass/compared:")
    for (code, label) in INVARIANTS:
        if code in per_inv:
            p_, c_ = per_inv[code]
            flag = "" if p_ == c_ else "  <-- MISMATCH"
            print(f"  {label:12s} ({code:3d}): {p_}/{c_}{flag}")

    if mismatches:
        print("\n*** MISMATCHES ***")
        for (D, label, code) in mismatches[:60]:
            print(f"  D={D} {label} ({code})")
    if z_only:
        print("\nover-Z-only (mod-p path did not produce output) [sample]:")
        for (D, label, code) in z_only[:40]:
            print(f"  D={D} {label} ({code})")
    if args.verbose and p_only:
        print("\nmod-p-only [sample]:")
        for (D, label, code) in p_only[:40]:
            print(f"  D={D} {label} ({code})")

    bad = len(mismatches) + len(z_only) + len(p_only)
    print("\n" + ("OK -- every mod-p polynomial matches the reduction of the over-Z polynomial."
                  if bad == 0 else f"FAIL -- {len(mismatches)} mismatches, "
                  f"{len(z_only)} over-Z-only, {len(p_only)} mod-p-only."))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
