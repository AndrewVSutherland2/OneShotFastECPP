#!/usr/bin/env python3
"""
Test 2 -- verify class invariants that PARI/GP's polclass does NOT support, by
mapping the roots of the invariant class polynomial back to j-invariants.

For a fundamental imaginary quadratic discriminant D (class number h) we:

  1. pick a word-size prime p that splits completely in the ring class field of
     D (equivalently, the Hilbert class polynomial H_D splits into h distinct
     linear factors mod p; equivalently 4p = t^2 + |D|v^2 is solvable);
  2. JH := the h roots of H_D mod p  (the CM j-invariants),  via PARI polclass;
  3. for each invariant `inv`: compute H^inv_D(X) mod p with classpoly, find its
     roots {x_i} mod p (these are the conjugate invariant values), and use the
     classpoly f->j modular polynomial (via the `invtoj` helper) to map each x_i
     to the j-invariant(s) J with Phi_inv(x_i, J) = 0 mod p;
  4. check that every x_i maps to at least one element of JH, and that the union
     of the mapped j-invariants lying in JH is exactly JH.

This pins down H^inv_D as a genuine class polynomial for D -- its roots are the
invariant values at the CM points and they map onto precisely the CM
j-invariants -- without needing PARI to know the invariant.  It is the only
independent check available for Atkin A_41/47/59/71, the Ramanujan t, the
single-eta quotients, and the w_{3,13} subfield/Fricke cases.

Usage:
    python3 verify_jinv.py [--maxd 500] [--jobs 16] [--pmin 100000000] [--verbose]
"""

import argparse
import os
import subprocess
import sys
import concurrent.futures
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CLASSPOLY = ROOT / "classpoly_v1.0.3" / "classpoly"
INVTOJ = ROOT / "classpoly_v1.0.3" / "invtoj"
PHI_DIR = ROOT / "phi_files"
WORK = ROOT / "work"
OUTDIR = WORK / "test2"

# Invariants to verify.  Most are NOT supported by PARI (Atkin 41..71, Ramanujan
# t/t^2, single-eta, the w_{3,13} subfield cases); a few PARI-supported ones are
# included as controls that exercise the Test-2 machinery on known-good data.
#   code, label, pari_supported(as a control?)
TEST_INVARIANTS = [
    (5,   "gamma_2",  True),    # control (jdeg 1)
    (6,   "w_{2,3}",  True),    # control (jdeg 2)
    (103, "A_3",      True),    # control (Atkin, jdeg 2)
    (131, "A_31",     True),    # control (Atkin, jdeg 2)
    (535, "w_{5,7}",  True),    # control (double-eta)
    (11,  "t",        False),   # Ramanujan (jdeg 1)
    (12,  "t^2",      False),   # Ramanujan
    (403, "w_3^12",   False),   # single-eta
    (405, "w_5^6",    False),   # single-eta
    (407, "w_7^4",    False),   # single-eta
    (413, "w_13^2",   False),   # single-eta
    (539, "w_{3,13}", False),   # double-eta (incl. subfield/Fricke cases)
    (141, "A_41",     False),   # Atkin
    (147, "A_47",     False),   # Atkin
    (159, "A_59",     False),   # Atkin -- key for the project
    (171, "A_71",     False),   # Atkin -- key for the project
]


def gp(script: str) -> str:
    # polclass needs a large PARI stack; set it on every invocation (gp reads
    # stdin line-by-line, so this must be its own leading statement).
    p = subprocess.run(["gp", "-q"], input='default(parisizemax,"8G");\n' + script,
                       capture_output=True, text=True, timeout=3600)
    if p.returncode != 0:
        sys.stderr.write(p.stderr)
        raise RuntimeError("gp failed")
    return p.stdout


def fundamental_discriminants(maxd):
    out = gp(f"for(D=-{maxd},-3, if(isfundamental(D), "
             f"print(D,\" \",qfbclassno(D))))\n")
    res = []
    for line in out.splitlines():
        if line.strip():
            d, h = line.split()
            res.append((int(d), int(h)))
    return res


def setup_disc(D, h, pmin):
    """Return (p, JH) where p splits completely and JH is the sorted Hilbert
    roots mod p; or (None, None) if no prime found in the search window."""
    # one statement per line: gp reads stdin line-by-line, so a statement must
    # not span newlines.
    script = (
        f"HD=polclass({D});h=poldegree(HD);pp=0;"
        f"forprime(q={pmin},{pmin}+10^7,if(#polrootsmod(HD,q)==h,pp=q;break));"
        f"if(pp==0,print(\"NOPRIME\"),v=vecsort(lift(polrootsmod(HD,pp)));"
        f"print1(\"P \",pp,\" JH\");for(i=1,#v,print1(\" \",v[i]));print())\n")
    out = gp(script).strip()
    if "NOPRIME" in out or not out.startswith("P "):
        return None, None
    parts = out.split()
    p = int(parts[1])
    assert parts[2] == "JH"
    JH = set(int(x) for x in parts[3:])
    return p, JH


def parse_modp_poly(path):
    """classpoly output file (coeffs already reduced mod p) -> body string for gp,
    and the polynomial degree."""
    terms = []
    deg = -1
    for line in path.read_text().splitlines():
        if "X^" not in line:
            continue
        t = line.replace("+", " ").strip()
        if t:
            terms.append(t)
            e = int(t.split("*X^")[1])
            deg = max(deg, e)
    if not terms:
        return None, -1
    return " + ".join(terms), deg


def run_classpoly_modp(D, inv, p, env):
    # classpoly self-isolates its CRT scratch (per-process dir under /tmp), so all
    # jobs can run concurrently from a shared cwd with no per-task setup.
    outfile = (OUTDIR / f"H_{-D}_{inv}_p.txt").resolve()
    if outfile.exists():
        outfile.unlink()
    subprocess.run([str(CLASSPOLY), str(D), str(inv), str(p), str(outfile), "-1"],
                   capture_output=True, text=True, env=env, cwd=str(OUTDIR),
                   timeout=600)
    if not outfile.exists():
        return None, -1
    return parse_modp_poly(outfile)


def invtoj_map(p, inv, xs, env):
    """Return dict x -> list of j's (roots of Phi_inv(x,J) mod p)."""
    proc = subprocess.run([str(INVTOJ), str(p), str(inv)],
                          input="\n".join(str(x) for x in xs) + "\n",
                          capture_output=True, text=True, env=env, timeout=120)
    out = {}
    for line in proc.stdout.splitlines():
        nums = line.split()
        if not nums:
            continue
        out[int(nums[0])] = [int(j) for j in nums[1:]]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--maxd", type=int, default=500)
    ap.add_argument("--jobs", type=int, default=16)
    ap.add_argument("--pmin", type=int, default=100_000_000)
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    for b in (CLASSPOLY, INVTOJ):
        if not b.exists():
            sys.exit(f"{b} not found; run `make` first.")
    OUTDIR.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ)
    env["CLASSPOLY_PHI_DIR"] = str(PHI_DIR)
    env["CLASSPOLY_H_DIR"] = str(OUTDIR)

    discs = fundamental_discriminants(args.maxd)
    print(f"{len(discs)} fundamental discriminants with |D| <= {args.maxd}")

    n_checked = n_pass = 0
    failures = []
    per_inv = {}   # code -> [pass, checked]
    skipped_noprime = 0

    for (D, h) in discs:
        p, JH = setup_disc(D, h, args.pmin)
        if p is None:
            skipped_noprime += 1
            continue
        if len(JH) != h:
            failures.append((D, h, "hilbert", f"|JH|={len(JH)} != h={h}"))
            continue

        # classpoly for every invariant (parallel), then map + check.
        polys = {}
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as ex:
            futs = {ex.submit(run_classpoly_modp, D, code, p, env): code
                    for (code, _l, _c) in TEST_INVARIANTS}
            for fut in concurrent.futures.as_completed(futs):
                polys[futs[fut]] = fut.result()

        # find roots of each invariant poly mod p in one gp call
        rootscript = []
        present = []
        for (code, _l, _c) in TEST_INVARIANTS:
            body, deg = polys[code]
            if body is None:
                continue
            present.append((code, deg))
            rootscript.append(
                f"f={body};print1(\"ROOTS {code}\");"
                f"v=lift(polrootsmod(f,{p}));for(i=1,#v,print1(\" \",v[i]));print()")
        roots = {}
        if rootscript:
            for line in gp("\n".join(rootscript) + "\n").splitlines():
                if line.startswith("ROOTS "):
                    parts = line.split()
                    roots[int(parts[1])] = [int(x) for x in parts[2:]]

        for (code, deg) in present:
            label = next(l for (c, l, _c) in TEST_INVARIANTS if c == code)
            xs = roots.get(code, [])
            slot = per_inv.setdefault(code, [0, 0])
            slot[1] += 1
            n_checked += 1

            # invariant poly must split completely mod p (p splits completely)
            if len(xs) != deg:
                failures.append((D, h, label, f"deg={deg} but {len(xs)} roots mod p"))
                continue
            mapping = invtoj_map(p, code, xs, env)
            covered = set()
            bad_x = None
            for x in xs:
                js = [j for j in mapping.get(x, []) if j in JH]
                if not js:
                    bad_x = x
                    break
                covered.update(js)
            if bad_x is not None:
                failures.append((D, h, label, f"invariant root {bad_x} maps to no CM j"))
                continue
            if covered != JH:
                failures.append((D, h, label,
                                 f"covered {len(covered)}/{len(JH)} CM j-invariants"))
                continue
            n_pass += 1
            slot[0] += 1
            if args.verbose:
                print(f"  D={D} h={h} {label}: deg={deg} -> covers all {h} j-invariants OK")

    print("\n" + "=" * 72)
    print("TEST 2 RESULTS  (invariant class poly roots -> j-invariants)")
    print("=" * 72)
    print(f"checks run : {n_checked}")
    print(f"  passed   : {n_pass}")
    print(f"  FAILED   : {len(failures)}")
    print(f"discriminants skipped (no split prime in window): {skipped_noprime}")
    print("\nper-invariant pass/checked:")
    for (code, label, ctrl) in TEST_INVARIANTS:
        if code in per_inv:
            p_, c_ = per_inv[code]
            tag = " [control]" if ctrl else ""
            flag = "" if p_ == c_ else "  <-- FAILURE"
            print(f"  {label:12s} ({code:3d}): {p_}/{c_}{tag}{flag}")
    if failures:
        print("\n*** FAILURES ***")
        for (D, h, label, msg) in failures[:60]:
            print(f"  D={D} h={h} {label}: {msg}")

    print("\n" + ("OK -- every invariant class polynomial maps onto the CM j-invariants."
                  if not failures else f"FAIL -- {len(failures)} checks failed."))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
