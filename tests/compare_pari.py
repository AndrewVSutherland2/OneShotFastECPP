#!/usr/bin/env python3
"""
Test 1 -- cross-check classpoly against PARI/GP's polclass.

For every fundamental imaginary quadratic discriminant D with |D| <= MAXD, and
every class invariant supported by BOTH classpoly and PARI/GP, we compute the
class polynomial H^inv_D(X) over Z with each tool and compare the coefficient
vectors exactly.

classpoly invariant codes and PARI invariant codes coincide for the shared
invariants (both derive from the Enge-Sutherland framework), except for the two
double-eta quotients w_{5,7} and w_{3,13}, where PARI uses 35/39 and classpoly
uses 535/539 (the 500+p1*p2 encoding).  The INVARIANTS table below records the
(pari_code, classpoly_code) pairs.

A pair (D, inv) is only *compared* when both tools produce a polynomial.  When
exactly one tool accepts (D, inv) we record a "validity difference" -- this is
expected for a handful of edge cases where the two tools impose slightly
different domain conditions (e.g. small-D exceptions, class-number bounds), and
is reported separately, not as a failure.

Usage:
    python3 compare_pari.py [--maxd 1000] [--jobs 24] [--verbose]
"""

import argparse
import os
import subprocess
import sys
import concurrent.futures
from pathlib import Path

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
ROOT = Path(__file__).resolve().parent.parent
CLASSPOLY = ROOT / "classpoly_v1.0.3" / "classpoly"
PHI_DIR = ROOT / "phi_files"
WORK = ROOT / "work"
OUTDIR = WORK / "test1"          # classpoly output files land here

# (pari_code, classpoly_code, human label).  Equal codes for everything except
# the two double-eta quotients that PARI numbers 35/39 and classpoly 535/539.
INVARIANTS = [
    (0,   0,   "j"),
    (1,   1,   "f (Weber)"),
    (2,   2,   "f^2"),
    (5,   5,   "gamma_2"),
    (6,   6,   "w_{2,3}"),
    (9,   9,   "w_{3,3}"),
    (10,  10,  "w_{2,5}"),
    (14,  14,  "w_{2,7}"),
    (15,  15,  "w_{3,5}"),
    (21,  21,  "w_{3,7}"),
    (23,  23,  "w_{2,3}^2"),
    (24,  24,  "w_{2,5}^2"),
    (26,  26,  "w_{2,13}"),
    (27,  27,  "w_{2,7}^2"),
    (28,  28,  "w_{3,3}^2"),
    (35,  535, "w_{5,7}"),
    (39,  539, "w_{3,13}"),
    (103, 103, "A_3"),
    (105, 105, "A_5"),
    (107, 107, "A_7"),
    (111, 111, "A_11"),
    (113, 113, "A_13"),
    (117, 117, "A_17"),
    (119, 119, "A_19"),
    (123, 123, "A_23"),
    (129, 129, "A_29"),
    (131, 131, "A_31"),
]


def gp(script: str) -> str:
    """Run a gp script (string) and return stdout."""
    p = subprocess.run(["gp", "-q"], input=script, capture_output=True,
                       text=True, timeout=3600)
    if p.returncode != 0:
        sys.stderr.write(p.stderr)
        raise RuntimeError("gp failed")
    return p.stdout


def fundamental_discriminants(maxd: int):
    """Return sorted list of fundamental D in [-maxd, -3] with class number."""
    out = gp(f"for(D=-{maxd},-3, if(isfundamental(D), "
             f"print(D,\" \",qfbclassno(D))))\n")
    res = []
    for line in out.splitlines():
        line = line.strip()
        if not line:
            continue
        d, h = line.split()
        res.append((int(d), int(h)))
    return res


# ---------------------------------------------------------------------------
# classpoly side
# ---------------------------------------------------------------------------
def run_classpoly(D: int, cp_inv: int, env: dict):
    """Run classpoly for (D, cp_inv); return coeff vector (high->low) or None.

    classpoly now isolates its own CRT scratch in a private per-process directory
    (mkdtemp under /tmp), so we can launch all jobs concurrently from a shared
    working directory with no per-task setup.
    """
    outfile = (OUTDIR / f"H_{-D}_{cp_inv}.txt").resolve()
    if outfile.exists():
        outfile.unlink()
    # classpoly D inv P filename verbosity  (P=0 -> over Z, verbosity -1 -> quiet)
    subprocess.run(
        [str(CLASSPOLY), str(D), str(cp_inv), "0", str(outfile), "-1"],
        capture_output=True, text=True, env=env, cwd=str(OUTDIR), timeout=600)
    if not outfile.exists():
        return None  # classpoly rejected (D, inv) or failed
    return parse_classpoly_file(outfile)


def parse_classpoly_file(path: Path):
    """Parse a classpoly output file into a coeff vector, highest degree first."""
    coeffs = {}
    for line in path.read_text().splitlines():
        if "X^" not in line:
            continue                      # skip I=.. and D=.. header lines
        term = line.replace("+", " ").strip()
        if not term:
            continue
        c_str, e_str = term.split("*X^")
        coeffs[int(e_str.strip())] = int(c_str.strip())
    if not coeffs:
        return None
    deg = max(coeffs)
    return [coeffs.get(e, 0) for e in range(deg, -1, -1)]


# ---------------------------------------------------------------------------
# PARI side -- compute all polclass vectors in a single gp process
# ---------------------------------------------------------------------------
def run_pari_batch(pairs):
    """pairs: list of (D, pari_inv).  Return dict (D,pari_inv) -> vec or None."""
    lines = ['default(parisizemax,"8G");']
    for (D, inv) in pairs:
        # print "R D inv OK c,c,..." or "R D inv FAIL"; print1 avoids line-wrapping
        lines.append(
            f"iferr(my(v=Vec(polclass({D},{inv})));"
            f"print1(\"R {D} {inv} OK \");"
            f"for(i=1,#v,print1(v[i],if(i<#v,\",\",\"\")));print(),"
            f"E,print(\"R {D} {inv} FAIL\"));")
    out = gp("\n".join(lines) + "\n")
    res = {}
    for line in out.splitlines():
        if not line.startswith("R "):
            continue
        parts = line.split(" ", 4)
        _, d, inv, status = parts[0], int(parts[1]), int(parts[2]), parts[3]
        if status == "FAIL":
            res[(d, inv)] = None
        else:
            vec = [int(x) for x in parts[4].split(",")]
            res[(d, inv)] = vec
    return res


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
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

    print(f"Enumerating fundamental discriminants with |D| <= {args.maxd} ...")
    discs = fundamental_discriminants(args.maxd)
    print(f"  {len(discs)} fundamental discriminants")

    # Build the full list of candidate (D, pari_inv, cp_inv) triples.
    triples = [(D, h, pi, ci, label)
               for (D, h) in discs
               for (pi, ci, label) in INVARIANTS]
    print(f"  {len(triples)} candidate (D, invariant) pairs to check")

    # --- classpoly: run all pairs in parallel ---
    print(f"Running classpoly ({args.jobs} parallel jobs) ...")
    cp_res = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as ex:
        futs = {ex.submit(run_classpoly, D, ci, env): (D, ci)
                for (D, _h, _pi, ci, _l) in triples}
        done = 0
        for fut in concurrent.futures.as_completed(futs):
            cp_res[futs[fut]] = fut.result()
            done += 1
            if done % 500 == 0:
                print(f"    {done}/{len(futs)} classpoly runs done")

    # --- PARI: one batched gp process ---
    print("Running PARI polclass (single gp process) ...")
    pari_pairs = sorted({(D, pi) for (D, _h, pi, _ci, _l) in triples})
    pari_res = run_pari_batch(pari_pairs)

    # --- compare ---
    n_compared = n_pass = 0
    mismatches = []
    only_pari = []      # PARI accepts, classpoly rejects
    only_cp = []        # classpoly accepts, PARI rejects
    per_inv = {}        # label -> [pass, compared]

    for (D, h, pi, ci, label) in triples:
        cv = cp_res.get((D, ci))
        pv = pari_res.get((D, pi))
        slot = per_inv.setdefault(label, [0, 0])
        if cv is not None and pv is not None:
            n_compared += 1
            slot[1] += 1
            if cv == pv:
                n_pass += 1
                slot[0] += 1
            else:
                mismatches.append((D, h, label, pi, ci, cv, pv))
        elif pv is not None and cv is None:
            only_pari.append((D, h, label, pi, ci))
        elif cv is not None and pv is None:
            only_cp.append((D, h, label, pi, ci))

    # --- report ---
    print("\n" + "=" * 72)
    print("RESULTS")
    print("=" * 72)
    print(f"compared (both tools produced a polynomial): {n_compared}")
    print(f"  exact matches : {n_pass}")
    print(f"  MISMATCHES    : {len(mismatches)}")
    print(f"validity differences (only one tool accepted): "
          f"{len(only_pari)} PARI-only, {len(only_cp)} classpoly-only")

    print("\nper-invariant pass/compared:")
    for (pi, ci, label) in INVARIANTS:
        if label in per_inv:
            p, c = per_inv[label]
            flag = "" if p == c else "  <-- MISMATCH"
            print(f"  {label:12s} (pari {pi:3d} / cp {ci:3d}): {p}/{c}{flag}")

    if mismatches:
        print("\n*** MISMATCHES ***")
        for (D, h, label, pi, ci, cv, pv) in mismatches[:50]:
            print(f"  D={D} h={h} {label} (pari {pi}/cp {ci})")
            print(f"      classpoly: {cv}")
            print(f"      pari:      {pv}")

    if args.verbose and only_pari:
        print("\nPARI-only (classpoly rejected) [sample]:")
        for (D, h, label, pi, ci) in only_pari[:40]:
            print(f"  D={D} h={h} {label} (pari {pi}/cp {ci})")
    if args.verbose and only_cp:
        print("\nclasspoly-only (PARI rejected) [sample]:")
        for (D, h, label, pi, ci) in only_cp[:40]:
            print(f"  D={D} h={h} {label} (pari {pi}/cp {ci})")

    # write a machine-readable summary
    summary = OUTDIR.parent / "test1_summary.txt"
    with open(summary, "w") as f:
        f.write(f"maxd {args.maxd}\n")
        f.write(f"discriminants {len(discs)}\n")
        f.write(f"compared {n_compared}\n")
        f.write(f"pass {n_pass}\n")
        f.write(f"mismatch {len(mismatches)}\n")
        f.write(f"only_pari {len(only_pari)}\n")
        f.write(f"only_cp {len(only_cp)}\n")
        for (D, h, label, pi, ci, cv, pv) in mismatches:
            f.write(f"MISMATCH D={D} h={h} {label} pari={pi} cp={ci}\n")
        for (D, h, label, pi, ci) in only_pari:
            f.write(f"ONLY_PARI D={D} h={h} {label} pari={pi} cp={ci}\n")
        for (D, h, label, pi, ci) in only_cp:
            f.write(f"ONLY_CP D={D} h={h} {label} pari={pi} cp={ci}\n")
    print(f"\nsummary written to {summary}")

    print("\n" + ("OK -- all compared polynomials agree."
                  if not mismatches else
                  f"FAIL -- {len(mismatches)} mismatches."))
    return 1 if mismatches else 0


if __name__ == "__main__":
    sys.exit(main())
