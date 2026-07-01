#!/usr/bin/env python3
"""
Correctness tests for the batched n^4-smoothness engine (component b).

  Test A (engine):   `smoothtest parts` smooth-part/cofactor for a batch of mixed
                     random integers must match PARI's factor-based smooth part,
                     over several smoothness bounds and batch sizes.

  Test B (gate):     end-to-end.  For a real prime p, feed the solvable (d,t,v)
                     from `dscan ... dump` into `smoothtest gate`; every reported
                     winner (p, order N=p+1-/+t, m, q_i) is checked *independently
                     with PARI*: m | N, m is n^4-smooth, L < m < min(Hasse, L*r),
                     and q_i are exactly the primes of m in (n^2, n^4).  Also the
                     gate decision itself (smooth part > L) is cross-checked
                     against PARI for the whole candidate set.

Usage: python3 test_smooth.py [--pbits 128] [--B 300000]
"""
import argparse, subprocess, sys, random
from pathlib import Path

HERE = Path(__file__).resolve().parent
SMOOTH = HERE / "smoothtest"
DSCAN = HERE / "dscan"


def gp(s):
    r = subprocess.run(["gp", "-q"], input=s, capture_output=True, text=True, timeout=3600)
    if r.returncode:
        sys.stderr.write(r.stderr); raise RuntimeError("gp failed")
    return r.stdout


def pari_smooth_part(nums, y):
    """List of (smooth_part, cofactor) for each n in nums, via PARI factor."""
    lines = [f"v={x};f=factor(v);s=prod(i=1,#f~,if(f[i,1]<={y},f[i,1]^f[i,2],1));print(s,\" \",v/s)"
             for x in nums]
    out = gp("\n".join(lines) + "\n").split()
    return [(int(out[2*i]), int(out[2*i+1])) for i in range(len(nums))]


def test_A():
    print("=== Test A: batched smooth-part extraction vs PARI ===")
    rng = random.Random(12345)
    ok = True
    for y, count in [(1000, 50), (10**6, 400), (2**28, 200)]:
        nums = []
        for _ in range(count):
            kind = rng.random()
            if kind < 0.3:                              # smooth-ish: product of small primes
                x = 1
                for _ in range(rng.randint(1, 25)):
                    x *= rng.choice([2, 3, 5, 7, 11, 13, 101, 9973, 99991])
                nums.append(x)
            elif kind < 0.6:                            # smooth * a big prime
                x = 1
                for _ in range(rng.randint(1, 12)):
                    x *= rng.choice([2, 3, 5, 7, 9973, 99991])
                big = int(gp(f"print(nextprime({rng.randint(10**8,10**12)}))"))
                nums.append(x * big)
            else:                                       # arbitrary
                nums.append(rng.randint(2, 10**rng.randint(3, 30)))
        ours = subprocess.run([str(SMOOTH), "parts", f"y={y}"],
                              input="\n".join(map(str, nums)) + "\n",
                              capture_output=True, text=True).stdout.split("\n")
        ours = [l.split() for l in ours if l.strip()]
        pari = pari_smooth_part(nums, y)
        bad = 0
        for i, x in enumerate(nums):
            os, oc = int(ours[i][1]), int(ours[i][2])
            ps, pc = pari[i]
            if os != ps or oc != pc or int(ours[i][0]) != x:
                bad += 1
                if bad <= 5:
                    print(f"  MISMATCH y={y} N={x}: ours=({os},{oc}) pari=({ps},{pc})")
        print(f"  y={y:>10} batch={count}: {'OK' if bad==0 else f'{bad} MISMATCHES'}")
        ok = ok and bad == 0
    return ok


def factor_pari(n):
    out = gp(f"f=factor({n});for(i=1,#f~,print(f[i,1],\" \",f[i,2]))")
    d = {}
    for line in out.splitlines():
        a, b = line.split(); d[int(a)] = int(b)
    return d


def test_B(pbits, B):
    print(f"\n=== Test B: end-to-end gate (pbits={pbits}, B={B}) ===")
    # dscan and smoothtest share gen_p(pbits,seed); drive both with the same pbits/seed
    # and recover the actual p from smoothtest's stderr.
    seed = 1
    dump = subprocess.run([str(DSCAN), f"pbits={pbits}", f"seed={seed}", f"B={B}", "dump"],
                          capture_output=True, text=True, timeout=1200).stdout
    ncand = 2 * len(dump.splitlines())
    print(f"  solvable D: {len(dump.splitlines())}  -> candidates: {ncand}")
    res = subprocess.run([str(SMOOTH), "gate", f"pbits={pbits}", f"seed={seed}"],
                         input=dump, capture_output=True, text=True, timeout=1200)
    # parse p, n, L from stderr
    p = L = Hass = None; n = n2 = n4 = None
    for line in res.stderr.splitlines():
        if line.startswith("p = "): p = int(line[4:])
        if line.startswith("n=") and "n^2" in line:
            parts = line.replace("=", " ").split()
            n = int(parts[1]); n2 = int(parts[3]); n4 = int(parts[5])
        if line.startswith("L = "): L = int(line.split()[2])
    Hass = p + 1 + int(gp(f"print(sqrtint(4*{p}))"))
    print(f"  p bit_length={p.bit_length()}  n={n} n^2={n2} n^4={n4}")
    print("  " + "\n  ".join(l for l in res.stderr.splitlines() if "winners" in l or "us/candidate" in l))

    wins = [l for l in res.stdout.splitlines() if l.startswith("WIN")]
    print(f"  winners reported: {len(wins)}")
    if any("SELFCHECK-FAIL" in w for w in wins):
        print("  C SELF-CHECK FAILED"); return False

    ok = True
    # independently verify each winner with PARI
    for w in wins:
        # WIN D=-d t=T order=N  m=M  smoothpart=S  q=[...]
        toks = w.replace("=", " ").replace("[", " ").replace("]", " ").split()
        d = int(toks[toks.index("D") + 1].lstrip("-")) if "D" in toks else None
        N = int(toks[toks.index("order") + 1])
        M = int(toks[toks.index("m") + 1])
        S = int(toks[toks.index("smoothpart") + 1])
        qidx = w.find("q=[")
        qs = w[qidx + 3:w.find("]", qidx)]
        qs = [int(x) for x in qs.split(",") if x.strip()]
        errs = []
        if N % M != 0: errs.append("m does not divide N")
        if not (L < M <= Hass): errs.append(f"m out of (L,Hasse]")
        fac = factor_pari(M)
        r = min(fac)                                    # least prime of m
        if not (M < L * r): errs.append("m >= L*r")
        if max(fac) > n4: errs.append(f"m not n^4-smooth (max {max(fac)})")
        expected_q = sorted([q for q in fac if n2 < q < n4])
        if expected_q != sorted(qs): errs.append(f"q_i mismatch: got {sorted(qs)} want {expected_q}")
        # S must be the true n^4-smooth part of N and exceed L
        ps, _ = pari_smooth_part([N], n4)[0]
        if ps != S: errs.append(f"reported smoothpart {S} != PARI {ps}")
        if not (S > L): errs.append("gate: S<=L but reported winner")
        if errs:
            ok = False
            print(f"  WINNER BAD (D=-{d}): {errs}")
    if wins and ok:
        print(f"  all {len(wins)} winners independently verified vs PARI: OK")

    # cross-check the GATE decision on the full candidate set vs PARI:
    # recompute candidate orders, compare (smooth part > L) our-vs-PARI.
    cand = []
    for line in dump.splitlines():
        dd, t, v = line.split()
        t = int(t)
        cand += [p + 1 - t, p + 1 + t]
    # our smooth parts via parts mode with y=n^4
    ours = subprocess.run([str(SMOOTH), "parts", f"y={n4}"],
                          input="\n".join(map(str, cand)) + "\n",
                          capture_output=True, text=True).stdout.split("\n")
    ours = {int(l.split()[0]): int(l.split()[1]) for l in ours if l.strip()}
    # PARI in one batch
    pari = pari_smooth_part(cand, n4)
    mism = 0
    for i, N in enumerate(cand):
        if ours.get(N) != pari[i][0]:
            mism += 1
            if mism <= 5: print(f"    smoothpart mismatch N={N}: ours={ours.get(N)} pari={pari[i][0]}")
    gate_ours = sum(1 for N in cand if ours.get(N, 0) > L)
    gate_pari = sum(1 for i in range(len(cand)) if pari[i][0] > L)
    print(f"  full-set smoothpart mismatches: {mism};  gate count ours={gate_ours} pari={gate_pari}")
    ok = ok and mism == 0 and gate_ours == gate_pari
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pbits", type=int, default=128)
    ap.add_argument("--B", type=int, default=300000)
    a = ap.parse_args()
    okA = test_A()
    okB = test_B(a.pbits, a.B)
    print("\n" + ("ALL OK" if (okA and okB) else "FAILED"))
    return 0 if (okA and okB) else 1


if __name__ == "__main__":
    sys.exit(main())
