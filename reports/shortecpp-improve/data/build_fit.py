#!/usr/bin/env python3
"""Build time of a winner against what is known of its discriminant before the build.

From every log of the c8a runs: the cm_method lines (degree of the class polynomial, its
decomposition m x n, CRT and root-finding seconds) with the size of the level they belong to.
The class number h(D) comes from PARI; the question is how well h alone predicts the cost:
the root-finding step at degree n = the larger part of the best split of h."""
import glob, json, math, os, re, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
CM = re.compile(r"cm_method D=-(\d+): H_D\^inv\(deg (\d+)(?: = (\d+) x (\d+))?, inv (\d+)\): ([\d.]+)s\s+roots?: ([\d.]+)s(?: \+ ([\d.]+)s)?")
BUILT = re.compile(r"level (\d+) bits: built D=-(\d+) h~\d+ in ([\d.]+)s")
rows = []
for d in ("rule_c8a", "yafu_c8a_1", "yafu_c8a_2", "yafu_c8a_3"):
    for fn in sorted(glob.glob(os.path.join(HERE, d, "*.log"))):
        pend = {}
        for ln in open(fn):
            m = CM.search(ln)
            if m:
                pend[int(m.group(1))] = m
                continue
            b = BUILT.search(ln)
            if b and int(b.group(2)) in pend:
                m = pend[int(b.group(2))]
                rows.append({"set": d, "file": os.path.basename(fn), "bits": int(b.group(1)), "D": int(m.group(1)), "deg": int(m.group(2)),
                             "m": int(m.group(3) or 1), "n": int(m.group(4) or m.group(2)), "inv": int(m.group(5)),
                             "poly_s": float(m.group(6)), "root_s": float(m.group(7)) + float(m.group(8) or 0), "build_s": float(b.group(3))})
                pend = {}
Ds = sorted({r["D"] for r in rows})
out = subprocess.run(["gp", "-q"], input="".join("print(qfbclassno(-%d))\n" % d for d in Ds), capture_output=True, text=True).stdout.split()
H = dict(zip(Ds, map(int, out)))


def split(h):
    """Larger part of the most balanced factorisation h = a * b."""
    best = h
    for a in range(1, int(math.isqrt(h)) + 1):
        if h % a == 0:
            best = h // a
    return best


def lpf(h):
    p, x, big = 2, h, 1
    while p * p <= x:
        while x % p == 0:
            big, x = p, x // p
        p += 1
    return max(big, x) if x > 1 else big


for r in rows:
    r["h"] = H[r["D"]]
    r["bal"] = split(r["h"])
    r["lpf"] = lpf(r["h"])
json.dump(rows, open(os.path.join(HERE, "build_fit.json"), "w"), indent=0)
if __name__ == "__main__":
    big = [r for r in rows if r["bits"] >= 500]
    print("%d builds, %d at 500 bits and more; degree = h(D) in %d of them; n = balanced split in %d, n = largest prime factor in %d"
          % (len(rows), len(big), sum(r["deg"] == r["h"] for r in rows), sum(r["n"] == r["bal"] for r in rows), sum(r["n"] == r["lpf"] for r in rows)))
    # least squares for log(root_s) = a + b log n + c log bits on builds with a root step that can be timed
    import numpy as np
    sel = [r for r in big if r["root_s"] >= 0.3 and r["n"] >= 200]
    A = np.array([[1.0, math.log(r["n"]), math.log(r["bits"] / 768.0)] for r in sel])
    y = np.array([math.log(r["root_s"]) for r in sel])
    co, res, _, _ = np.linalg.lstsq(A, y, rcond=None)
    pred = A @ co
    print("root seconds ~ %.4f * n^%.2f * (bits/768)^%.2f on %d builds; rms log error %.2f (factor %.2f)"
          % (math.exp(co[0]), co[1], co[2], len(sel), float(np.sqrt(np.mean((y - pred) ** 2))), math.exp(float(np.sqrt(np.mean((y - pred) ** 2))))))
    for key in ("bal", "lpf"):
        A2 = np.array([[1.0, math.log(max(r[key], 2)), math.log(r["bits"] / 768.0)] for r in sel])
        c2 = np.linalg.lstsq(A2, y, rcond=None)[0]
        e = float(np.sqrt(np.mean((y - A2 @ c2) ** 2)))
        print("  from h alone via %s: %.4f * x^%.2f * (bits/768)^%.2f, rms log error %.2f (factor %.2f)" % (key, math.exp(c2[0]), c2[1], c2[2], e, math.exp(e)))
    print("largest builds:")
    for r in sorted(big, key=lambda r: -r["build_s"])[:14]:
        print("  %4d bits D=-%-13d h=%-7d deg %-7d = %5d x %-6d inv %-3d poly %6.1fs roots %6.1fs build %6.1fs  (bal %d, lpf %d)  %s/%s"
              % (r["bits"], r["D"], r["h"], r["deg"], r["m"], r["n"], r["inv"], r["poly_s"], r["root_s"], r["build_s"], r["bal"], r["lpf"], r["set"], r["file"]))
    for lo, hi in ((500, 700), (700, 800), (800, 1000), (1000, 1100)):
        b = sorted(r["build_s"] for r in big if lo <= r["bits"] < hi)
        if b:
            print("%d-%d bits: %d builds, median %.1f s, mean %.1f s, max %.1f s; over 20 s: %d" % (lo, hi, len(b), b[len(b) // 2], sum(b) / len(b), b[-1], sum(x > 20 for x in b)))
