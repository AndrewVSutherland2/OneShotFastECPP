#!/usr/bin/env python3
"""Scan small discriminants for decompositions whose V has a repeated root (s = e_1 fails)."""
import os, subprocess, sys, re
lo, hi, inv = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
S = os.path.dirname(os.path.abspath(__file__))
gp = "forstep(d=%d,%d,1,if((d%%4==0||d%%4==3)&&d>4,my(D=-d);if(isfundamental(D)||1,my(h=qfbclassno(D));if(h>=4&&!isprime(h),print(D,\" \",h)))));\n" % (lo, hi)
out = subprocess.run(["gp", "-q"], input=gp, capture_output=True, text=True).stdout.split("\n")
rev = tot = 0
for ln in out:
    if not ln.strip(): continue
    D, h = map(int, ln.split())
    f = os.path.join(S, "s_%d_%s.txt" % (os.getpid(), inv))
    for n in [x for x in range(2, h) if h % x == 0]:
        e = dict(os.environ); e["CLASSPOLY_TOWER_N"] = str(n)
        r = subprocess.run(["classpoly", str(D), inv, "0", f, "0"], capture_output=True, text=True, env=e, cwd=S)
        o = r.stdout + r.stderr
        if "Decomposition:" in o: tot += 1
        if "repeated coset value" in o:
            rev += 1
            print("REPEATED D=%d h=%d n=%d inv=%s: %s" % (D, h, n, inv, [l for l in o.splitlines() if "repeated" in l][0][:160]), flush=True)
    if os.path.exists(f): os.remove(f)
print("scan %d..%d inv=%s: %d decompositions, %d reverted" % (lo, hi, inv, tot, rev))
