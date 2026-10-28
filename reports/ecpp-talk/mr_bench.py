#!/usr/bin/env python3
"""mr_bench.py -- verification cost of the repo's certificates in Miller-Rabin-round
equivalents (the "Cost to verify" table of talk.tex).  Both sides use CPython's own
big-integer arithmetic (the verifiers' ladders vs. pow()), so the time ratio tracks the
modular-multiplication ratio.  Run from anywhere:  python3 reports/ecpp-talk/mr_bench.py
"""
import os, sys, time, csv
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "ecpp")); sys.path.insert(0, ROOT)
import vshort2, voneshot
cnt = {'mul': 0}
def counted(mod):
    od, oa = mod.xdbl, mod.xadd
    def xd(*a):
        cnt['mul'] += 6; return od(*a)
    def xa(*a):
        cnt['mul'] += 6; return oa(*a)
    mod.xdbl, mod.xadd = xd, xa
def mr_time(p, reps=20):
    d = p - 1; s = 0
    while d % 2 == 0: d //= 2; s += 1
    t = time.perf_counter()
    for a in range(2, 2 + reps):
        x = pow(a, d, p)
        for _ in range(s - 1):
            if x == 1 or x == p - 1: break
            x = x * x % p
    return (time.perf_counter() - t) / reps
def best(f, reps=5):
    b = 1e9
    for _ in range(reps):
        t = time.perf_counter(); ok = f(); b = min(b, time.perf_counter() - t); assert ok
    return b
rows = list(csv.reader(open(os.path.join(ROOT, 'certs/short2/certs.csv'))))
print("short ECPP (vshort2.py, revised format)")
for r in rows:
    seq = [int(v) for v in r if v.strip()]
    p = seq[0]
    if p.bit_length() not in (34, 266, 333, 499, 665, 798, 1030, 1329): continue
    tv = best(lambda: vshort2.verify(seq)); tm = mr_time(p)
    od, oa = vshort2.xdbl, vshort2.xadd
    counted(vshort2); cnt['mul'] = 0; vshort2.verify(seq); muls = cnt['mul']
    vshort2.xdbl, vshort2.xadd = od, oa
    print(f"  {p.bit_length():5d} bits  levels={(len(seq)-1)//3}  verify={tv*1e3:8.2f} ms  MR round={tm*1e3:7.3f} ms  ratio={tv/tm:6.1f}  ladder-mults={muls} = {muls/p.bit_length():.1f} n")
print("one-shot (voneshot.py)")
for f in ['25519', 'p256', 'k256', '1e100p267']:
    toks = open(os.path.join(ROOT, f'certs/{f}.txt')).read().split()
    vals = [int(t) for t in toks]
    p = vals[0]
    tv = best(lambda: voneshot.verify(vals[0], vals[1], vals[2], vals[3], tuple(vals[4:])))
    tm = mr_time(p)
    print(f"  {f:10s} {p.bit_length():4d} bits verify={tv*1e3:7.2f} ms  MR={tm*1e3:6.3f} ms  ratio={tv/tm:5.1f}")
