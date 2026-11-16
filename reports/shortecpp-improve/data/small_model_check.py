"""Model check at a small size: work every live order of a scan through the whole ladder, sweep by
sweep, and count the winners; compare with the model's win probabilities summed over the orders, in
total, per sweep, by coverage class, by tail size, and among the candidates the rule ranks best.

Every sweep is run in full (8 curves at a time, on what is left of the tail), which is what the
scheduler's state "j sweeps done" means.

usage: small_model_check.py yafu=<binary> [c=90] [B=8e6] [threads=24]     (c: the prime is nextprime(10^c))"""
import sys, os, math, concurrent.futures as cf
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "ecpp"))
import numpy as np
import shortECPP as s
import shortsim as ss
import shortsched as sd
from math import isqrt, log2, ceil

kv = dict(a.split("=", 1) for a in sys.argv[1:])
c10 = int(kv.get("c", 90)); B = int(float(kv.get("B", 8e6))); thr = int(kv.get("threads", 24))
s.YAFU_BIN = kv["yafu"]; s.ENGINE = "yafu"
p = int(s.gp("print(nextprime(10^%d))" % c10).split()[-1]) if "p" not in kv else int(kv["p"])
n = p.bit_length(); nb = n
s.V2_BCAP = ceil(n * n / log2(n)); s.V2_RADLIM = ceil(n / log2(n)); s.V2_QFLOOR = s.V2_BCAP ** 2
n2 = n * n
small = s.sieve_primes(n2)
rt = isqrt(p); L = rt + 1 + isqrt(4 * rt); lam = log2(L)
y0b, B0, _ = s.level_params(nb, thr, s.V2_BCAP)
if "y0" in kv: y0b = int(kv["y0"])
s.set_strip_bound(y0b)
recs = s.run_dscan(p, B, 0, thr, s.fb_cap(B))
cands = [s.Cand(d, t, v, N) for (d, t, v) in recs for N in (p + 1 - t, p + 1 + t) if N % 4 == 0]
strip = s.smooth_strip(sorted(set(c.N for c in cands)), 8, s.V2_BCAP)
ss.set_grid(nb)
engine = "yafu"
em, sec, mcost, quantum = ss.engine_model(engine, 1.0, s.ECPP_DIR)
tiers = sd.default_tiers(engine, nb)
if "top" in kv: tiers = [t for t in sd.default_tiers(engine) if t[1] <= int(float(kv["top"]))]
lad = ss.Ladder(tiers, em, y0bits=y0b, cost=mcost)
scale = min(1.0, (nb / 500.0) ** 8)
tab = ss.IndexTable(lad, lam, rlo=min(30.0, 3.0 * scale), rhi=3.0e7 * max(scale, 1e-3), threads=8)
print("p %d bits, D=%d, y0=2^%d, %d orders, tiers %s" % (nb, ss.D, y0b, len(cands), sorted(set((b, c) for (_, b, c) in tiers))))
live = []; win0 = 0; TB0 = {}
winners = set()
def accept(c, f):
    if f * c.s <= L: return False
    return bool(s.accept_q(c, f, p, L, n2, small))
for c in cands:
    S, c.tail, c.s, prp = strip[c.N]
    if c.tail == 1 or c.s * c.tail <= L: continue
    if prp:
        if c.tail * c.tail < p and s.is_prp(c.tail, 24) and accept(c, c.tail):
            win0 += 1
        continue
    W = log2(c.s) if c.s > 1 else 0.0
    if c.s.bit_length() > s.V2_RADLIM:
        W = float(ss.coverage_cells(s.factor_small(c.s, small), bcap=s.V2_BCAP, radlim=s.V2_RADLIM).sum())
    c.cls = int(round(W)); c.tb = int(round(log2(c.tail))); TB0[c.N] = c.tb
    live.append(c)
# model: probability of a win somewhere in the ladder = value at the top reward / reward
g = len(tab.rgrid) - 1
pm = np.array([tab.V0[g][min(c.cls, tab.cmax), min(c.tb, ss.NB - 1)] / tab.rgrid[g] for c in live], dtype=float)
# the same by summing the per-sweep win probabilities along the model's path is not available in closed form;
# V0 at the top reward counts wins minus cost/reward, so add back the expected cost/reward
cm = np.array([tab.C0[g][min(c.cls, tab.cmax), min(c.tb, ss.NB - 1)] / tab.rgrid[g] for c in live], dtype=float)
print("intake winners %d; %d live composite; model: expected winners with the whole ladder %.2f (+cost term %.3f)"
      % (win0, len(live), pm.sum(), cm.sum()))
# first-sweep win probabilities
p1 = sum(tab.pwin(c.cls, 0, c.tb) for c in live)
print("model: expected wins of the first sweep alone %.2f" % p1)
runner = s.EcmRunner("yafu")
alive = list(live)
PW0 = {c.N: tab.pwin(c.cls, 0, c.tb) for c in live}
R0 = {c.N: tab.r(c.cls, 0, c.tb) for c in live}
CL = {c.N: (c.cls, c.tb) for c in live}
WONAT = {}
for j, (m, b1, curves) in enumerate(tiers):
    exp = sum(tab.pwin(c.cls, j, c.tb) for c in alive)
    n0 = len(alive)
    nw = 0
    for rep in range(curves // 8):
        k = max(1, min(32, int(2.0 / max(mcost(b1, nb) * 8, 1e-4))))
        units = [alive[i:i + k] for i in range(0, len(alive), k)]
        with cf.ThreadPoolExecutor(max_workers=thr) as ex:
            for u, res in zip(units, ex.map(lambda u: runner.run([c.tail for c in u], b1, 8, 600), units)):
                for c in u:
                    if c.tail not in res: continue
                    found, cof, cprp = res[c.tail]
                    won = False
                    for f in found:
                        if f > 1 and c.tail % f == 0:
                            c.tail //= f
                            if s.is_prp(f, 12) and accept(c, f): won = True
                    if c.tail > 1 and c.tail * c.tail < p and s.is_prp(c.tail, 24) and accept(c, c.tail): won = True
                    if won: winners.add(c.N); nw += 1; c.dead = True; WONAT[c.N] = j
                    elif c.tail == 1 or c.s * c.tail <= L or s.is_prp(c.tail, 12): c.dead = True
        alive = [c for c in alive if not c.dead]
    for c in alive:
        c.tb = int(round(log2(max(c.tail, 2))))
    print("sweep %d B1=%d x %d on %d: %d winners (model expected %.2f for this sweep), %d left" % (j, b1, curves, n0, nw, exp, len(alive)), flush=True)
print("TOTAL winners %d of %d live (model %.2f)" % (len(winners), len(live), pm.sum() + cm.sum()))
# breakdown by class
import collections
bc = collections.defaultdict(lambda: [0, 0, 0.0])
for c, a, b in zip(live, pm, cm):
    k = min(c.cls // 5 * 5, 40)
    bc[k][0] += 1; bc[k][1] += (c.N in winners); bc[k][2] += a + b
for k in sorted(bc): print("class %2d-%2d: %5d orders, %3d winners, model %.2f" % (k, k + 4, bc[k][0], bc[k][1], bc[k][2]))
bt = collections.defaultdict(lambda: [0, 0, 0.0])
for c, a, b in zip(live, pm, cm):
    k = (TB0[c.N] - int(lam)) // 20 * 20
    bt[k][0] += 1; bt[k][1] += (c.N in winners); bt[k][2] += a + b
for k in sorted(bt): print("tail bits - lam in [%3d,%3d): %5d orders, %3d winners, model %.2f" % (k, k + 20, bt[k][0], bt[k][1], bt[k][2]))

order = sorted(live, key=lambda c: R0[c.N])
for k in (50, 100, 200, 400, 800, 1600, len(order)):
    top = order[:k]
    print("best %4d by r: first-sweep wins expected %.2f, got %d; won in any sweep %d; r up to %.2f" % (len(top), sum(PW0[c.N] for c in top), sum(1 for c in top if WONAT.get(c.N) == 0), sum(1 for c in top if c.N in WONAT), R0[top[-1].N]))
print("top 30 by r: (class, tail bits - lam, pwin0, r, won at)")
for c in order[:30]:
    print("  cls %2d  tb-lam %4d  pw0 %.3f  r %.3f  %s" % (CL[c.N][0], CL[c.N][1] - int(lam), PW0[c.N], R0[c.N], WONAT.get(c.N, "-")))
