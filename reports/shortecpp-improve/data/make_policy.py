#!/usr/bin/env python3
"""Build the stopping policy for the 768-bit experiment and print what the model
expects of it.
usage: make_policy.py <p.txt> <compact.txt> <fac> <price> <out policy.json> [calib.json]
<compact.txt>: "log2(tail) prp s" per order; <fac>: factorisation of each s."""
import sys, os, math, json, time, pickle, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', 'ecpp'))
import shortsim as ss
pfile, compact, facf, price, outp = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4]), sys.argv[5]
calib = sys.argv[6] if len(sys.argv) > 6 else None
p = int(open(pfile).read()); n = p.bit_length()
rt = math.isqrt(p); L = rt + 1 + math.isqrt(4 * rt); lam = math.log2(L)
bcap = math.ceil(n * n / math.log2(n)); radlim = math.ceil(n / math.log2(n))
t0 = time.time()
tb, cov, prp, lns = [], [], [], []
for line, fl in zip(open(compact), open(facf)):
    w = line.split()
    sfac = {int(a): int(b) for a, b in (tok.split("^") for tok in fl.split())}
    tb.append(float(w[0])); prp.append(w[1] == "1"); lns.append(math.log(max(int(w[2]), 1)))
    cov.append(ss.coverage_cells(sfac, bcap=bcap, radlim=radlim))
pool = ss.Pool(tb, np.array(cov), lam, prp, lns)
print("%d orders loaded in %.0fs; lam=%.1f; mean coverage %.1f bits; median tail %.0f bits; prime tails %.2f%%"
      % (len(tb), time.time() - t0, lam, pool.w.mean(), np.median(pool.t), 100 * pool.prp.mean()))
# loaded core-seconds per curve on a 768-bit modulus (c8a.48xlarge, gmp-ecm 7.0.5, default B2)
sec768 = {2000: 0.0089, 11000: 0.0422, 50000: 0.172, 250000: 0.772, 1000000: 3.16, 3000000: 8.64,
          10000000: 27.0, 30000000: 80.0}
def cost(b1, tbits):
    return sec768[b1] * (np.asarray(tbits, dtype=float) / 768.0) ** 1.5
spec = [("ecm", 2000, 12, 1), ("ecm", 11000, 20, 1), ("ecm", 50000, 16, 2), ("ecm", 250000, 12, 4),
        ("ecm", 1000000, 12, 6), ("ecm", 3000000, 12, 8), ("ecm", 10000000, 12, 8), ("ecm", 30000000, 12, 4)]
tiers = []
for (m, b1, c, k) in spec:
    tiers += [(m, b1, c)] * k
em = ss.EcmModel(calib)
lad = ss.Ladder(tiers, em, y0bits=26, cost=cost)
cstar, pol = ss.optimal(lad, pool, price)
ev = ss.evaluate(lad, pool, pol)
N = len(tb)
pw = ev['win'].sum() + ev['win0']; ce = ev['cost'].sum()
print("price %.2f core-s/order: model cost per win %.1f core-h; P(win)/order %.3e (%.0f orders/win); ECM %.2f core-s/order"
      % (price, cstar / 3600, pw, 1 / pw, ce))
print("expected for these %d orders: %.2f wins (%.2f at the strip), %.1f core-h of ECM" % (N, pw * N, ev['win0'] * N, ce * N / 3600))
print("tier  B1        c   E[run]   E[wins]  E[core-s]")
for j, (m, b1, c) in enumerate(tiers):
    print("%3d  %-9d %2d  %8.1f  %7.4f  %9.0f" % (j, b1, c, ev['run'][j] * N, ev['win'][j] * N, ev['cost'][j] * N))
go = {}
for c, g in pol.items():
    per = []
    for j in range(len(tiers)):
        a = g[j].astype(int)
        d = np.diff(np.concatenate([[0], a, [0]]))
        per.append([[int(x), int(y)] for x, y in zip(np.nonzero(d == 1)[0], np.nonzero(d == -1)[0])])
    go[int(c)] = per
json.dump({"tiers": tiers, "cls": [int(x) for x in pool.cls], "go": go, "price": price,
           "sec_per_curve": {str(k): v for k, v in sec768.items()},
           "model": {"run": (ev['run'] * N).tolist(), "win": (ev['win'] * N).tolist(), "cost": (ev['cost'] * N).tolist(),
                     "win0": ev['win0'] * N, "cstar_core_h": cstar / 3600}}, open(outp, "w"))
pickle.dump(pool, open(outp + ".pool.pkl", "wb"))
print("policy written to", outp)
