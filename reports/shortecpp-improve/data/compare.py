#!/usr/bin/env python3
"""Observed (policy_run.json) against the model's prediction for the same policy.
usage: compare.py <policy.json> <policy_run.json> [calib.json]   (with calib: re-predict using the measured pi)"""
import sys, os, json, math, pickle
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', 'ecpp'))
import shortsim as ss
pol = json.load(open(sys.argv[1])); run = json.load(open(sys.argv[2]))
calib = sys.argv[3] if len(sys.argv) > 3 else None
tiers = [tuple(t) for t in pol["tiers"]]
model = pol["model"]
if calib:
    pool = pickle.load(open(sys.argv[1] + ".pool.pkl", "rb"))
    sec = {int(k): v for k, v in pol["sec_per_curve"].items()}
    cost = lambda b1, tb: sec[b1] * (np.asarray(tb, dtype=float) / 768.0) ** 1.5
    lad = ss.Ladder(tiers, ss.EcmModel(calib), y0bits=26, cost=cost)
    go = {}
    for c, per in pol["go"].items():
        arrs = []
        for iv in per:
            a = np.zeros(ss.NB, dtype=bool)
            for lo, hi in iv: a[lo:hi] = True
            arrs.append(a)
        go[int(c)] = arrs
    ev = ss.evaluate(lad, pool, go)
    N = len(pool.t)
    model = {"run": (ev["run"] * N).tolist(), "win": (ev["win"] * N).tolist(), "cost": (ev["cost"] * N).tolist(), "win0": ev["win0"] * N}
st = {s["tier"]: s for s in run["stats"]}
b1s = []
for t in tiers:
    if t[1] not in b1s: b1s.append(t[1])
print("B1         sweeps   tails run (model/obs)     wins (model/obs)    core-s (model/obs)")
TM = TO = CM = CO = 0.0
rows = []
for b1 in b1s:
    js = [j for j, t in enumerate(tiers) if t[1] == b1]
    rm = sum(model["run"][j] for j in js); ro = sum(st.get(j, {}).get("run", 0) for j in js)
    wm = sum(model["win"][j] for j in js); wo = sum(st.get(j, {}).get("wins", 0) for j in js)
    cm = sum(model["cost"][j] for j in js); co = sum(st.get(j, {}).get("cpu_s", 0) for j in js)
    done = all(j in st for j in js)
    TM += wm; TO += wo; CM += cm; CO += co
    rows.append((b1, len(js), rm, ro, wm, wo, cm, co))
    print("%-9d  %3d   %10.0f / %-8d   %7.2f / %-3d   %10.0f / %-8.0f %s" % (b1, len(js), rm, ro, wm, wo, cm, co, "" if done else "(incomplete)"))
w0 = st.get(-1, {}).get("wins", 0)
print("strip: wins model %.3f obs %d" % (model.get("win0", 0), w0))
print("TOTAL wins: model %.2f, observed %d;  ECM core-h: model %.1f, observed %.1f" % (TM + model.get("win0", 0), TO + w0, CM / 3600, CO / 3600))
if TO + w0:
    print("observed ECM core-h per winner: %.2f   (model %.2f)" % (CO / 3600 / (TO + w0), CM / 3600 / (TM + model.get("win0", 0))))
json.dump({"rows": rows, "model_wins": TM + model.get("win0", 0), "obs_wins": TO + w0, "model_core_h": CM / 3600, "obs_core_h": CO / 3600}, open(sys.argv[2] + ".cmp.json", "w"))
