#!/usr/bin/env python3
"""Fill index.template.html with the numbers computed from data/ and the model.

usage: make_report.py <pool_1030.pkl> <pool_768.pkl> [extra.json]
The two pools are ss.Pool pickles of the real orders (1030 bits, original
format, 2^32 strip; 768 bits, revised format, 2^26 strip); they are too large
to commit (see section 8 of the report for how they are produced).
"""
import json
import math
import os
import pickle
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "ecpp"))
import shortsim as ss  # noqa: E402

DATA = os.path.join(HERE, "data")
CALIB = os.path.join(DATA, "pi_calib.json")
HAVE_CALIB = os.path.exists(CALIB)
pool1030 = pickle.load(open(sys.argv[1], "rb"))
pool768 = pickle.load(open(sys.argv[2], "rb"))
EXTRA = json.load(open(sys.argv[3])) if len(sys.argv) > 3 else {}
em = ss.EcmModel(CALIB if HAVE_CALIB else None)
R = {}


def fmt(x, d=0):
    return ("{:,.%df}" % d).format(x)


def sweeps(spec):
    out = []
    for (m, b1, c, k) in spec:
        out += [(m, b1, c)] * k
    return out


# ---------------------------------------------------------------- 4.1 calibration vs the log
camp = [("ecm", 2000, 12), ("ecm", 11000, 20), ("ecm", 50000, 32), ("ecm", 250000, 48),
        ("ecm", 1000000, 72), ("ecm", 3000000, 96)]
caps = [0, 0, 12000, 4000, 1200, 400]
NPASS = 26958
obs = {"run": [25411, 24761, 12000, 4000, 1200, 400], "split": [7792, 6271, 2326, 683, 189, 59],
       "wall": [74.6, 259.4, 675.8, 1392.7, 2106.8, 2392.0]}
ccost = ss.make_cost(dict(ss.SEC1030_CAMPAIGN, **{}), exponent=0.0)
lad = ss.Ladder(camp, em, y0bits=32, cost=ccost)


def admit_ref(j, liv, tails):
    if not caps[j]:
        return liv
    idx = np.log(np.maximum(pool1030.lns[liv], 1.1)) - 0.03 * np.maximum(0, tails - pool1030.lam)
    return liv[np.argsort(-idx, kind="stable")[:caps[j]]]


rng = np.random.default_rng(1)
tot = None
NREP = 8
for rep in range(NREP):
    idx = np.arange(NPASS) if rep == 0 else rng.choice(len(pool1030.t), NPASS, replace=False)
    st, _ = ss.emulate(lad, pool1030, rng, admit_ref, idx=idx)
    arr = np.array([[s["run"], s["split"], s["dead"]] for s in st], dtype=float)
    tot = arr if tot is None else tot + arr
tot /= NREP
fr = [1, 1] + [c / NPASS for c in caps[2:]]
ev = ss.evaluate(lad, pool1030, ss.reference_policy(lad, pool1030, fr))
rows = ["<table><tr><th>Tier (B1, curves)</th><th>tails run: model</th><th>log</th>"
        "<th>tails split: model</th><th>log</th><th>expected winners</th><th>core-h in the log</th></tr>"]
for j, (m, b1, c) in enumerate(camp):
    rows.append("<tr><td>%s, %d</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%.3f</td><td>%.0f</td></tr>"
                % (fmt(b1), c, fmt(tot[j, 0]), fmt(obs["run"][j]), fmt(tot[j, 1]), fmt(obs["split"][j]),
                   ev["win"][j] * NPASS, obs["wall"][j] * 192 / 3600))
rows.append("</table>")
R["CALIB_TABLE"] = "\n".join(rows)
refwins = (ev["win"].sum() + ev["win0"]) * NPASS
R["REF_WINS"] = "%.2f, at %.0f core-hours of ECM" % (refwins, ev["cost"].sum() * NPASS / 3600)

# ---------------------------------------------------------------- 4.2 decomposition at 1030 bits
c2 = dict(ss.SEC1030_CAMPAIGN)
c2.update({10000000: 36.0, 30000000: 100.0})
cost = ss.make_cost(c2, exponent=1.5)
variants = [
    ("Campaign tiers and caps (the reference)", camp, "ref"),
    ("Campaign tiers, best stopping rule", camp, "opt"),
    ("Same B1, in sweeps of 12&ndash;20 curves", sweeps([("ecm", 2000, 12, 1), ("ecm", 11000, 20, 1), ("ecm", 50000, 16, 2),
                                                        ("ecm", 250000, 12, 4), ("ecm", 1000000, 12, 6), ("ecm", 3000000, 12, 8)]), "opt"),
    ("More sweeps at 10<sup>6</sup> and 3&middot;10<sup>6</sup>", sweeps([("ecm", 2000, 12, 1), ("ecm", 11000, 20, 1), ("ecm", 50000, 16, 2),
                                                                          ("ecm", 250000, 12, 4), ("ecm", 1000000, 12, 10), ("ecm", 3000000, 12, 24)]), "opt"),
    ("Adding B1 = 10<sup>7</sup> and 3&middot;10<sup>7</sup> (cost and &pi; extrapolated)",
     sweeps([("ecm", 2000, 12, 1), ("ecm", 11000, 20, 1), ("ecm", 50000, 16, 2), ("ecm", 250000, 12, 4),
             ("ecm", 1000000, 12, 6), ("ecm", 3000000, 12, 8), ("ecm", 10000000, 12, 8), ("ecm", 30000000, 12, 8)]), "opt"),
]
prices = (8.8, 1.0, 0.3)
rows = ["<table><tr><th>Schedule</th>" + "".join("<th>supply at %.1f core-s per order</th>" % p for p in prices) + "</tr>"]
keep = {}
for name, tiers, kind in variants:
    ladv = ss.Ladder(tiers, em, y0bits=32, cost=cost)
    cells = []
    for price in prices:
        if kind == "ref":
            pol = ss.reference_policy(ladv, pool1030, fr)
        else:
            _, pol = ss.optimal(ladv, pool1030, price)
        e = ss.evaluate(ladv, pool1030, pol)
        pw = e["win"].sum() + e["win0"]
        ce = e["cost"].sum()
        cells.append("<td>%.0f <span class=\"mut\">(%s orders, %.1f core-s of ECM each)</span></td>"
                     % ((price + ce) / pw / 3600, fmt(1 / pw), ce))
        keep[(name, price)] = (e, tiers, pw, ce)
    rows.append("<tr><td>%s</td>%s</tr>" % (name, "".join(cells)))
rows.append("</table>")
R["DECOMP_TABLE"] = "\n".join(rows) + (
    "\n<p class=\"mut\">8.8 core-s per order is what supply cost in the August 1030-bit run; 1.0 and 0.3 bracket "
    "what it costs now (scan 0.1&ndash;0.25, strip under 0.01, plus overheads).  In brackets: orders bought per "
    "winner and ECM spent per order.</p>")
e, tiers, pw, ce = keep[(variants[2][0], 1.0)]
b1s = []
for t in tiers:
    if t[1] not in b1s:
        b1s.append(t[1])
frac = {b1: float(np.mean([e["run"][j] for j, t in enumerate(tiers) if t[1] == b1])) for b1 in b1s}
R["SELECT_SENTENCE"] = ("at 1.0 core-s per order, %.0f%% of fresh orders get the B1&nbsp;=&nbsp;2,000 sweep, %.0f%% reach 11,000, "
                        "%.1f%% reach 50,000 and %.1f%% reach 250,000"
                        % (100 * frac[2000], 100 * frac[11000], 100 * frac[50000], 100 * frac[250000]))
wshare = {b1: sum(e["win"][j] for j, t in enumerate(tiers) if t[1] == b1) / pw for b1 in b1s}
R["DEPTH_SENTENCE"] = ("in the third row of the table %.0f%% of the winners come from the 10<sup>6</sup> sweeps and %.0f%% from "
                       "3&middot;10<sup>6</sup>" % (100 * wshare[1000000], 100 * wshare[3000000]))

# ---------------------------------------------------------------- 5 experiment
pol = json.load(open(os.path.join(DATA, "exp768_policy.json")))
run = json.load(open(os.path.join(DATA, "exp768_policy_run.json")))
tiers = [tuple(t) for t in pol["tiers"]]
sec = {int(k): v for k, v in pol["sec_per_curve"].items()}


def cost768(b1, tb):
    return sec[b1] * (np.asarray(tb, dtype=float) / 768.0) ** 1.5


go = {}
for c, per in pol["go"].items():
    arrs = []
    for iv in per:
        a = np.zeros(ss.NB, dtype=bool)
        for lo, hi in iv:
            a[lo:hi] = True
        arrs.append(a)
    go[int(c)] = arrs
N = len(pool768.t)
models = {"theory": pol["model"]}
if HAVE_CALIB:
    lad768 = ss.Ladder(tiers, em, y0bits=26, cost=cost768)
    e7 = ss.evaluate(lad768, pool768, go)
    models["measured"] = {"run": (e7["run"] * N).tolist(), "win": (e7["win"] * N).tolist(),
                          "cost": (e7["cost"] * N).tolist(), "win0": e7["win0"] * N}
st = {s["tier"]: s for s in run["stats"]}
b1s = []
for t in tiers:
    if t[1] not in b1s:
        b1s.append(t[1])
hdr = "<table><tr><th>B1</th><th>sweeps</th><th>tails run: predicted</th><th>observed</th><th>winners: predicted</th>"
if HAVE_CALIB:
    hdr += "<th>re-predicted with measured &pi;</th>"
hdr += "<th>observed</th><th>ECM core-s: predicted</th><th>observed</th></tr>"
rows = [hdr]
T = {"rm": 0, "ro": 0, "wm": 0, "wm2": 0, "wo": 0, "cm": 0, "co": 0}
for b1 in b1s:
    js = [j for j, t in enumerate(tiers) if t[1] == b1]
    rm = sum(models["theory"]["run"][j] for j in js)
    ro = sum(st.get(j, {}).get("run", 0) for j in js)
    wm = sum(models["theory"]["win"][j] for j in js)
    wm2 = sum(models["measured"]["win"][j] for j in js) if HAVE_CALIB else 0
    wo = sum(st.get(j, {}).get("wins", 0) for j in js)
    cm = sum(models["theory"]["cost"][j] for j in js)
    co = sum(st.get(j, {}).get("cpu_s", 0) for j in js)
    if ro == 0 and rm < 0.5:
        continue
    for k, v in (("rm", rm), ("ro", ro), ("wm", wm), ("wm2", wm2), ("wo", wo), ("cm", cm), ("co", co)):
        T[k] += v
    rows.append("<tr><td>%s</td><td>%d</td><td>%s</td><td>%s</td><td>%.2f</td>%s<td>%d</td><td>%s</td><td>%s</td></tr>"
                % (fmt(b1), len(js), fmt(rm), fmt(ro), wm, ("<td>%.2f</td>" % wm2) if HAVE_CALIB else "", wo, fmt(cm), fmt(co)))
rows.append("<tr><th>Total</th><th></th><th>%s</th><th>%s</th><th>%.1f</th>%s<th>%d</th><th>%s (%.1f core-h)</th><th>%s (%.1f core-h)</th></tr>"
            % (fmt(T["rm"]), fmt(T["ro"]), T["wm"], ("<th>%.1f</th>" % T["wm2"]) if HAVE_CALIB else "", T["wo"],
               fmt(T["cm"]), T["cm"] / 3600, fmt(T["co"]), T["co"] / 3600))
rows.append("</table>")
R["EXP_TABLE"] = "\n".join(rows)
R["EXP_N"] = fmt(N)
R["EXP_WINS"] = "%d" % T["wo"]
R["EXP_COREH"] = "%.1f" % (T["co"] / 3600)
R["EXP_PER_WIN"] = "%.1f" % (T["co"] / 3600 / T["wo"])
R["EXP_MODEL_WINS"] = "%.1f" % T["wm"]
R["EXP_MODEL_COREH"] = "%.1f" % (T["cm"] / 3600)
R["EXP_MODEL_WINS2"] = ("%.1f" % T["wm2"]) if HAVE_CALIB else ""
R["EXP_WALL"] = fmt(run["elapsed"])
R["PI_SOURCE"] = ("measured for gmp-ecm on the box (section 8, <code>data/pi_calib.json</code>) where the measurement has "
                  "at least 8 successes, and otherwise the two-stage Dickman function for orders near 2<sup>b</sup>/24, "
                  "rescaled to join the measured cells") if HAVE_CALIB else (
    "the two-stage Dickman function for a curve order near 2<sup>b</sup>/24 with gmp-ecm's default B2")
R["PI_FILE"] = ", result <code>data/pi_calib.json</code>" if HAVE_CALIB else " (not yet run to completion)"
R.update(EXTRA)

tpl = open(os.path.join(HERE, "index.template.html")).read()
body = open(os.path.join(HERE, "sections.html")).read()
# sections.html holds the prose blocks, each introduced by a line "#### NAME"
blocks = {}
cur = None
for line in body.splitlines():
    if line.startswith("#### "):
        cur = line[5:].strip()
        blocks[cur] = []
    elif cur:
        blocks[cur].append(line)
for k, v in blocks.items():
    R.setdefault(k, "\n".join(v))
for _ in range(3):                       # blocks may contain tokens themselves
    for k, v in R.items():
        tpl = tpl.replace("@@%s@@" % k, str(v))
left = sorted(set(x for x in __import__("re").findall(r"@@([A-Z0-9_]+)@@", tpl)))
if left:
    print("unfilled tokens:", left)
open(os.path.join(HERE, "index.html"), "w").write(tpl)
json.dump({k: v for k, v in R.items() if len(str(v)) < 200}, open(os.path.join(DATA, "report_numbers.json"), "w"), indent=1)
print("index.html written;", "measured pi" if HAVE_CALIB else "theory pi")
