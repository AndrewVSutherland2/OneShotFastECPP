#!/usr/bin/env python3
"""Fill index.template.html with the numbers computed from data/ and the model.

usage: make_report.py <pool_1030.pkl> <pool_768.pkl> [extra.json]
The two pools are ss.Pool pickles of the real orders (1030 bits, original
format, 2^32 strip; 768 bits, revised format, 2^26 strip); they are too large
to commit (see section 8 of the report for how they are produced).
"""
import json
import re
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

# ---- sections 7 and 8: class-group decomposition timings and the yafu ECM benchmark (data/tower_*, data/yafu_*)
def cm_rows(fn):
    """cm_method timing log -> {|D|: dict}"""
    out, cur, retry = {}, None, 0
    if not os.path.exists(fn):
        return out
    for ln in open(fn):
        if "retrying" in ln:
            retry = 1
        m = re.match(r"\[H_D\^inv\(deg (\d+) = (\d+) x (\d+), inv (\d+)\): ([\d.]+)s\s+roots: ([\d.]+)s \+ ([\d.]+)s", ln)
        if m:
            cur = {"deg": int(m.group(1)), "m": int(m.group(2)), "n": int(m.group(3)), "inv": int(m.group(4)),
                   "crt": float(m.group(5)), "root": float(m.group(6)) + float(m.group(7))}
        m = re.match(r"\[H_D\^inv\(deg (\d+), inv (\d+)\): ([\d.]+)s\s+root: ([\d.]+)s", ln)
        if m:
            cur = {"deg": int(m.group(1)), "m": 0, "n": 0, "inv": int(m.group(2)), "crt": float(m.group(3)), "root": float(m.group(4))}
        m = re.match(r"D=(\d+) wall=([\d.]+) user=([\d.]+)", ln)
        if m and float(m.group(2)) > 0 and cur:
            cur.update({"wall": float(m.group(2)), "user": float(m.group(3)), "retry": retry})
            out[int(m.group(1))] = cur
            cur, retry = None, 0
    return out

INVNAME = {1: "f", 14: "w<sub>2,7</sub>", 26: "w<sub>2,13</sub>", 129: "A<sub>29</sub>", 147: "A<sub>47</sub>", 159: "A<sub>59</sub>",
           171: "A<sub>71</sub>", 521: "w<sub>3,7</sub>", 535: "w<sub>5,7</sub>"}
tw4 = cm_rows(os.path.join(DATA, "tower_cm_192.log"))
tw3 = cm_rows(os.path.join(DATA, "tower_cm_192_pass3.log"))
plain = cm_rows(os.path.join(DATA, "p768_cm_timing.log"))
if tw4:
    rows = []
    for d in sorted(tw4):
        r, b = tw4[d], plain.get(d)
        rows.append("<tr><td>&minus;%s</td><td>%s</td><td>%s</td><td>%d &times; %d</td><td>%.1f</td><td>%.1f</td><td><b>%.1f</b>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>" % (
            fmt(d), fmt(r["deg"]), INVNAME.get(r["inv"], str(r["inv"])), r["m"], r["n"], r["crt"], r["root"], r["wall"],
            " *" if r["retry"] else "", fmt(r["user"]), ("%.0f" % b["wall"]) if b else "", fmt(b["user"]) if b else ""))
    R["TOWER_TABLE"] = ("<table><tr><th>D</th><th>degree</th><th>invariant</th><th>m &times; n</th><th>CRT, s</th><th>roots, s</th>"
                        "<th>wall, s</th><th>core-s</th><th>before: wall, s</th><th>core-s</th></tr>" + "".join(rows) + "</table>")
    walls = sorted(r["wall"] for r in tw4.values())
    R["TOWER_N"] = len(tw4)
    R["TOWER_WALL_TOTAL"] = fmt(sum(walls))
    R["TOWER_WALL_MEDIAN"] = "%.1f" % walls[len(walls) // 2]
    R["TOWER_WALL_MAX"] = "%.0f" % walls[-1]
    R["TOWER_COREH"] = "%.1f" % (sum(r["user"] for r in tw4.values()) / 3600)
    R["TOWER_COREH_EACH"] = "%.2f" % (sum(r["user"] for r in tw4.values()) / 3600 / len(tw4))
    both = [d for d in tw4 if d in plain]
    R["TOWER_BEFORE_WALL"] = fmt(sum(plain[d]["wall"] for d in both))
    R["TOWER_AFTER_WALL"] = fmt(sum(tw4[d]["wall"] for d in both))
    R["TOWER_BEFORE_COREH"] = "%.0f" % (sum(plain[d]["user"] for d in both) / 3600)
    R["TOWER_AFTER_COREH"] = "%.1f" % (sum(tw4[d]["user"] for d in both) / 3600)
    R["TOWER_RETRIES"] = sum(r["retry"] for r in tw4.values())
    if tw3:
        common = [d for d in tw4 if d in tw3]
        R["TOWER_ROOT_SERIAL"] = fmt(sum(tw3[d]["root"] for d in common))
        R["TOWER_ROOT_PAR"] = fmt(sum(tw4[d]["root"] for d in common))
yt = os.path.join(DATA, "yafu_test.json")
if os.path.exists(yt):
    res = json.load(open(yt))["res"]
    def med(v):
        v = sorted(v)
        return v[len(v) // 2]
    rows = []
    for b1 in sorted(set(r["b1"] for r in res)):
        g = med([r["secs"] / 8 for r in res if r["kind"] == "cost" and r["eng"] == "gmp" and r["b1"] == b1])
        y = med([r["secs"] / 8 for r in res if r["kind"] == "cost" and r["eng"] == "yafu" and r["b1"] == b1])
        rows.append("<tr><td>%s</td><td>%.3f</td><td>%.3f</td><td>%.1f</td></tr>" % (fmt(b1), g, y, g / y))
    R["YAFU_COST_TABLE"] = ("<table><tr><th>B1</th><th>gmp-ecm 7.0.5, core-s per curve</th><th>yafu, core-s per curve</th><th>ratio</th></tr>"
                            + "".join(rows) + "</table>")
    rows, tg, ty, tn = [], 0, 0, 0
    for b1 in sorted(set(r["b1"] for r in res)):
        for b in sorted(set(r["bits"] for r in res if r["kind"] == "success" and r["b1"] == b1)):
            cell = [r for r in res if r["kind"] == "success" and r["b1"] == b1 and r["bits"] == b]
            g = sum(r["found"] for r in cell if r["eng"] == "gmp")
            y = sum(r["found"] for r in cell if r["eng"] == "yafu")
            n = len(cell) // 2
            tg, ty, tn = tg + g, ty + y, tn + n
            rows.append("<tr><td>%s</td><td>%d</td><td>%d of %d</td><td>%d of %d</td></tr>" % (fmt(b1), b, g, n, y, n))
    R["YAFU_SUCCESS_TABLE"] = ("<table><tr><th>B1</th><th>factor, bits</th><th>gmp-ecm found</th><th>yafu found</th></tr>" + "".join(rows) + "</table>")
    R["YAFU_FOUND_GMP"], R["YAFU_FOUND_YAFU"], R["YAFU_FOUND_N"] = tg, ty, tn
    R["YAFU_BAD"] = sum(r["bad"] for r in res)
yc = os.path.join(DATA, "yafu_cost.json")
if os.path.exists(yc):
    res = json.load(open(yc))
    rows = []
    for nb in sorted(set(r["bits"] for r in res)):
        for b1 in sorted(set(r["b1"] for r in res if r["bits"] == nb)):
            c = {(r["eng"], r["mode"]): r for r in res if r["bits"] == nb and r["b1"] == b1}
            if ("gmp", "") in c and ("yafu", "1") in c and ("yafu", "2") in c:
                g, y1, y2 = c[("gmp", "")]["sec_per_curve"], c[("yafu", "1")]["sec_per_curve"], c[("yafu", "2")]["sec_per_curve"]
                rows.append("<tr><td>%d</td><td>%s</td><td>%.4f</td><td>%.4f</td><td>%.1f</td><td>%.4f</td></tr>" % (nb, fmt(b1), g, y1, g / y1, y2))
    R["YAFU_GRID_TABLE"] = ("<table><tr><th>input, bits</th><th>B1</th><th>gmp-ecm</th><th>yafu, B2 as gmp-ecm</th><th>ratio</th>"
                            "<th>yafu, B2 from memory limit</th></tr>" + "".join(rows) + "</table>")

# ---------------------------------------------------------------- 9 the stopping rule in the prover
rr = os.path.join(DATA, "rule_runs.json")
if os.path.exists(rr):
    runs = json.load(open(rr))
    r768 = [r for r in runs if r["set"] == "rule768_local"]
    rows = []
    for k, r in enumerate(r768):
        rows.append("<tr><td>%d</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%.0f</td><td>%s</td><td>%s</td><td>%.1f</td><td>%.0f</td></tr>"
                    % (k, fmt(r["first_chunk_D"]), fmt(r["orders"]), fmt(r["scan_s"]), fmt(r["ecm_s"]), r["search_wall"],
                       "%.2f" % max(r["expwin"], 0.0) if r["expwin"] is not None else "&ndash;", fmt(r["D"]).replace("-", "&minus;"),
                       r["build_s"], r["chain_s"]))
    R["RULE_TABLE"] = ("<div class=\"wide\"><table><tr><th>run</th><th>solvable D below 6&middot;10<sup>7</sup></th><th>orders bought</th>"
                       "<th>scan, thread-s</th><th>ECM, thread-s</th><th>search, s</th><th>wins the model expected</th>"
                       "<th>top-level D</th><th>build, s</th><th>whole chain, s</th></tr>" + "".join(rows) + "</table></div>")
    n = len(r768)
    work = [r["scan_s"] + r["ecm_s"] for r in r768]
    chain = sorted(r["chain_s"] for r in r768)
    ew = [max(r["expwin"], 0.0) for r in r768 if r["expwin"] is not None]
    R["RULE_N"] = n
    R["RULE_CHAIN_MEAN"] = "%.0f" % (sum(chain) / n)
    R["RULE_CHAIN_MEDIAN"] = "%.0f" % chain[n // 2]
    R["RULE_CHAIN_MIN"] = "%.0f" % chain[0]
    R["RULE_CHAIN_MAX"] = "%.0f" % chain[-1]
    R["RULE_WORK_MEAN"] = fmt(sum(work) / n)
    R["RULE_WORK_MEAN_H"] = "%.1f" % (sum(work) / n / 3600)
    R["RULE_WORK_MEDIAN_H"] = "%.1f" % (sorted(work)[n // 2] / 3600)
    R["RULE_SEARCH_MEAN"] = "%.0f" % (sum(r["search_wall"] for r in r768) / n)
    R["RULE_EW_N"] = len(ew)
    R["RULE_EW_MEAN"] = "%.2f" % (sum(ew) / len(ew))
    R["RULE_EW_SE"] = "%.2f" % (1.0 / len(ew) ** 0.5)        # of the mean of that many Exp(1) draws
    R["RULE_EW_SUM"] = "%.1f" % sum(ew)
    R["RULE_D_MIN"] = fmt(min(r["first_chunk_D"] for r in r768))
    R["RULE_D_MAX"] = fmt(max(r["first_chunk_D"] for r in r768))
    big = [r for r in runs if r["set"] == "rule1030_local"]
    if big:
        b = big[0]
        R["RULE1030_ORDERS"] = fmt(b["orders"])
        R["RULE1030_SCAN_H"] = "%.1f" % (b["scan_s"] / 3600)
        R["RULE1030_ECM_H"] = "%.1f" % (b["ecm_s"] / 3600)
        R["RULE1030_WORK_H"] = "%.1f" % ((b["scan_s"] + b["ecm_s"]) / 3600)
        R["RULE1030_SEARCH"] = fmt(b["search_wall"])
        R["RULE1030_SEARCH_H"] = "%.2f" % (b["search_wall"] / 3600)
        R["RULE1030_EW"] = "%.2f" % b["expwin"]
        R["RULE1030_D"] = fmt(b["D"]).replace("-", "&minus;")
        R["RULE1030_BUILD"] = "%.0f" % b["build_s"]
        R["RULE1030_CHAIN"] = fmt(b["chain_s"])
        R["RULE1030_CHAIN_H"] = "%.2f" % (b["chain_s"] / 3600)
        R["RULE1030_R_H"] = "%.0f" % (b["R"] / 3600)
        R["RULE1030_FACTOR"] = "%.2f" % b["cost_factor"]
        R["RULE1030_C8A_H"] = "%.0f" % ((b["scan_s"] + b["ecm_s"]) / 3600 / b["cost_factor"])
        cj = os.path.join(DATA, "rule1030_local", "p1030.log")
        mm = __import__("re").search(r"RUN wall=([\d.]+) user=([\d.]+) sys=([\d.]+)", open(cj).read())
        R["RULE1030_CPU_H"] = "%.1f" % ((float(mm.group(2)) + float(mm.group(3))) / 3600)



# ---------------------------------------------------------------- 10 the rule on 192 cores
def run_stats(a):
    """Summary of a list of run records."""
    n = len(a)
    work = sorted(r["scan_s"] + r["ecm_s"] for r in a)
    chain = sorted(r["chain_s"] for r in a)
    ew = [max(r["expwin"], 0.0) for r in a if r["expwin"] is not None]
    return {"n": n, "chain_mean": sum(chain) / n, "chain_med": chain[n // 2], "chain_min": chain[0], "chain_max": chain[-1],
            "wall_mean": sum(r["search_wall"] for r in a) / n, "wall_med": sorted(r["search_wall"] for r in a)[n // 2],
            "work_mean": sum(work) / n, "work_med": work[n // 2], "orders_mean": sum(r["orders"] for r in a) / n,
            "scan_mean": sum(r["scan_s"] for r in a) / n, "ecm_mean": sum(r["ecm_s"] for r in a) / n,
            "ew_mean": sum(ew) / max(len(ew), 1), "ew_n": len(ew), "ew": ew,
            "build_mean": sum(r["build_s"] for r in a) / n, "build_max": max(r["build_s"] for r in a),
            "cpu_mean": sum(r["cpu_s"] or 0 for r in a) / n, "ok": all(r["verified"] for r in a),
            "machine_mean": sum(r["search_wall"] * (r["threads"] or 0) for r in a) / n}


def run_group(runs, sname, prefix):
    return [r for r in runs if r["set"] == sname and os.path.basename(r["file"]).startswith(prefix)]


def run_rows(a, unit="core"):
    rows = []
    for r in a:
        k = int(re.search(r"_(\d+)\.log$", r["file"]).group(1))
        rows.append("<tr><td>%d</td><td>%s</td><td>%s</td><td>%s</td><td>%.1f</td><td>%s</td><td>%s</td><td>%.1f</td><td>%.1f</td><td>%.2f</td></tr>"
                    % (k, fmt(r["orders"]), fmt(r["scan_s"]), fmt(r["ecm_s"]), r["search_wall"],
                       "%.2f" % max(r["expwin"], 0.0) if r["expwin"] is not None else "&ndash;", fmt(r["D"]).replace("-", "&minus;"),
                       r["build_s"], r["chain_s"], (r["cpu_s"] or 0) / 3600))
    return ("<div class=\"wide\"><table><tr><th>prime</th><th>orders bought</th><th>scan, %s-s</th><th>ECM, %s-s</th><th>search, s</th>"
            "<th>wins the model expected</th><th>top-level D</th><th>build, s</th><th>whole chain, s</th><th>CPU, %s-h</th></tr>"
            % (unit, unit, unit) + "".join(rows) + "</table></div>")


if os.path.exists(rr):
    ga, gb, gc = (run_group(runs, "rule_c8a", p) for p in ("a768_", "b1030_", "c768_"))
    if ga:
        sa, sl = run_stats(ga), run_stats(r768)
        g768 = ga + gc
        s2 = run_stats(g768)
        R["C8A_TABLE_768"] = run_rows(g768)
        R["C8A_N768"] = s2["n"]
        cmp_rows = []
        cols = [("16-core box, 32 threads (9.3)", sl), ("192 cores, the same 13 primes", sa)]
        if gc:
            cols.append(("192 cores, 13 more primes", run_stats(gc)))
            cols.append(("192 cores, all %d" % s2["n"], s2))
        def line(name, f):
            return "<tr><td>%s</td>%s</tr>" % (name, "".join("<td>%s</td>" % f(c[1]) for c in cols))
        R["C8A_CMP"] = ("<table><tr><th>768 bits</th>" + "".join("<th>%s</th>" % c[0] for c in cols) + "</tr>"
                        + line("whole chain, s: mean (median)", lambda x: "<b>%.0f</b> (%.0f)" % (x["chain_mean"], x["chain_med"]))
                        + line("slowest chain, s", lambda x: "%.0f" % x["chain_max"])
                        + line("top-level search, s: mean (median)", lambda x: "%.0f (%.0f)" % (x["wall_mean"], x["wall_med"]))
                        + line("search work, core-h: mean (median)", lambda x: "%.2f (%.2f)" % (x["work_mean"] / 3600, x["work_med"] / 3600))
                        + line("of which scan", lambda x: "%.2f" % (x["scan_mean"] / 3600))
                        + line("orders bought, mean", lambda x: fmt(x["orders_mean"]))
                        + line("build of the top-level curve, s: mean (max)", lambda x: "%.1f (%.0f)" % (x["build_mean"], x["build_max"]))
                        + line("wins the model expected at the win, mean", lambda x: "%.2f" % x["ew_mean"])
                        + "</table>")
        for k, v in (("A", sa), ("ALL", s2)):
            R["C8A_%s_CHAIN_MEAN" % k] = "%.0f" % v["chain_mean"]
            R["C8A_%s_CHAIN_MED" % k] = "%.0f" % v["chain_med"]
            R["C8A_%s_CHAIN_MIN" % k] = "%.0f" % v["chain_min"]
            R["C8A_%s_CHAIN_MAX" % k] = "%.0f" % v["chain_max"]
            R["C8A_%s_WALL_MEAN" % k] = "%.0f" % v["wall_mean"]
            R["C8A_%s_WORK_H" % k] = "%.2f" % (v["work_mean"] / 3600)
            R["C8A_%s_MACHINE_H" % k] = "%.2f" % (v["machine_mean"] / 3600)
            R["C8A_%s_UTIL" % k] = "%.0f" % (100 * v["work_mean"] / v["machine_mean"])
            R["C8A_%s_CPU_H" % k] = "%.2f" % (v["cpu_mean"] / 3600)
            R["C8A_%s_EW" % k] = "%.2f" % v["ew_mean"]
        R["C8A_LOCAL_WORK_H"] = "%.2f" % (sl["work_mean"] / 3600)
        R["C8A_ALL_BUILD_MEAN"] = "%.1f" % s2["build_mean"]
        R["C8A_ALL_BUILD_MAX"] = "%.0f" % s2["build_max"]
        if gc:
            sc = run_stats(gc)
            R["C8A_C_WORK_H"] = "%.2f" % (sc["work_mean"] / 3600)
            R["C8A_C_EW"] = "%.2f" % sc["ew_mean"]
        # the wins the model expected at the first win, all real 768-bit runs: Exp(1) if the model is right
        ew = sorted(sl["ew"] + s2["ew"])
        n = len(ew)
        ks = max(max(abs((i + 1) / n - (1 - math.exp(-x))), abs(i / n - (1 - math.exp(-x)))) for i, x in enumerate(ew))
        R["EW_ALL_N"] = n
        R["EW_ALL_MEAN"] = "%.2f" % (sum(ew) / n)
        R["EW_ALL_SE"] = "%.2f" % (1 / n ** 0.5)
        R["EW_ALL_KS"] = "%.2f" % ks
        R["EW_ALL_KS5"] = "%.2f" % (1.36 / n ** 0.5)
        R["EW_ALL_LOW"] = sum(1 for x in ew if x < 0.12)
        R["EW_ALL_LOW_EXP"] = "%.1f" % (n * (1 - math.exp(-0.12)))
        # one value per prime: the runs on this box
        ew = sorted(s2["ew"])
        n = len(ew)
        ks = max(max(abs((i + 1) / n - (1 - math.exp(-x))), abs(i / n - (1 - math.exp(-x)))) for i, x in enumerate(ew))
        R["EW_C8A_MEAN"] = "%.2f" % (sum(ew) / n)
        R["EW_C8A_SE"] = "%.2f" % (1 / n ** 0.5)
        R["EW_C8A_KS"] = "%.2f" % ks
        R["EW_C8A_KS5"] = "%.2f" % (1.36 / n ** 0.5)
        R["EW_LOW_N"] = sum(1 for x in ew if x < 0.12)
        R["EW_LOW_EXP"] = "%.1f" % (n * (1 - math.exp(-0.12)))
        # the same primes on the driver box (gmp-ecm, yafu) and here: which discriminant won
        def pidx(r):
            return int(re.search(r"(\d+)(?:_first_run)?\.log$", r["file"]).group(1))
        by = {}
        for r in r768 + ga + run_group(runs, "yafu_local", "y768_"):
            by.setdefault(pidx(r), []).append(r["D"])
        two = sum(1 for k, v in by.items() if len(v) >= 2 and max(v.count(d) for d in v) >= 2)
        three = sum(1 for k, v in by.items() if len(v) >= 3 and max(v.count(d) for d in v) >= 3)
        R["SAME_D_2"] = two
        R["SAME_D_3"] = three
    if gb:
        R["C8A_TABLE_1030"] = run_rows(gb)
        sb = run_stats(gb)
        R["C8A_B_N"] = sb["n"]
        R["C8A_B_CHAIN_MEAN"] = fmt(sb["chain_mean"])
        R["C8A_B_WORK_H"] = "%.1f" % (sb["work_mean"] / 3600)
        R["C8A_B_CPU_H"] = "%.1f" % (sb["cpu_mean"] / 3600)
        R["C8A_B_WALL_MEAN"] = fmt(sb["wall_mean"])
        R["C8A_B_EW"] = "%.2f" % sb["ew_mean"]
        R["C8A_B_BUILD_MAX"] = fmt(sb["build_max"])
        R["C8A_B0_CHAIN"] = "%.0f" % gb[0]["chain_s"]
        R["C8A_B0_CPU_H"] = "%.1f" % ((gb[0]["cpu_s"] or 0) / 3600)
        R["C8A_B0_EW"] = "%.2f" % gb[0]["expwin"]
        R["C8A_B0_SEARCH"] = "%.0f" % gb[0]["search_wall"]
        ratios = {}
        for r in gb:
            for b1, v in r["cost_ratio"].items():
                ratios.setdefault(int(b1), []).append(v)
        R["C8A_B_RATIOS"] = ", ".join("%s: %.2f" % (fmt(b1), sum(v) / len(v)) for b1, v in sorted(ratios.items()))
        rmid = [sum(v) / len(v) for b1, v in ratios.items() if 11000 <= b1 <= 3000000]
        R["C8A_B_RATIO_LO"] = "%.2f" % min(rmid)
        R["C8A_B_RATIO_HI"] = "%.2f" % max(rmid)
        R["C8A_B_CHAINS"] = ", ".join("%s" % fmt(r["chain_s"]) for r in gb)
        R["C8A_B_CPUS"] = ", ".join("%.1f" % ((r["cpu_s"] or 0) / 3600) for r in gb)
        R["C8A_B_WORKS"] = ", ".join("%.1f" % ((r["scan_s"] + r["ecm_s"]) / 3600) for r in gb)
        R["C8A_B_BUILDS"] = ", ".join("%.0f" % r["build_s"] for r in gb)
        R["C8A_B_SCAN_SHARE"] = "%.0f" % (100 * sum(r["scan_s"] for r in gb) / sum(r["scan_s"] + r["ecm_s"] for r in gb))
        R["C8A_B_BUILD_SHARE"] = "%.0f" % (100 * sum(r["build_s"] for r in gb) / sum(r["chain_s"] for r in gb))
    call = [r for r in runs if r["set"] == "rule_c8a"]
    if call:
        R["C8A_NCERT"] = len(call)
        R["C8A_NCERT_OK"] = sum(1 for r in call if r["verified"])
        R["C8A_TOTAL_CPU_H"] = "%.0f" % (sum(r["cpu_s"] or 0 for r in call) / 3600)
        R["C8A_TOTAL_WALL_MIN"] = "%.0f" % (sum(r["wall_s"] or 0 for r in call) / 60)
    if gc:
        won = {}
        exp = {}
        for r in gc:
            if r["win_b1"] is not None:
                won[r["win_b1"]] = won.get(r["win_b1"], 0) + 1
            for b1, v in r["expwin_b1"].items():
                exp[int(b1)] = exp.get(int(b1), 0.0) + v
        b1s = [b1 for b1 in sorted(set(won) | set(exp)) if won.get(b1, 0) or exp.get(b1, 0.0) >= 0.05]
        R["C8A_WINB1"] = ("<table><tr><th>B1 of the sweep</th>" + "".join("<th>%s</th>" % fmt(b1) for b1 in b1s) + "<th>all</th></tr>"
                          "<tr><td>first wins, of %d</td>" % len(gc) + "".join("<td>%d</td>" % won.get(b1, 0) for b1 in b1s)
                          + "<td>%d</td></tr>" % sum(won.values())
                          + "<tr><td>wins the model expected from the sweeps run</td>" + "".join("<td>%.1f</td>" % exp.get(b1, 0.0) for b1 in b1s)
                          + "<td>%.1f</td></tr></table>" % sum(exp.values()))

    # build time against the decomposition of the class polynomial, the 29 runs on the c8a
    def lpf(x):
        f, l = 2, 1
        while f * f <= x:
            while x % f == 0:
                l, x = f, x // f
            f += 1
        return max(l, x) if x > 1 else l
    cms = sorted([r for r in call if r.get("cm")], key=lambda r: (r["bits"], lpf(r["cm"]["n"])))
    if cms:
        rows = []
        for r in cms:
            c = r["cm"]
            if lpf(c["n"]) < 150 and r["build_s"] < 8:
                continue
            rows.append("<tr><td>%d</td><td>%s</td><td>%s = %d &times; %s</td><td>%s</td><td>%.1f</td><td>%.1f</td><td>%.1f</td><td>%.1f</td></tr>"
                        % (r["bits"], fmt(r["D"]).replace("-", "&minus;"), fmt(c["deg"]), c["m"], fmt(c["n"]), fmt(lpf(c["n"])),
                           c["poly_s"], c["root1_s"] + c["root2_s"], r["build_s"], r["search_wall"]))
        small = [r for r in cms if lpf(r["cm"]["n"]) < 150 and r["build_s"] < 8]
        R["C8A_BUILD_TABLE"] = ("<table><tr><th>bits</th><th>D</th><th>degree = m &times; n</th><th>largest prime in n</th>"
                                "<th>class polynomial, s</th><th>roots, s</th><th>build, s</th><th>search, s</th></tr>"
                                + "".join(rows) + "</table>")
        R["C8A_BUILD_SMALL_N"] = len(small)
        R["C8A_BUILD_SMALL_MAX"] = "%.1f" % max(r["build_s"] for r in small)
        big = [r for r in cms if lpf(r["cm"]["n"]) >= 1000]
        for nb in (769, 1030):
            v = [r["cm"]["root2_s"] / r["cm"]["n"] * 1000 for r in big if r["bits"] == nb]
            if v:
                R["C8A_BUILD_MS_%d" % nb] = "%.0f to %.0f" % (min(v), max(v))


# ---------------------------------------------------------------- 11 yafu as the ECM engine
ypi = os.path.join(HERE, "..", "..", "ecpp", "ecm_pi_yafu.json")
if os.path.exists(ypi) and HAVE_CALIB:
    yc = {(c["b1"], c["bits"]): c for c in json.load(open(ypi))}
    gc_ = {(c["b1"], c["bits"]): c for c in json.load(open(CALIB)) if c["method"] == "ecm"}
    emy = ss.EcmModel(ypi, b2=ss.B2_YAFU)
    rows = []
    for b1 in sorted({k[0] for k in yc}):
        # three factor sizes per B1: where gmp-ecm's measured success per curve is nearest 10%, 1% and 0.2%
        bits_all = sorted(b for (bb, b) in yc if bb == b1)
        good = [b for b in bits_all if yc[(b1, b)]["hits"] >= 8]
        cells = []
        for target in (0.1, 0.01, 0.002):
            cand = [b for b in good if (b1, b) in gc_ and gc_[(b1, b)]["hits"] >= 8] or good
            if not cand:
                continue
            b = min(cand, key=lambda x: abs(math.log(max(yc[(b1, x)]["hits"] / yc[(b1, x)]["curves"], 1e-9) / target)))
            if b not in [c[0] for c in cells]:
                cells.append((b, yc[(b1, b)], gc_.get((b1, b))))
        for (b, y, g) in sorted(cells):
            py = y["hits"] / y["curves"]
            pg = g["hits"] / g["curves"] if g and g["curves"] and g["hits"] >= 8 else None
            rows.append("<tr><td>%s</td><td>%d</td><td>%s / %s</td><td>%.4f</td><td>%s</td><td>%s</td></tr>"
                        % (fmt(b1), b, fmt(y["hits"]), fmt(y["curves"]), py, "%.4f" % pg if pg else "&ndash;",
                           "%.2f" % (py / pg) if pg else "&ndash;"))
    R["YPI_TABLE"] = ("<table><tr><th>B1</th><th>factor, bits</th><th>yafu: hits / curves</th><th>yafu, per curve</th>"
                      "<th>gmp-ecm, per curve</th><th>ratio</th></tr>" + "".join(rows) + "</table>")
    # overall ratio yafu / gmp-ecm over the cells where both have 8 hits or more, per B1 (weighted by yafu's hits)
    rr_ = []
    for b1 in sorted({k[0] for k in yc}):
        num = den = 0.0
        for (bb, b), y in yc.items():
            g = gc_.get((bb, b))
            if bb == b1 and g and g["hits"] >= 8 and y["hits"] >= 8 and g["hits"] < g["primes"] * 50:
                r_ = (y["hits"] / y["curves"]) / (g["hits"] / g["curves"])
                if y["hits"] / y["curves"] < 0.3:
                    num += math.log(r_) * y["hits"]
                    den += y["hits"]
        if den:
            rr_.append((b1, math.exp(num / den)))
    R["YPI_RATIOS"] = ", ".join("%s: %.2f" % (fmt(b1), v) for b1, v in rr_)
    R["YPI_CURVES"] = fmt(sum(c["curves"] for c in yc.values()))
    R["YPI_CELLS"] = len(yc)
    R["YPI_CPU_H"] = "%.0f" % (sum(c["cpu_s"] for c in yc.values()) / 3600)
    # cost per curve at 1030 bits, c8a reference, and what a 16-curve sweep buys per core-second
    crow = []
    for b1 in sorted(ss.SEC1030_YAFU):
        crow.append("<tr><td>%s</td><td>%.4g</td><td>%.4g</td><td>%.1f</td><td>%s</td><td>%s</td></tr>"
                    % (fmt(b1), ss.SEC1030_LOADED[b1], ss.SEC1030_YAFU[b1], ss.SEC1030_LOADED[b1] / ss.SEC1030_YAFU[b1],
                       fmt(ss.B2_DEFAULT[b1]), fmt(ss.B2_YAFU[b1])))
    R["YCOST_TABLE"] = ("<table><tr><th>B1</th><th>gmp-ecm, core-s per curve</th><th>yafu</th><th>ratio</th>"
                        "<th>B2, gmp-ecm</th><th>B2, yafu</th></tr>" + "".join(crow) + "</table>")
if os.path.exists(rr):
    gy, gy10 = run_group(runs, "yafu_local", "y768_"), run_group(runs, "yafu_local", "y1030_")
    if gy:
        def idx(r):
            m_ = re.search(r"(\d+)(?:_first_run)?\.log$", r["file"])
            return int(m_.group(1))
        gl = {idx(r): r for r in r768}
        rows = []
        for r in gy:
            k = idx(r)
            g = gl.get(k)
            rows.append("<tr><td>%d</td><td>%s</td><td>%s</td><td>%.0f</td><td>%.0f</td><td>%s</td><td>%s</td><td>%.0f</td><td>%s</td><td>%.1f</td><td>%.0f</td></tr>"
                        % (k, fmt(g["orders"]) if g else "&ndash;", fmt(g["scan_s"] + g["ecm_s"]) if g else "&ndash;",
                           g["search_wall"] if g else 0, g["chain_s"] if g else 0,
                           fmt(r["orders"]), fmt(r["scan_s"] + r["ecm_s"]), r["search_wall"],
                           "%.2f" % max(r["expwin"], 0.0) if r["expwin"] is not None else "&ndash;", r["build_s"], r["chain_s"]))
        R["YENG_LOCAL_TABLE"] = ("<div class=\"wide\"><table><tr><th rowspan=\"2\">prime</th><th colspan=\"4\">gmp-ecm (section 9.3)</th>"
                                 "<th colspan=\"6\">yafu</th></tr><tr><th>orders</th><th>search, thread-s</th><th>search, s</th><th>chain, s</th>"
                                 "<th>orders</th><th>search, thread-s</th><th>search, s</th><th>wins the model expected</th><th>build, s</th>"
                                 "<th>chain, s</th></tr>" + "".join(rows) + "</table></div>")
        sy, sg = run_stats(gy), run_stats([gl[idx(r)] for r in gy if idx(r) in gl])
        R["YL_N"] = sy["n"]
        R["YL_OK"] = "all" if sy["ok"] else "NOT ALL"
        R["YL_CHAIN_MEAN"] = "%.0f" % sy["chain_mean"]
        R["YL_CHAIN_MED"] = "%.0f" % sy["chain_med"]
        R["YL_CHAIN_MIN"] = "%.0f" % sy["chain_min"]
        R["YL_CHAIN_MAX"] = "%.0f" % sy["chain_max"]
        R["YL_WORK_H"] = "%.2f" % (sy["work_mean"] / 3600)
        R["YL_WORK_MED_H"] = "%.2f" % (sy["work_med"] / 3600)
        R["YL_WALL_MEAN"] = "%.0f" % sy["wall_mean"]
        R["YL_EW"] = "%.2f" % sy["ew_mean"]
        R["GL_CHAIN_MEAN"] = "%.0f" % sg["chain_mean"]
        R["GL_CHAIN_MED"] = "%.0f" % sg["chain_med"]
        R["GL_WORK_H"] = "%.2f" % (sg["work_mean"] / 3600)
        R["GL_WORK_MED_H"] = "%.2f" % (sg["work_med"] / 3600)
        R["GL_WALL_MEAN"] = "%.0f" % sg["wall_mean"]
        R["YL_SCAN_SHARE"] = "%.0f" % (100 * sy["scan_mean"] / sy["work_mean"])
        R["YL_BUILD_MEAN"] = "%.1f" % sy["build_mean"]
        ratios = {}
        for r in gy + gy10:
            for b1, v in r["cost_ratio"].items():
                ratios.setdefault(int(b1), []).append(v)
        R["YL_RATIOS"] = ", ".join("%s: %.2f" % (fmt(b1), sum(v) / len(v)) for b1, v in sorted(ratios.items()))
        won, exp = {}, {}
        for r in gy:
            if r["win_b1"] is not None:
                won[r["win_b1"]] = won.get(r["win_b1"], 0) + 1
            for b1, v in r["expwin_b1"].items():
                exp[int(b1)] = exp.get(int(b1), 0.0) + v
        R["YL_WON"] = ", ".join("%s: %d" % (fmt(b1), n_) for b1, n_ in sorted(won.items()))
        R["YL_EXP"] = ", ".join("%s: %.1f" % (fmt(b1), v) for b1, v in sorted(exp.items()) if v >= 0.05)
    if gy10:
        b = gy10[0]
        g = big[0] if big else None
        R["Y1030_CHAIN"] = fmt(b["chain_s"])
        R["Y1030_CHAIN_MIN"] = "%.0f" % (b["chain_s"] / 60)
        R["Y1030_SEARCH"] = fmt(b["search_wall"])
        R["Y1030_ORDERS"] = fmt(b["orders"])
        R["Y1030_SCAN_H"] = "%.1f" % (b["scan_s"] / 3600)
        R["Y1030_ECM_H"] = "%.1f" % (b["ecm_s"] / 3600)
        R["Y1030_WORK_H"] = "%.1f" % ((b["scan_s"] + b["ecm_s"]) / 3600)
        R["Y1030_CPU_H"] = "%.1f" % ((b["cpu_s"] or 0) / 3600)
        R["Y1030_EW"] = "%.2f" % b["expwin"]
        R["Y1030_BUILD"] = "%.0f" % b["build_s"]
        R["Y1030_D"] = fmt(b["D"]).replace("-", "&minus;")
        R["Y1030_FACTOR"] = "%.2f" % b["cost_factor"]

yl = os.path.join(DATA, "yafu_cost2_local.json")
if os.path.exists(yl):
    t = {(r["bits"], r["b1"]): r["sec_per_curve"] for r in json.load(open(yl))}
    sizes = sorted({k[0] for k in t})
    b1s = [b for b in sorted({k[1] for k in t}) if b <= 1000000]
    R["YSIZE_TABLE"] = ("<table><tr><th>number, bits</th>" + "".join("<th>B1 = %s</th>" % fmt(b) for b in b1s) + "</tr>"
                        + "".join("<tr><td>%d</td>%s</tr>" % (sz, "".join("<td>%.4f</td>" % t[(sz, b)] for b in b1s)) for sz in sizes)
                        + "</table>")

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
