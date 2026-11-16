#!/usr/bin/env python3
"""Summarise the prover runs with the stopping-rule scheduler into rule_runs.json.

One record per run log: the top level's search line (orders bought, scan and ECM core-seconds, wall, the reward R
of the rule at the end, the wins the model expected from the work done), the build, and the whole chain.

Sets (directories of logs):
  rule768_local/   13 primes of 768 bits on the 16-core driver box (32 threads), gmp-ecm
  rule1030_local/  nextprime(10^310) on the same box
  rule_c8a/        c8a.48xlarge, 192 cores, gmp-ecm: a768_k (the same 13 primes), b1030_k (three primes of 1030
                   bits), c768_k (13 more primes of 768 bits)
  yafu_c8a/        the same box with yafu as the ECM engine: y768_k, y1030_k
  yafu_local/      the driver box with yafu: y768_k
"Core-seconds" on the driver box are thread-seconds of its 32 hardware threads."""
import glob
import json
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
SETS = ["rule768_local", "rule1030_local", "rule_c8a", "yafu_c8a_1", "yafu_c8a_2", "yafu_c8a_3", "yafu_local"]
SEARCH = re.compile(r"search: (\d+) orders, scan (\d+) \+ ecm (\d+) core-s(?: \((\d+) above the rule\))?, ([\d.]+)s wall, (\d+) units, "
                    r"(\d+) splits, rule R = (\d+) core-s(?:; the model expected (-?[\d.]+) wins from this work)?"
                    r"(?:; measured/modelled ECM cost by B1: ([^;\n]*))?")


def key(fn):
    """a768_3.log -> ("a768", 3); p12.log -> ("p", 12); p0_first_run.log -> ("p", 0)."""
    b = os.path.basename(fn)
    m = re.match(r"([a-z]+\d*)_(\d+)\.log$", b) or re.match(r"([a-z]+)(\d+)", b)
    return (m.group(1), int(m.group(2))) if m else (b, 0)


def parse(fn, sname):
    txt = open(fn).read()
    mb = re.search(r"proving p0 \((\d+) bits\)", txt)
    if not mb:
        return None
    bits = int(mb.group(1))
    m = SEARCH.search(txt)
    if not m:
        return None
    b = re.search(r"level %d bits: built D=-(\d+) h~(\d+) in ([\d.]+)s" % bits, txt)
    c = re.search(r"chain complete: (\d+) levels, ([\d.]+)s", txt)
    f = re.search(r"(gmp-ecm|yafu) here: ([\d.]+) s per curve at B1=\d+ on \d+ bits with (\d+) processes -> cost factor ([\d.]+)", txt)
    fd = re.search(r"dscan\[0,[^)]*\)[^:]*: \+(\d+) D", txt)
    run = re.search(r"RUN \S+ wall=([\d.]+) user=([\d.]+) sys=([\d.]+)", txt)
    ratios = {}
    for tok in (m.group(10) or "").split():
        k, _, v = tok.partition(":")
        try:
            ratios[str(int(float(k)))] = float(v)
        except ValueError:
            pass
    nsearch = len(SEARCH.findall(txt))
    # the top level's cm_method line: degree of the class polynomial, its decomposition m x n, the times
    cm = re.search(r"cm_method D=-(\d+): H_D\^inv\(deg (\d+)(?: = (\d+) x (\d+))?, inv (\d+)\): ([\d.]+)s\s+roots?: ([\d.]+)s(?: \+ ([\d.]+)s)?",
                   txt[m.start():])
    cmrec = None
    if cm and b and int(cm.group(1)) == int(b.group(1)):
        cmrec = {"deg": int(cm.group(2)), "m": int(cm.group(3) or 1), "n": int(cm.group(4) or cm.group(2)), "inv": int(cm.group(5)),
                 "poly_s": float(cm.group(6)), "root1_s": float(cm.group(7)), "root2_s": float(cm.group(8) or 0)}
    line = txt[m.start():txt.find("\n", m.start())]
    eb = re.search(r"wins expected by B1: ([^;]*); won at B1=(\w+)", line)
    ewb = {}
    if eb:
        for tok in eb.group(1).split():
            k, _, v = tok.partition(":")
            ewb[str(int(float(k)))] = float(v)
    return {"set": sname, "file": os.path.relpath(fn, HERE), "bits": bits, "orders": int(m.group(1)), "scan_s": int(m.group(2)),
            "ecm_s": int(m.group(3)), "over_s": int(m.group(4)) if m.group(4) else None, "search_wall": float(m.group(5)),
            "units": int(m.group(6)), "splits": int(m.group(7)), "R": int(m.group(8)),
            "expwin": float(m.group(9)) if m.group(9) else None, "cost_ratio": ratios,
            "expwin_b1": ewb, "win_b1": (int(eb.group(2)) if eb and eb.group(2).isdigit() else None),
            "D": -int(b.group(1)) if b else None, "build_s": float(b.group(3)) if b else None, "cm": cmrec,
            "levels": int(c.group(1)) if c else None, "chain_s": float(c.group(2)) if c else None,
            "engine": ("yafu" if f.group(1) == "yafu" else "gmp") if f else "gmp",
            "threads": int(f.group(3)) if f else None, "cost_factor": float(f.group(4)) if f else None,
            "first_chunk_D": int(fd.group(1)) if fd else None,
            "cpu_s": float(run.group(2)) + float(run.group(3)) if run else None,
            "wall_s": float(run.group(1)) if run else None, "rule_levels": nsearch,
            "verified": "verified twice" in txt}


def main():
    rows = []
    for sname in SETS:
        for fn in sorted(glob.glob(os.path.join(HERE, sname, "*.log")), key=key):
            r = parse(fn, sname)
            if r:
                rows.append(r)
    json.dump(rows, open(os.path.join(HERE, "rule_runs.json"), "w"), indent=1)
    for r in rows:
        print("%-28s %4d bits %-4s %3s thr %6d orders  scan %6d + ecm %7d core-s  search %7.1fs  E[wins] %s  build %5.1fs  chain %7.1fs  %s"
              % (r["file"], r["bits"], r["engine"], r["threads"] or "?", r["orders"], r["scan_s"], r["ecm_s"], r["search_wall"],
                 "%5.2f" % r["expwin"] if r["expwin"] is not None else "  n/a", r["build_s"] or 0, r["chain_s"] or 0,
                 "verified" if r["verified"] else "NOT VERIFIED"))
    groups = {}
    for r in rows:
        g = key(r["file"])[0]
        groups.setdefault((r["set"], g, r["bits"] // 100), []).append(r)
    for (sname, g, _), a in groups.items():
        n = len(a)
        w = sorted(r["scan_s"] + r["ecm_s"] for r in a)
        ch = sorted(r["chain_s"] or 0 for r in a)
        e = [r["expwin"] for r in a if r["expwin"] is not None]
        print("%s %s (%d bits): %d runs, search work mean %.0f core-s (%.2f core-h) median %.2f core-h, search wall mean %.1f s, "
              "chain mean %.1f s median %.1f s, expected wins at the win: mean %.2f over %d"
              % (sname, g, a[0]["bits"], n, sum(w) / n, sum(w) / n / 3600, w[n // 2] / 3600, sum(r["search_wall"] for r in a) / n,
                 sum(ch) / n, ch[n // 2], sum(e) / max(len(e), 1), len(e)))


if __name__ == "__main__":
    main()
