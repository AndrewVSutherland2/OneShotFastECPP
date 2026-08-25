#!/usr/bin/env python3
"""add_record.py -- fold newly proved short ECPPs into the repository.

usage: add_record.py <cert1.txt> [<cert2.txt> ...]

Each cert file holds one flat short-ECPP sequence (whitespace-separated); a
sibling <cert>.json (as written by OneShotFastECPP's shortECPP.py) supplies
timing metadata when present.  Every certificate is re-verified here with
vsmallECPP.verify before anything is written.  certs.csv gains one line per
certificate (kept sorted by p0); a README <details> block is printed for each
so it can be pasted/committed.
"""

import json
import os
import sys
from math import isqrt

HERE = os.environ.get("SPP_DIR", "")
if not HERE:
    for c in ("/home/claude/ShortPrimalityProofs", os.path.expanduser("~/ShortPrimalityProofs")):
        if any(os.path.exists(os.path.join(c, v))
               for v in ("vshortECPP.py", "vsmallECPP.py")):
            HERE = c
            break
def _find_verifier():
    for name in ("vshortECPP.py", "vsmallECPP.py"):
        path = os.path.join(HERE, name)
        if HERE and os.path.exists(path):
            return path
    sys.exit("add_record.py: cannot find the ShortPrimalityProofs repository "
             "(clone it next to this repo or set SPP_DIR) -- refusing to guess "
             "where certs.csv lives")

import importlib.util
_spec = importlib.util.spec_from_file_location("spp_verifier", _find_verifier())
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)
verify = _mod.verify


def load_cert(path):
    seq = [int(t) for t in open(path).read().split()]
    meta = {}
    jp = path + ".json"
    if os.path.exists(jp):
        meta = json.load(open(jp))
    tp = path.rsplit(".txt", 1)[0] + ".time"     # GNU time -v of the prover run
    if os.path.exists(tp):
        cpu = 0.0
        for line in open(tp, errors="replace"):
            if "User time (seconds):" in line or "System time (seconds):" in line:
                cpu += float(line.split(":")[-1])
        if cpu:
            meta["cpu_seconds"] = cpu
    return seq, meta


def details_block(seq, meta):
    p0 = seq[0]
    digits = len(str(p0))
    c = digits - 1
    off = p0 - 10 ** c
    levels = (len(seq) - 1) // 3
    secs = meta.get("seconds")
    threads = meta.get("threads")
    cpu = meta.get("cpu_seconds") or (secs * threads if secs and threads else None)
    timing = ""
    if cpu:
        ch = cpu / 3600.0
        cost = "~%d core-hours" % (round(ch / 10.0) * 10 if ch >= 100 else round(ch))
        if secs and threads:
            wall = ("%.1f hours" % (secs / 3600.0) if secs >= 5400
                    else "%.0f minutes" % (secs / 60.0))
            cost += " (%s wall time on %s cores)" % (wall, threads)
        timing = cost + "; "
    hdr = ("$p=10^{%d}+%d$,&nbsp; <a href=\"https://math.mit.edu/~drew/\">AVS</a> and Claude Code "
           "(Fable 5) via <a href=\"https://github.com/AndrewVSutherland2/OneShotFastECPP\">"
           "OneShotFastECPP/shortECPP.py</a> (%s%d levels)."
           % (c, off, timing, levels))
    return ("<details>\n<summary>%s</summary>\n\n```\n%s\n```\n</details>"
            % (hdr, " ".join(map(str, seq))))


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(2)
    csv_path = os.path.join(HERE, "certs.csv")
    rows = []
    if os.path.exists(csv_path):
        rows = [[int(t) for t in line.split(",")]
                for line in open(csv_path).read().splitlines() if line.strip()]
    have = {r[0] for r in rows}
    blocks = []
    for path in sys.argv[1:]:
        seq, meta = load_cert(path)
        ok = verify(seq)
        print("%s: p0 has %d digits, %d levels, verify=%s"
              % (path, len(str(seq[0])), (len(seq) - 1) // 3, ok))
        if not ok:
            sys.exit("REFUSING to add an unverified certificate: %s" % path)
        if seq[0] in have:
            print("  (already in certs.csv; replacing)")
            rows = [r for r in rows if r[0] != seq[0]]
        rows.append(seq)
        have.add(seq[0])
        blocks.append(details_block(seq, meta))
    rows.sort(key=lambda r: r[0])
    with open(csv_path, "w") as f:
        for r in rows:
            f.write(",".join(map(str, r)) + "\n")
    print("\ncerts.csv now has %d certificates.\n" % len(rows))
    print("README <details> blocks:\n")
    for b in blocks:
        print(b)
        print()


if __name__ == "__main__":
    main()
