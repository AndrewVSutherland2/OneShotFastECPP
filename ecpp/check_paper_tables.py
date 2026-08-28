#!/usr/bin/env python3
"""check_paper_tables.py -- pin the paper's short-ECPP results tables to the
tracked certificate table.

Parses every row of tab:shortresults and tab:shortrecords in
reports/ecpp-varieties/ecpp-varieties.tex (rows have the shape
"$10^{c}+off$ & bits & levels & ...") and compares (c, offset, bits, levels)
against certs/short2/certs.csv: all 31 chains must appear exactly once with
matching bit lengths and level counts, and the level total must match.  The
timing columns are measurement data (attested by the run logs), not derivable
from the CSV, so they are not checked here.

Exit 0 = tables agree with the CSV; exit 1 = any drift, printed per row.
"""
import os, re, sys

ROOT = os.path.join(os.path.dirname(__file__), "..")
TEX = os.path.join(ROOT, "reports/ecpp-varieties/ecpp-varieties.tex")
CSV = os.path.join(ROOT, "certs/short2/certs.csv")
TIMCSV = os.path.join(ROOT, "certs/short2/timings.csv")

truth = {}
total_levels = 0
for line in open(CSV):
    if not line.strip():
        continue
    ints = [int(t) for t in line.split(',')]
    p0 = ints[0]
    c = len(str(p0)) - 1
    lev = (len(ints) - 1) // 3
    truth[c] = (p0 - 10 ** c, p0.bit_length(), lev)
    total_levels += lev

ROW = re.compile(r"\$10\^\{(\d+)\}\+(\d+)\$\s*&\s*(\d+)\s*&\s*(\d+)\s*&")
tex = open(TEX).read()
blocks = [b for b in re.findall(r"\\begin\{table\}.*?\\end\{table\}", tex, re.S)
          if "tab:shortresults" in b or "tab:shortrecords" in b]
assert len(blocks) == 2, f"expected the two short-ECPP tables, found {len(blocks)}"
seen = {}
bad = 0
for m in ROW.finditer("\n".join(blocks)):
    c, off, bits, lev = (int(g) for g in m.groups())
    if c in seen:
        print(f"c={c}: appears twice in the tables")
        bad += 1
    seen[c] = (off, bits, lev)
    if c not in truth:
        print(f"c={c}: row has no chain in certs.csv")
        bad += 1
    elif truth[c] != (off, bits, lev):
        print(f"c={c}: table says offset +{off}, {bits} bits, {lev} levels; "
              f"certs.csv says +{truth[c][0]}, {truth[c][1]} bits, {truth[c][2]} levels")
        bad += 1
for c in sorted(set(truth) - set(seen)):
    print(f"c={c}: chain in certs.csv but missing from the tables")
    bad += 1
# the records table's wall/threads columns must match the tracked timing data
tim = {}
for line in open(TIMCSV):
    line = line.strip()
    if not line or line.startswith("#") or line.startswith("c,"):
        continue
    f = line.split(",")
    tim[int(f[0])] = (int(f[3]), int(f[4]))
RREC = re.compile(r"\$10\^\{(\d+)\}\+\d+\$\s*&\s*\d+\s*&\s*\d+\s*&\s*([\d{},]+) s(?:\$[^$]*\$)?\s*&\s*(\d+)")
rec_block = [b for b in blocks if "tab:shortrecords" in b][0]
nrec = 0
for m in RREC.finditer(rec_block):
    c = int(m.group(1))
    wall = int(m.group(2).replace("{,}", "").replace(",", ""))
    thr = int(m.group(3))
    nrec += 1
    if c not in tim:
        print(f"c={c}: records row has no entry in timings.csv"); bad += 1
    elif tim[c] != (wall, thr):
        print(f"c={c}: table wall/threads {wall}/{thr} vs timings.csv {tim[c]}"); bad += 1
if nrec != len(tim):
    print(f"records table has {nrec} timed rows, timings.csv has {len(tim)}"); bad += 1
if not bad:
    print(f"{len(seen)} table rows match certs.csv exactly "
          f"({total_levels} levels in total); {nrec} record walls match timings.csv")
sys.exit(1 if bad else 0)
