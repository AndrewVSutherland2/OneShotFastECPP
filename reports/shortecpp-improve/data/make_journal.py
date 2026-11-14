#!/usr/bin/env python3
"""Turn the winners of policy_run.json into a shortECPP winner journal, so the
prover builds the top level from them (smallest |D| first) and finishes the chain.
usage: make_journal.py <expdir> <journal.jsonl>"""
import json, os, sys
E, out = sys.argv[1], sys.argv[2]
p = int(open(os.path.join(E, "p768.txt")).read())
run = json.load(open(os.path.join(E, "policy_run.json")))
want = {w["i"]: w for w in run["winners"]}
recs = []
with open(os.path.join(E, "strip.txt")) as f, open(os.path.join(E, "orders_d.txt")) as g:
    for i, (line, dl) in enumerate(zip(f, g)):
        if i in want:
            N = int(line.split()[0]); d, v = map(int, dl.split()); q = int(want[i]["q"])
            t = abs(p + 1 - N)
            assert t * t + d * v * v == 4 * p and N % q == 0
            recs.append({"p": str(p), "d": d, "t": t, "v": v, "N": str(N), "q": str(q)})
recs.sort(key=lambda r: r["d"])
os.makedirs(os.path.dirname(out), exist_ok=True)
with open(out, "w") as f:
    for r in recs:
        f.write(json.dumps(r) + "\n")
print("%d winners journalled; |D| from %d to %d" % (len(recs), recs[0]["d"], recs[-1]["d"]) if recs else "no winners")
