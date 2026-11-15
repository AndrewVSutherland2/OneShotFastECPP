#!/usr/bin/env python3
"""Sweep: for each (D, inv) list the candidate subgroups, compute the decomposition over Z for each (or a sample),
check it against the plain class polynomial with PARI, and record bound vs actual height."""
import os, re, subprocess, sys, json, random
S = os.path.dirname(os.path.abspath(__file__))
ENV = dict(os.environ)
def run(D, inv, P, out, env_extra, verb=1):
    e = dict(ENV); e.update(env_extra)
    r = subprocess.run(["classpoly", str(D), str(inv), str(P), out, str(verb)], capture_output=True, text=True, env=e, cwd=S)
    return r.stdout + r.stderr
def main():
    cases = [tuple(map(int, a.split(','))) for a in sys.argv[1:]]
    res = []
    for D, inv in cases:
        plain = os.path.join(S, "sw_p_%d_%d.txt" % (-D, inv))
        o = run(D, inv, 0, plain, {})
        if not os.path.exists(plain):
            print("D=%d inv=%d: plain failed: %s" % (D, inv, o[-200:])); continue
        mh = re.search(r"degree (\d+), height (\d+) \((\d+)\)", o)
        ainv = int(open(plain).readline().split('=')[1])
        sq = "square-root" in o
        tw = os.path.join(S, "sw_t_%d_%d.txt" % (-D, inv))
        o = run(D, ainv, 0, tw, {"CLASSPOLY_TOWER": "1"})
        cands = [(int(a), int(b), cv, float(c)) for a, b, cv, c in re.findall(r"^\s+n=(\d+) m=(\d+) c=\(([\d,]+)\): (\d+) bits", o, re.M)]
        pres = re.findall(r"^\s+(\d+)\^(\d+) =", o, re.M)
        print("D=%d inv=%d deg=%s plain height %s (bound %s) sqrt=%d cands=%d" % (D, ainv, mh.group(1), mh.group(2), mh.group(3), sq, len(cands)), flush=True)
        random.seed(D)
        pick = cands if len(cands) <= 8 else random.sample(cands, 8)
        for n, m, cv, b in pick:
            if os.path.exists(tw): os.remove(tw)
            o = run(D, ainv, 0, tw, {"CLASSPOLY_TOWER_C": cv})
            mt = re.search(r"degree (\d+), height (\d+) \((\d+)\)", o)
            if not os.path.exists(tw) or not mt:
                bad = os.path.exists(tw + ".bad")
                print("   n=%d m=%d c=%s: FAILED (%s) %s" % (n, m, cv, "height margin exceeded" if bad else "no output", o.strip().splitlines()[-1][:200])); 
                res.append({"D": D, "inv": ainv, "n": n, "m": m, "ok": False, "bad": bad}); 
                if bad: os.remove(tw + ".bad")
                continue
            c = subprocess.run([sys.executable, os.path.join(S, "chk.py"), tw, plain], capture_output=True, text=True).stdout.split()
            import re as _re; _t = _re.findall(r"(\d{6,}) ([01])\b", " ".join(c)); ok = len(_t) >= 3 and all(x[1] == "1" for x in _t)
            print("   n=%d m=%d c=%s: height %s (bound %s, plain %s) check %s" % (n, m, cv, mt.group(2), mt.group(3), mh.group(2), "OK" if ok else "FAIL " + " ".join(c)), flush=True)
            res.append({"D": D, "inv": ainv, "n": n, "m": m, "c": cv, "ok": ok, "height": int(mt.group(2)), "bound": int(mt.group(3)), "plain": int(mh.group(2)), "plain_bound": int(mh.group(3)), "sqrt": sq})
        for f in (plain, tw):
            if os.path.exists(f): os.remove(f)
    json.dump(res, open(os.path.join(S, "sweep_%d.json" % os.getpid()), "w"))
    nok = sum(1 for r in res if r["ok"])
    print("SWEEP %d/%d OK" % (nok, len(res)))
main()
