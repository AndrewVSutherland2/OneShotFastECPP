#!/usr/bin/env python3
"""Sweep: for each (D, inv) list the candidate subgroups (boxes "c=(..)" and lattices "c=(..)/I"), compute the
decomposition over Z for a sample of them, check it against the plain class polynomial with PARI, and record
bound vs actual height.

usage: tower_sweep.py [--lat] [--max=8] D,inv [D,inv ...]      (--lat: lattice candidates only)
Needs classpoly on PATH (setenv.sh) and tower_chk.py next to this file (as chk.py or tower_chk.py)."""
import os, re, subprocess, sys, json, random
S = os.path.dirname(os.path.abspath(__file__))
CHK = os.path.join(S, "chk.py") if os.path.exists(os.path.join(S, "chk.py")) else os.path.join(S, "tower_chk.py")
ENV = dict(os.environ)
def run(D, inv, P, out, env_extra, verb=1):
    e = dict(ENV); e.update(env_extra)
    r = subprocess.run(["classpoly", str(D), str(inv), str(P), out, str(verb)], capture_output=True, text=True, env=e, cwd=S)
    return r.stdout + r.stderr
def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    latonly = "--lat" in sys.argv
    mx = max([int(a[6:]) for a in sys.argv if a.startswith("--max=")] or [8])
    cases = [tuple(map(int, a.split(','))) for a in args]
    res = []
    for D, inv in cases:
        plain = os.path.join(S, "sw_p_%d_%d_%d.txt" % (-D, inv, os.getpid()))
        o = run(D, inv, 0, plain, {})
        if not os.path.exists(plain):
            print("D=%d inv=%d: plain failed: %s" % (D, inv, o[-200:])); continue
        mh = re.search(r"degree (\d+), height (\d+) \((\d+)\)", o)
        ainv = int(open(plain).readline().split('=')[1])
        sq = "square-root" in o
        tw = os.path.join(S, "sw_t_%d_%d_%d.txt" % (-D, inv, os.getpid()))
        o = run(D, ainv, 0, tw, {"CLASSPOLY_TOWER": "1"})
        cands = [(int(a), int(b), cv, float(c)) for a, b, cv, c in re.findall(r"^\s+n=(\d+) m=(\d+) c=\(([\d,]+\)?(?:/\d+)?)\)?: (\d+) bits", o, re.M)]
        cands = [(n, m, cv.replace(")", ""), b) for n, m, cv, b in cands]
        nlat = sum(1 for c in cands if "/" in c[2])
        pres = " ".join(re.findall(r"relative order (\d+)", o)[-8:])
        print("D=%d inv=%d deg=%s plain height %s (bound %s) sqrt=%d cands=%d (%d lattices)" % (D, ainv, mh.group(1), mh.group(2), mh.group(3), sq, len(cands), nlat), flush=True)
        if latonly:
            cands = [c for c in cands if "/" in c[2]]
        random.seed(D)
        pick = cands if len(cands) <= mx else random.sample(cands, mx)
        for n, m, cv, b in pick:
            if os.path.exists(tw): os.remove(tw)
            o = run(D, ainv, 0, tw, {"CLASSPOLY_TOWER_C": cv})
            mt = re.search(r"degree (\d+), height (\d+) \((\d+)\)", o)
            md = re.search(r"Decomposition: F of degree (\d+) \(cosets\), G of degree (\d+)", o)
            if not os.path.exists(tw) or not mt or not md:
                bad = os.path.exists(tw + ".bad")
                rep = "repeated coset value" in o
                print("   n=%d m=%d c=%s: %s %s" % (n, m, cv, "REVERTED (repeated coset value)" if rep else "FAILED (%s)" % ("height margin exceeded" if bad else "no output"), "" if rep else o.strip().splitlines()[-1][:200]))
                res.append({"D": D, "inv": ainv, "n": n, "m": m, "c": cv, "ok": rep, "bad": bad, "reverted": rep})
                if bad: os.remove(tw + ".bad")
                continue
            same = (int(md.group(1)), int(md.group(2))) == (m, n)
            c = subprocess.run([sys.executable, CHK, tw, plain], capture_output=True, text=True).stdout.split()
            _t = re.findall(r"(\d{6,}) ([01])\b", " ".join(c)); ok = same and len(_t) >= 3 and all(x[1] == "1" for x in _t)
            print("   n=%d m=%d c=%s: height %s (bound %s, plain %s) check %s" % (n, m, cv, mt.group(2), mt.group(3), mh.group(2), "OK" if ok else "FAIL " + " ".join(c)), flush=True)
            res.append({"D": D, "inv": ainv, "n": n, "m": m, "c": cv, "ok": ok, "height": int(mt.group(2)), "bound": int(mt.group(3)), "plain": int(mh.group(2)), "plain_bound": int(mh.group(3)), "sqrt": sq})
        for f in (plain, tw):
            if os.path.exists(f): os.remove(f)
    json.dump(res, open(os.path.join(S, "sweep_%d.json" % os.getpid()), "w"))
    nok = sum(1 for r in res if r["ok"])
    print("SWEEP %d/%d OK (%d reverted for a repeated coset value)" % (nok, len(res), sum(1 for r in res if r.get("reverted"))))
main()
