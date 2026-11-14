#!/usr/bin/env python3
"""ECM engine benchmark driver; runs on the spot box after build.sh.

Phases (results accumulate in out/results.json, raw output under out/raw/):
  0  GMP micro-benchmark (system / k8 / zen builds), one pinned core
  A  per-curve stage 1 / stage 2 times, light load (24 pinned jobs at a time)
  B  throughput with every core busy, 1030 bits
  C  launch overhead of one-process-per-call at B1=2000
  D  avx-ecm stage 1 + gmp-ecm resume: does it find factors at the expected rate
"""
import concurrent.futures as cf
import json
import multiprocessing
import os
import queue
import random
import re
import resource
import shutil
import subprocess
import sys
import time

R = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(R, "out")
RAW = os.path.join(OUT, "raw")
os.makedirs(RAW, exist_ok=True)
NCPU = os.cpu_count()
DEBLIB = R + "/inst/deb/usr/lib/x86_64-linux-gnu"
ENG = {"deb705": (R + "/inst/deb/usr/bin/ecm", {"LD_LIBRARY_PATH": DEBLIB}),
       "ecm-k8": (R + "/inst/ecm-k8", {}),
       "ecm-zen": (R + "/inst/ecm-zen", {}),
       "ecm-zentuned": (R + "/inst/ecm-zentuned", {})}
AVX = {"avx-icl": R + "/inst/avx-icl", "avx-native": R + "/inst/avx-native"}
B1S = [2000, 11000, 50000, 250000, 1000000, 3000000]
BITS = [768, 1030, 1263]
RES = {"host": os.uname().nodename, "ncpu": NCPU}
PHASES = sys.argv[1] if len(sys.argv) > 1 else "0ABCD"


def log(*a):
    print(time.strftime("[%H:%M:%S]"), *a, flush=True)


def save():
    with open(os.path.join(OUT, "results.json.tmp"), "w") as f:
        json.dump(RES, f, indent=1)
    os.replace(os.path.join(OUT, "results.json.tmp"), os.path.join(OUT, "results.json"))


def env(extra):
    e = dict(os.environ)
    e.update(extra)
    return e


def N(name):
    return open("%s/inputs/%s.txt" % (R, name)).read().strip()


def med(xs):
    xs = sorted(xs)
    return xs[len(xs) // 2] if xs else None


def ecm_cmd(eng, curves, b1, extra=()):
    binp, e = ENG[eng]
    return [binp, "-v", "-c", str(curves)] + list(extra) + [str(b1)], env(e)


def parse_ecm(out):
    s1 = [int(x) for x in re.findall(r"Step 1 took (\d+)ms", out)]
    s2 = [int(x) for x in re.findall(r"Step 2 took (\d+)ms", out)]
    b2 = re.search(r"B2=(\d+)", out)
    return s1, s2, (int(b2.group(1)) if b2 else None)


# ---------------------------------------------------------------- phase 0
def phase0():
    for g in ("sys", "k8", "zen"):
        b = "%s/inst/mb-%s" % (R, g)
        if os.path.exists(b):
            r = subprocess.run(["taskset", "-c", "4", b], capture_output=True, text=True)
            RES.setdefault("micro", {})[g] = r.stdout
            log("micro", g, "\n" + r.stdout)
    save()


# ---------------------------------------------------------------- phase A
def run_gmp_light(eng, bits, b1, core, s1only=False):
    nc = {2000: 40, 11000: 16, 50000: 8, 250000: 4, 1000000: 2, 3000000: 2}[b1]
    cmd, e = ecm_cmd(eng, nc, b1)
    if s1only:
        cmd.append(str(b1))
    t = time.time()
    r = subprocess.run(["taskset", "-c", str(core)] + cmd, input=N("N%d" % bits) + "\n",
                       capture_output=True, text=True, env=e)
    wall = time.time() - t
    tag = "A_%s_%d_%d%s" % (eng, bits, b1, "_s1" if s1only else "")
    open(os.path.join(RAW, tag + ".txt"), "w").write(r.stdout + r.stderr)
    s1, s2, b2 = parse_ecm(r.stdout)
    return {"eng": eng, "bits": bits, "b1": b1, "curves": nc, "s1_ms": med(s1), "s2_ms": med(s2),
            "s1_mean": sum(s1) / max(1, len(s1)), "s2_mean": sum(s2) / max(1, len(s2)),
            "n": len(s1), "b2": b2, "wall_per_curve": wall / nc, "rc": r.returncode}


def run_avx_light(var, bits, b1, core):
    nc = {2000: 64, 11000: 32}.get(b1, 8)
    d = os.path.join(RAW, "A_%s_%d_%d" % (var, bits, b1))
    shutil.rmtree(d, ignore_errors=True)
    os.makedirs(d)
    t = time.time()
    r = subprocess.run(["taskset", "-c", str(core), AVX[var], N("N%d" % bits), str(nc), str(b1), "1", str(b1)],
                       capture_output=True, text=True, cwd=d)
    wall = time.time() - t
    open(os.path.join(d, "stdout.txt"), "w").write(r.stdout + r.stderr)
    st1 = [float(x) for x in re.findall(r"Stage 1 took ([\d.]+) seconds", r.stdout)]
    init = re.search(r"Initialization took ([\d.]+)", r.stdout)
    mb = re.search(r"MAXBITS = (\d+)", r.stdout)
    res = {"eng": var, "bits": bits, "b1": b1, "curves": nc, "rc": r.returncode,
           "s1_ms": 1000 * sum(st1) / nc if st1 else None, "batches": len(st1),
           "init_s": float(init.group(1)) if init else None, "maxbits": int(mb.group(1)) if mb else None,
           "wall": wall}
    # stage 2 of the saved residues with gmp-ecm (first 8 only: timing, not yield)
    sf = os.path.join(d, "save_b1.txt")
    if os.path.exists(sf):
        lines = open(sf).read().splitlines()[:8]
        open(os.path.join(d, "save8.txt"), "w").write("\n".join(lines) + "\n")
        res["save_lines"] = len(lines)
        binp, e = ENG["ecm-zen"]
        r2 = subprocess.run(["taskset", "-c", str(core), binp, "-v", "-resume", "save8.txt", str(b1)],
                            capture_output=True, text=True, cwd=d, env=env(e), stdin=subprocess.DEVNULL)
        open(os.path.join(d, "resume.txt"), "w").write(r2.stdout + r2.stderr)
        s1, s2, b2 = parse_ecm(r2.stdout)
        res.update({"resume_s2_ms": med(s2), "resume_n": len(s2), "resume_s1_ms": med(s1), "resume_b2": b2,
                    "resume_rc": r2.returncode})
    return res


def phaseA():
    cores = queue.Queue()
    for c in range(0, NCPU, max(1, NCPU // 24)):
        cores.put(c)

    def wrap(fn, *a):
        c = cores.get()
        try:
            return fn(*a, c)
        finally:
            cores.put(c)

    jobs = []
    with cf.ThreadPoolExecutor(max_workers=cores.qsize()) as ex:
        for b1 in sorted(B1S, reverse=True):       # longest first
            for bits in BITS:
                for eng in ENG:
                    jobs.append(ex.submit(wrap, run_gmp_light, eng, bits, b1))
                for var in AVX:
                    if os.path.exists(AVX[var]):
                        jobs.append(ex.submit(wrap, run_avx_light, var, bits, b1))
        out = []
        for j in cf.as_completed(jobs):
            try:
                r = j.result()
            except Exception as exc:               # keep going: one bad engine must not sink the run
                log("phase A job failed:", repr(exc))
                continue
            out.append(r)
            log("A", json.dumps(r))
            RES["A"] = out
            save()


def tcurve(eng, b1, what, bits=1030):
    for r in RES.get("A", []):
        if r["eng"] == eng and r["b1"] == b1 and r["bits"] == bits:
            if eng in AVX:
                return (r.get("s1_ms") or 0) / 1000.0
            s1, s2 = r.get("s1_mean") or 0, r.get("s2_mean") or 0
            return (s1 + (s2 if what == "full" else 0)) / 1000.0
    return None


# ---------------------------------------------------------------- phase B
def loaded(tag, mk):
    """Start one pinned process per core (mk(core) -> (cmd, env, cwd, stdin_text)), wait for all."""
    procs = []
    t = time.time()
    for c in range(NCPU):
        cmd, e, cwd, inp = mk(c)
        p = subprocess.Popen(["taskset", "-c", str(c)] + cmd, stdin=subprocess.PIPE, stdout=subprocess.DEVNULL,
                             stderr=subprocess.DEVNULL, env=e, cwd=cwd, text=True)
        if inp is not None:
            p.stdin.write(inp)
        p.stdin.close()
        procs.append(p)
    tl = time.time() - t
    rcs = [p.wait() for p in procs]
    wall = time.time() - t
    return {"tag": tag, "wall": wall, "launch_s": tl, "bad_rc": sum(1 for x in rcs if x not in (0, 2, 6, 8, 10, 14))}


def phaseB():
    out = []
    target = 40.0
    n1030 = N("N1030") + "\n"
    for b1 in (50000, 250000, 1000000):
        for eng in ENG:
            for what in ("full", "s1"):
                if what == "s1" and eng not in ("deb705", "ecm-zentuned"):
                    continue
                t1 = tcurve(eng, b1, what)
                if not t1:
                    continue
                nc = max(1, round(target / t1))
                cmd, e = ecm_cmd(eng, nc, b1)
                cmd = [x for x in cmd if x != "-v"] + ([str(b1)] if what == "s1" else [])
                r = loaded("%s %s B1=%d" % (eng, what, b1), lambda c: (cmd, e, None, n1030))
                r.update({"eng": eng, "what": what, "b1": b1, "curves": nc,
                          "curve_core_s": r["wall"] / nc, "light_curve_s": t1})
                out.append(r)
                log("B", json.dumps(r))
                RES["B"] = out
                save()
        for var in AVX:
            t1 = tcurve(var, b1, "s1")
            if not t1:
                continue
            nc = max(8, 8 * round(target / t1 / 8))
            base = os.path.join(RAW, "B_%s_%d" % (var, b1))
            shutil.rmtree(base, ignore_errors=True)

            def mk(c, var=var, nc=nc, b1=b1, base=base):
                d = os.path.join(base, str(c))
                os.makedirs(d)
                return ([AVX[var], n1030.strip(), str(nc), str(b1), "1", str(b1)], None, d, None)
            r = loaded("%s s1 B1=%d" % (var, b1), mk)
            r.update({"eng": var, "what": "s1", "b1": b1, "curves": nc,
                      "curve_core_s": r["wall"] / nc, "light_curve_s": t1})
            out.append(r)
            log("B", json.dumps(r))
            RES["B"] = out
            save()
            # stage 2 of those residues under full load with gmp-ecm
            if var == "avx-native" or "avx-native" not in [x["eng"] for x in out]:
                binp, e = ENG["ecm-zentuned"]
                k = min(nc, 16)

                def mk2(c, base=base, k=k, b1=b1):
                    d = os.path.join(base, str(c))
                    lines = open(os.path.join(d, "save_b1.txt")).read().splitlines()[:k]
                    open(os.path.join(d, "savek.txt"), "w").write("\n".join(lines) + "\n")
                    return ([binp, "-resume", "savek.txt", str(b1)], env(e), d, None)
                try:
                    r = loaded("resume-s2 of %s B1=%d" % (var, b1), mk2)
                    r.update({"eng": "ecm-zentuned", "what": "s2-resume", "b1": b1, "curves": k,
                              "curve_core_s": r["wall"] / k})
                    out.append(r)
                    log("B", json.dumps(r))
                    RES["B"] = out
                    save()
                except Exception as exc:
                    log("resume phase failed:", repr(exc))


# ---------------------------------------------------------------- phase C
def _peel(args):
    binp, e, nums, c, b1 = args
    r = subprocess.run([binp, "-q", "-c", str(c), str(b1)], input="\n".join(nums) + "\n",
                       capture_output=True, text=True, env=e)
    return len(r.stdout)


def cpu_children():
    u = resource.getrusage(resource.RUSAGE_CHILDREN)
    return u.ru_utime + u.ru_stime


def phaseC():
    rng = random.Random(12345)
    binp, e0 = ENG["deb705"]
    e = env(e0)
    out = []
    for (b1, c, M) in ((2000, 12, 24000), (11000, 20, 6000)):
        nums = [str(rng.getrandbits(1030) | (1 << 1029) | 1) for _ in range(M)]
        for mode in ("thread-per-call", "thread-chunk8", "procpool-per-call"):
            c0, t = cpu_children(), time.time()
            if mode == "thread-per-call":
                with cf.ThreadPoolExecutor(max_workers=NCPU) as ex:
                    list(ex.map(_peel, [(binp, e, [n], c, b1) for n in nums]))
            elif mode == "thread-chunk8":
                with cf.ThreadPoolExecutor(max_workers=NCPU) as ex:
                    list(ex.map(_peel, [(binp, e, nums[i:i + 8], c, b1) for i in range(0, M, 8)]))
            else:
                with multiprocessing.Pool(NCPU) as pool:
                    list(pool.imap_unordered(_peel, [(binp, e, [n], c, b1) for n in nums], chunksize=4))
                    pool.close()
                    pool.join()
            wall = time.time() - t
            cpu = cpu_children() - c0
            r = {"mode": mode, "b1": b1, "c": c, "M": M, "wall": wall, "child_cpu": cpu,
                 "util": cpu / (wall * NCPU), "calls_per_s": M / wall}
            out.append(r)
            log("C", json.dumps(r))
            RES["C"] = out
            save()


# ---------------------------------------------------------------- phase D
def phaseD():
    out = {}
    b1 = 50000
    for var in AVX:
        if not os.path.exists(AVX[var]):
            continue
        d = os.path.join(RAW, "D_" + var)
        shutil.rmtree(d, ignore_errors=True)
        os.makedirs(d)
        r = subprocess.run([AVX[var], N("F20"), "256", str(b1), "32", str(b1)], capture_output=True, text=True, cwd=d)
        open(os.path.join(d, "stdout.txt"), "w").write(r.stdout + r.stderr)
        res = os.path.join(d, "ecm_results.txt")
        s1hits = len(open(res).read().splitlines()) if os.path.exists(res) else 0
        binp, e = ENG["ecm-zen"]
        r2 = subprocess.run([binp, "-v", "-resume", "save_b1.txt", str(b1)], capture_output=True, text=True,
                            cwd=d, env=env(e), stdin=subprocess.DEVNULL)
        open(os.path.join(d, "resume.txt"), "w").write(r2.stdout + r2.stderr)
        out[var] = {"curves": 256, "stage1_result_lines": s1hits,
                    "resume_step1_hits": len(re.findall(r"Factor found in step 1", r2.stdout)),
                    "resume_step2_hits": len(re.findall(r"Factor found in step 2", r2.stdout)),
                    "resume_curves": len(re.findall(r"Step 2 took|Factor found in step 1", r2.stdout)),
                    "rc": r.returncode, "resume_rc": r2.returncode}
        log("D", var, json.dumps(out[var]))
    # reference: 256 independent gmp-ecm curves on the same number
    binp, e0 = ENG["ecm-zen"]
    e = env(e0)
    f20 = N("F20") + "\n"

    def one(_):
        r = subprocess.run([binp, "-v", "-c", "1", str(b1)], input=f20, capture_output=True, text=True, env=e)
        return (1 if "Factor found in step 1" in r.stdout else 0, 1 if "Factor found in step 2" in r.stdout else 0)
    with cf.ThreadPoolExecutor(max_workers=64) as ex:
        hits = list(ex.map(one, range(256)))
    out["gmp-ecm"] = {"curves": 256, "step1_hits": sum(h[0] for h in hits), "step2_hits": sum(h[1] for h in hits)}
    log("D gmp-ecm", json.dumps(out["gmp-ecm"]))
    RES["D"] = out
    save()


if __name__ == "__main__":
    if os.path.exists(os.path.join(OUT, "results.json")):
        try:
            RES.update(json.load(open(os.path.join(OUT, "results.json"))))
        except Exception:
            pass
    for ph, fn in (("0", phase0), ("A", phaseA), ("D", phaseD), ("B", phaseB), ("C", phaseC)):
        if ph in PHASES:
            log("=== phase", ph)
            try:
                fn()
            except Exception as exc:
                log("phase", ph, "FAILED:", repr(exc))
    log("BENCH_DONE")
