#!/usr/bin/env python3
"""Cost per curve of yafu's AVX-ECM and of gmp-ecm 7.0.5 across B1 and input size, all cores busy
(one single-threaded process of 8 curves per core, input with no factor in reach).
usage: ycost.py <outdir>"""
import concurrent.futures as cf, json, os, random, re, shutil, subprocess, sys, tempfile, time, resource
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
OUT = sys.argv[1]
YAFU = os.path.expanduser("~/yb/yafu/yafu")
BENCH = os.path.expanduser("~/ecmbench")
ECM = BENCH + "/inst/deb/usr/bin/ecm"
ENV = dict(os.environ); ENV["LD_LIBRARY_PATH"] = BENCH + "/inst/deb/usr/lib/x86_64-linux-gnu"
NCPU = os.cpu_count()

def is_prime(n, rng):
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % q == 0: return n == q
    d, s = n - 1, 0
    while d % 2 == 0: d //= 2; s += 1
    for _ in range(16):
        a = rng.randrange(2, n - 1); x = pow(a, d, n)
        if x in (1, n - 1): continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1: break
        else: return False
    return True

def rand_prime(bits, rng):
    while True:
        n = rng.getrandbits(bits) | (1 << (bits - 1)) | 1
        if is_prime(n, rng): return n

def run(eng, n, b1, mode):
    t = time.time()
    if eng == "gmp":
        r = subprocess.run([ECM, "-c", "8", str(b1)], input="%d\n" % n, capture_output=True, text=True, env=ENV)
        m = re.search(r"Using B1=\d+, B2=(\d+)", r.stdout)
    else:
        d = tempfile.mkdtemp(prefix="y", dir=OUT + "/tmp")
        r = subprocess.run([YAFU, "-B1ecm", str(b1), "-threads", "1", "-avxecm_stg2", mode, "-v"], input="ecm(%d,8)\n" % n, capture_output=True, text=True, cwd=d)
        shutil.rmtree(d, ignore_errors=True)
        m = re.search(r"B2=([0-9]+)", r.stdout)
    return time.time() - t, (m.group(1) if m else "")

def main():
    os.makedirs(OUT + "/tmp", exist_ok=True)
    rng = random.Random(7)
    res = []
    for nbits in (768, 1030, 1330):
        hard = [rand_prime(nbits // 2, rng) * rand_prime(nbits - nbits // 2, rng) for _ in range(8)]
        for b1 in (2000, 11000, 50000, 250000, 1000000, 3000000):
            for eng, mode in (("gmp", ""), ("yafu", "1"), ("yafu", "2")):
                if b1 == 3000000 and nbits != 1030: continue
                with cf.ThreadPoolExecutor(max_workers=NCPU) as ex:
                    v = list(ex.map(lambda i: run(eng, hard[i % 8], b1, mode), range(NCPU)))
                ts = sorted(x[0] for x in v)
                rec = {"bits": nbits, "b1": b1, "eng": eng, "mode": mode, "sec_per_curve": ts[len(ts) // 2] / 8, "b2": v[0][1]}
                res.append(rec); print(rec, flush=True)
                json.dump(res, open(OUT + "/ycost.json", "w"))
    # process start-up: yafu and gmp-ecm on a trivial job
    t = time.time(); subprocess.run([YAFU, "-threads", "1"], input="2+2\n", capture_output=True, text=True, cwd=OUT + "/tmp"); print("yafu start-up %.3fs" % (time.time() - t))
    print("YCOST_DONE", flush=True)

if __name__ == "__main__":
    _s, _h = resource.getrlimit(resource.RLIMIT_NOFILE); resource.setrlimit(resource.RLIMIT_NOFILE, (min(_h, 65536), _h))
    main()
